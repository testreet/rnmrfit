library(cowplot)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

source('simulation.R')

#------------------------------------------------------------------------
# Defining lorenz peaks

n <- 2^8
x <- seq(0, 1, length.out = n)
f <- gen_lorenz

# Singlet
y_s <- f(x, 0.5, 0.01)
peaks <- '0.5 s'

d <- gen_nmrdata(x, y_s, nsr = 0.1, phase = 15, baseline = 0.2,
                 procs = list(ssb = 1), apod = FALSE)

s <- nmrscaffold_1d(peaks, d, baseline.degree = 2,
                    n.knots = 1, include.phase = TRUE)
fit1 <- nmrfit_1d(s, include.convolution = FALSE)

# Singlet + Doublet
y_d <- f(x, 0.5, 0.01) + 0.5*f(x, 0.45, 0.02) + 0.5*f(x, 0.55, 0.02)
peaks <- list(a = '0.5 s', b = '0.5 d 10')

d <- gen_nmrdata(x, y_d, nsr = 0.1, phase = 15, baseline = 0.2,
                 procs = list(ssb = 1), apod = FALSE)

s <- nmrscaffold_1d(peaks, d, baseline.degree = 2, 
                    n.knots = 1, include.phase = TRUE)
fit2 <- nmrfit_1d(s, include.convolution = FALSE)


# Doublet + doublet 
y_dd <- 0.5*f(x, 0.57, 0.01) + 0.5*f(x, 0.67, 0.01) + 
            f(x, 0.45, 0.02) + f(x, 0.55, 0.02)
peaks <- list(a = '0.62 d 10', b = '0.5 d 10')

d <- gen_nmrdata(x, y_dd, nsr = 0.1, phase = 15, baseline = 0.2,
                 procs = list(ssb = 1), apod = FALSE)

s <- nmrscaffold_1d(peaks, d, baseline.degree = 2, 
                    n.knots = 1, include.phase = TRUE)
fit3 <- nmrfit_1d(s, include.convolution = FALSE)

#------------------------------------------------------------------------
# Defining functions for data generation

f_s <- function(nsr, phase, apod) {
  peaks <- '0.5 s'
  d <- gen_nmrdata(x, y_s, nsr = nsr, phase, baseline = 0.2,
                   procs = list(ssb = 1), apod = apod)
  s <- nmrscaffold_1d(peaks, d, baseline.degree = 2,
                  n.knots = 1, include.phase = TRUE)
  fit <<- nmrfit_1d(s, include.convolution = FALSE)
  fit
}

f_d <- function(nsr, phase, apod) {
  peaks <- list(a = '0.5 s', b = '0.5 d 10')
  d <- gen_nmrdata(x, y_d, nsr = nsr, phase, baseline = 0.2,
                   procs = list(ssb = 1), apod = apod)
  s <- nmrscaffold_1d(peaks, d, baseline.degree = 2, 
                      n.knots = 1, include.phase = TRUE)
  fit <<- nmrfit_1d(s, include.convolution = FALSE)
  fit
}

f_dd <- function(nsr, phase, apod) {
  peaks <- list(a = '0.62 d 10', b = '0.5 d 10')
  d <- gen_nmrdata(x, y_dd, nsr = nsr, phase, baseline = 0.2,
                   procs = list(ssb = 1), apod = apod)
  s <- nmrscaffold_1d(peaks, d, baseline.degree = 2, 
                      n.knots = 1, include.phase = TRUE)
  fit <<- nmrfit_1d(s, include.convolution = FALSE)
  fit
}

#------------------------------------------------------------------------
# The actual simulation

if ( TRUE ) {
  d <- expand.grid(s = 1:100, nsr = c(0.0316, 0.1, 0.316), 
                   phase = c(10, 20, 30))
  fits_s <- mutate(d, fit = pmap(list(nsr, phase, FALSE), f_s),
                   peaks = 'Singlet', true = pi)
  fits_d <- mutate(d, fit = pmap(list(nsr, phase, FALSE), f_d),
                   peaks = 'Singlet + Doublet', true = 3*pi)
  fits_dd <- mutate(d, fit = pmap(list(nsr, phase, FALSE), f_dd),
                    peaks = 'Doublet + Doublet', true = 5*pi)

  fits <- rbind(fits_s, fits_d, fits_dd)
  
  save(fits, file = 'lorenz_error.rda')
}

load('lorenz_error.rda')

#------------------------------------------------------------------------
# Post-processing

# Areas
areas <- fits %>%
           mutate(area = map(fit, calc_area)) %>%
           select(-fit) %>%
           unnest(area) %>%
           group_by(peaks, true, nsr, phase, s) %>%
           summarize(area = sum(area))

stats <- areas %>%
           mutate(error = (area - true)/true) %>%
           group_by(peaks, nsr, phase) %>%
           summarize(error = median(abs(error))*100, 
                     cv = sd(area)/mean(area)*100) %>%
           ungroup()


# Converting to db
stats <- stats %>%
          mutate(snr =round(10*log10(1/stats$nsr)),
                 peaks = factor(as.character(peaks),
                                levels = c('Singlet',
                                           'Singlet + Doublet',
                                           'Doublet + Doublet')),
                 phase = as.factor(phase))

# Extracting example lineshape
examples <- fits %>%
              filter(s == 1) %>%
              mutate(lineshape = map(fit, function(x) {x@nmrdata@processed})) %>%
              unnest(lineshape) %>%
              mutate(intensity = Re(intensity),
                     peaks = factor(as.character(peaks),
                                levels = c('Singlet',
                                           'Singlet + Doublet',
                                           'Doublet + Doublet')))

#------------------------------------------------------------------------
# Plotting

theme_set(theme_bw(18))

f_cols <- scales::seq_gradient_pal('grey80', 'cornflowerblue')
cols <- f_cols(seq(0, 1, length.out = 3))

ex.1 <- filter(examples, nsr < 0.05, phase == 10)
p.ex.1 <- ggplot(ex.1, aes(x = direct.shift, y = intensity)) +
          geom_line() +
          xlab('') +
          ylab('Relative intensity') +
          facet_wrap(~ peaks, nrow = 1)

ex.2 <- filter(examples, nsr > 0.15, phase == 30)
p.ex.2 <- ggplot(ex.2, aes(x = direct.shift, y = intensity)) +
          geom_line() +
          xlab('Relative chemical shift') +
          ylab('Relative intensity') +
          facet_wrap(~ peaks, nrow = 1) +
          theme(legend.position = 'bottom',
                strip.text.x = element_blank())

p.error <- ggplot(stats, aes(x = snr, y = error, 
                             colour = phase, shape = phase)) +
           geom_point() +
           geom_line() +
           ylab('Median absolute error (%)') + 
           xlab('') +
           scale_colour_manual(values = cols) +
           facet_wrap(~ peaks, nrow = 1) +
           theme(legend.position = 'none')

p.cv <- ggplot(stats, aes(x = snr, y = cv, 
                          colour = phase, shape = phase)) +
           geom_point() +
           geom_line() +
           ylab('Coefficient of variance (%)') + 
           xlab('Signal to noise ratio (dB)') +
           scale_colour_manual('Phase error (°)', values = cols) +
           scale_shape_discrete('Phase error (°)') +
           facet_wrap(~ peaks, nrow = 1) +
           theme(legend.position = 'bottom',
                 strip.text.x = element_blank())

p1 <- plot_grid(p.ex.1, p.ex.2, p.error, p.cv, ncol = 1)
ggsave('lorenz_error.pdf', width = 12, height = 12, units = 'in')

