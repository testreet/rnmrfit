# Code used to simulate fit data
library(cowplot)
library(ggplot2)
library(purrr)
library(nmrfit)
library(RcppFaddeeva)

# Setting seed for simulations
set.seed(1111)

#------------------------------------------------------------------------
# Peak functions

# Lorenz
gen_lorenz <- function(x, a, b) {
  z <- (x - a)/b
  (1 + complex(im = z))/(1 + z^2)
}

#------------------------------------------------------------------------
# Converts peaks into formatted nmrdata

gen_nmrdata <- function(x, y, nsr = 0, phase = 0, baseline = 0, 
                        apodization = FALSE, 
                        acqus = list(), procs = list()) {

  # Building acqus and procs lists
  d.acqus <- list(td = 2*n, sw.h = n, sfo1 = 100) 
  for ( i in names(acqus) ) d.acqus[i] <- acqus[i]

  d.procs <- list(si = n, lb = 0, gb = 0, ssb = 0, tm1 = 0, tm2 = 0)
  for ( i in names(procs) ) d.procs[i] <- procs[i]

  # Turning peak into NMRData
  d <- data.frame(direct.shift = x, intensity = y)
  out <- new("NMRData1D", processed = d, procs = d.procs, acqus = d.acqus)

  # Applying apodization
  if ( apodization ) {
    out <- set_convolution(out)

    intensity <- out@processed$intensity
    convolution <- out@convolution
    intensity <- convolve(intensity, convolution, type = 'open')
    n.edge <- (length(convolution) - 1)/2
    index <- (1 + n.edge):(length(x) + n.edge)
    out@processed$intensity <- intensity[index]/max(Re(intensity[index]))
  }

  # Going back to x and y data
  d <- processed(out)
  x <- d$direct.shift
  y <- d$intensity

  # Adding baseline
  x.baseline <- runif(1)
  y.baseline <- baseline*max(Re(y))*(x - x.baseline)^2
  y <- y + complex(re = y.baseline, im = y.baseline)

  # Adding phase
  y <- phase_spectrum(y, (rbinom(1,1,.5)*2-1)*phase)

  # Adding noise
  s = nsr*sd(Re(y))
  noise <- complex(re = rnorm(n, 0, s), im = rnorm(n, 0, s))
  y <- y + noise

  # Turning peak into NMRData
  d <- data.frame(direct.shift = x, intensity = y)
  out <- new("NMRData1D", processed = d, procs = d.procs, acqus = d.acqus)

  out
}

#------------------------------------------------------------------------
# Example lorenz peaks

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

# Doublet + singlet
y_d <- 0.5*f(x, 0.6, 0.01) + f(x, 0.45, 0.02) + f(x, 0.55, 0.02)
peaks <- list(a = '0.51 d 10', b = '0.61 s')

d <- gen_nmrdata(x, y_d, nsr = 0.1, phase = 15, baseline = 0.2,
                 procs = list(ssb = 1), apod = FALSE)

s <- nmrscaffold_1d(peaks, d, baseline.degree = 2, 
                    n.knots = 1, include.phase = TRUE)
fit2 <- nmrfit_1d(s, include.convolution = FALSE)


# Doublet + doublet 
y_dd <- 0.5*f(x, 0.57, 0.01) + 0.5*f(x, 0.67, 0.01) + 
            f(x, 0.45, 0.02) + f(x, 0.55, 0.02)
peaks <- list(a = '0.64 d 10', b = '0.52 d 10')

d <- gen_nmrdata(x, y_dd, nsr = 0.1, phase = 15, baseline = 0.2,
                 procs = list(ssb = 1), apod = FALSE)

s <- nmrscaffold_1d(peaks, d, baseline.degree = 2, 
                    n.knots = 1, include.phase = TRUE)
fit3 <- nmrfit_1d(s, include.convolution = FALSE)

#------------------------------------------------------------------------
# Comprehensive comparison

if ( TRUE ) {
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
    peaks <- list(a = '0.5 d 10', b = '0.6 s')
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

  if ( TRUE ) {
    d <- expand.grid(s = 1:100, nsr = c(0.0316, 0.1, 0.316), 
                     phase = c(10, 20, 30))
    fits_s <- mutate(d, fit = pmap(list(nsr, phase, FALSE), f_s),
                     peaks = 'Singlet', true = pi)
    fits_d <- mutate(d, fit = pmap(list(nsr, phase, FALSE), f_d),
                     peaks = 'Doublet + Singlet', true = 4.5*pi)
    fits_dd <- mutate(d, fit = pmap(list(nsr, phase, FALSE), f_dd),
                      peaks = 'Doublet + Doublet', true = 5*pi)

    fits <- rbind(fits_s, fits_d, fits_dd)
    
    save(fits, file = 'fits.rda')
  }
}

load('fits.rda')

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
                                           'Doublet + Singlet',
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
                                           'Doublet + Singlet',
                                           'Doublet + Doublet')))

#------------------------------------------------------------------------
# Plots

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

#------------------------------------------------------------------------
# Initial guess

# Generating data
n <- 2^8
x <- seq(0, 1, length.out = n)
f <- gen_lorenz
y_s <- f(x, 0.5, 0.03)

# Function for generating sequence of initial guesses from 0 to 1 with
# with a concentration around 0.5
f_seq <- function(n, p) {
  x <- seq(-1, 1, length.out = n)
  ifelse(x < 0, -x^p, x^p)/3 + 0.5
}

if ( TRUE ) {
  f_n <- function(nsr, position) {
    peaks <- sprintf('%f s', position)
    d <- gen_nmrdata(x, y_s, nsr = nsr, phase = 0, baseline = 0,
                     procs = list(ssb = 1), apod = FALSE)
    s <- nmrscaffold_1d(peaks, d, baseline.degree = NULL,
                        include.phase = FALSE)
    fit <<- nmrfit_1d(s, include.convolution = FALSE)
    fit
  }

  f_b <- function(nsr, position) {
    peaks <- sprintf('%f s', position)
    d <- gen_nmrdata(x, y_s, nsr = nsr, phase = 0, baseline = 0.2,
                     procs = list(ssb = 1), apod = FALSE)
    s <- nmrscaffold_1d(peaks, d, baseline.degree = 2,
                        n.knots = 1, include.phase = FALSE)
    fit <<- nmrfit_1d(s, include.convolution = FALSE)
    fit
  }

  f_p <- function(nsr, position) {
    peaks <- sprintf('%f s', position)
    d <- gen_nmrdata(x, y_s, nsr = nsr, phase = 30, baseline = 0,
                     procs = list(ssb = 1), apod = FALSE)
    s <- nmrscaffold_1d(peaks, d, baseline.degree = NULL,
                        n.knots = 1, include.phase = TRUE)
    fit <<- nmrfit_1d(s, include.convolution = FALSE)
    fit
  }

  f_bp <- function(nsr, position) {
    peaks <- sprintf('%f s', position)
    d <- gen_nmrdata(x, y_s, nsr = nsr, phase = 30, baseline = .2,
                     procs = list(ssb = 1), apod = FALSE)
    s <- nmrscaffold_1d(peaks, d, baseline.degree = 2,
                        n.knots = 1, include.phase = TRUE)
    fit <<- nmrfit_1d(s, include.convolution = FALSE)
    fit
  }

  if (TRUE) {
  d <- expand.grid(s = 1:1000, guess = f_seq(20, 2))
  fits_n <- mutate(d, fit = pmap(list(.316, guess), f_n),
                   error = 'No Error')
  fits_b <- mutate(d, fit = pmap(list(.316, guess), f_b),
                   error = 'Baseline')
  fits_p <- mutate(d, fit = pmap(list(.316, guess), f_p),
                   error = 'Phase')
  fits_bp <- mutate(d, fit = pmap(list(.316, guess), f_bp),
                   error = 'Baseline + Phase')


  fits <- rbind(fits_n, fits_b, fits_p, fits_bp)
  
  save(fits, file = 'fits2.rda')
  }
}

load('fits2.rda')

# Positions
stats <- fits %>%
           mutate(position = map(fit, function(x) {x@peaks$position})) %>%
           select(-fit) %>%
           unnest(position) %>%
           filter(abs((position - 0.5))/0.5 < 0.01) %>% 
           group_by(error, guess) %>%
           summarize(percent = n())

# Extracting example lineshape
examples <- fits %>%
              filter(s == 1, error == 'No Error', guess == guess[1]) %>%
              mutate(lineshape = map(fit, function(x) {x@nmrdata@processed})) %>%
              unnest(lineshape) %>%
              mutate(Real = Re(intensity),
                     Imaginary = Im(intensity)) %>%
              gather('component', 'intensity', Real, Imaginary)

#------------------------------------------------------------------------
# Plots

theme_set(theme_bw(18))

f_cols <- scales::seq_gradient_pal('grey80', 'cornflowerblue')
cols <- f_cols(seq(0, 1, length.out = 3))

p.ex <- ggplot(examples, 
               aes(x = direct.shift, y = intensity, colour = component)) +
        geom_line() +
        xlab('Relative chemical shift') +
        ylab('Relative intensity') +
        scale_color_manual('Component', values = c('grey', 'black')) +
        theme(legend.position = 'top')


p.guess <- ggplot(stats, aes(x = guess, y = percent, colour = error)) +
           geom_point() +
           geom_line() +
           ylab('Convergence (%)') + 
           xlab('') +
           scale_colour_brewer(palette='Dark2') +
           theme(legend.position = 'top')

p2 <- plot_grid(p.ex, p.guess, ncol = 1)


