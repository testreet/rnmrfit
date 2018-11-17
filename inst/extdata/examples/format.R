# Uses the nmrfit package to prepare set of ready to load examples.
library(ggplot2)
library(rnmrfit)
library(purrr)

# Peak list
peaks <- list('C2'='31.62 s',
              'C8'='31.90 s',
              'C7'='31.94 s')

# Line-broadening 2 Hz
nmrdata1 <- nmrdata_1d('c13-inept', number = 1)

d1 <- filter_1d(nmrdata1, 30.7, 33.3)

scaffold <- nmrscaffold_1d(peaks, d1)
fit1 <- nmrfit_1d(scaffold)

print(fit1)
nmr.chol <- nmrdata1

# Adding phase and baseline error to the data
d <- processed(nmrdata1)
nmrdata1@processed$intensity <- phase_spectrum(d$intensity, c(15, 15))

d <- processed(nmrdata1)
x <- seq(0, 40*pi, length.out = nrow(d))
y <- 0.05*max(Re(d$intensity))
nmrdata1@processed$intensity <- d$intensity + sin(x)*y

d2 <- filter_1d(nmrdata1, 30.7, 33.3)

scaffold <- nmrscaffold_1d(peaks, d2)
fit2 <- nmrfit_1d(scaffold)

print(fit2)
nmr.chol.bad <- nmrdata1

# Sine apodized 
nmrdata3 <- nmrdata_1d('c13-inept', number = 2)
nmrdata3 <- set_convolution(nmrdata3)

d3 <- filter_1d(nmrdata3, 30.7, 33.3, round.up = TRUE)

scaffold <- nmrscaffold_1d(peaks, include.phase = FALSE, d3)
fit3 <- nmrfit_1d(scaffold)

print(fit3)
nmr.chol.apod <- nmrdata3

save(nmr.chol, nmr.chol.bad, nmr.chol.apod, file = 'nmr_chol.rda')

#------------------------------------------------------------------------
# Static plot of the fits
theme_set(theme_bw(18))

# Left column -- raw
types <- c('Normal', 'Baseline + Phase Error', 'Sine Apodization')
d.fit <- data_frame(type = factor(types, levels = types),
                    fit = list(fit1, fit2, fit3))

d.raw <- d.fit %>%
           mutate(lineshape = map(fit, function(x) {x@nmrdata@processed})) %>%
           unnest(lineshape) %>%
           mutate(Real = Re(intensity),
                  Imaginary = Im(intensity)) %>%
           gather('component', 'intensity', Real, Imaginary) %>%
           filter(direct.shift > 31.5, direct.shift < 32.1)

p.raw <- ggplot(d.raw, 
                aes(x = direct.shift, y = intensity, colour = component)) +
         geom_line() +
         facet_wrap(~ type, ncol = 1) +
         xlab('Chemical shift (ppm)') +
         ylab('Relative intensity') +
         scale_color_manual('Component', values = c('grey', 'black')) +
         theme(legend.position = 'top')

# Right column -- fit
types <- c('Normal', 'Baseline + Phase Error', 'Sine Apodization')
d.fit <- data_frame(type = factor(types, levels = types),
                    fit = list(fit1, fit2, fit3))

d.fit <- d.fit %>%
           mutate(lineshape = map(fit, calc_lineshape)) %>%
           unnest(lineshape) %>%
           group_by(type, direct.shift) %>%
           summarize(intensity = sum(intensity)) %>%
           ungroup() %>%
           mutate(Real = Re(intensity),
                  Imaginary = Im(intensity)) %>%
           gather('component', 'intensity', Real, Imaginary) %>%
           filter(direct.shift > 31.5, direct.shift < 32.1)

p.fit <- ggplot(d.fit, 
                aes(x = direct.shift, y = intensity, colour = component)) +
         geom_line() +
         facet_wrap(~ type, ncol = 1) +
         xlab('Chemical shift (ppm)') +
         ylab('Relative intensity') +
         scale_color_manual('Component', values = c('grey', 'black')) +
         theme(legend.position = 'top')

# Combined
d.combined <- rbind(cbind(d.raw, fit = 'Raw'),
                    cbind(d.fit, fit = 'Fit')) %>%
              mutate(fit = factor(as.character(fit),
                                  levels = c('Raw', 'Fit')),
                     component = factor(as.character(component),
                                        levels = c('Real', 'Imaginary')))

p.combined <- ggplot(d.combined, 
                     aes(x = direct.shift, y = intensity, 
                         colour = fit, linetype = fit)) +
              geom_line() +
              facet_grid(type ~ component, scale = 'free_y') +
              xlab('Chemical shift (ppm)') +
              ylab('Relative intensity') +
              scale_color_manual('Data', values = c('grey', 'black')) +
              scale_linetype_discrete('Data') +
              theme(legend.position = 'bottom')

ggsave('example.pdf', width = 12, height = 9, units = 'in')






