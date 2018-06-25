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
y_s <- f(x, 0.5, 0.03)

#------------------------------------------------------------------------
# Defining functions for data generation

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

#------------------------------------------------------------------------
# The actual simulation

if (FALSE) {

  # Function for generating sequence of initial guesses from 0 to 1 with
  # with a concentration around 0.5
  f_seq <- function(n, p) {
  x <- seq(-1, 1, length.out = n)
  ifelse(x < 0, -x^p, x^p)/3 + 0.5
  }

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
  
  save(fits, file = 'lorenz_guess.rda')
}

#------------------------------------------------------------------------
# Post-processing

load('lorenz_guess.rda')

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
# Plotting

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
