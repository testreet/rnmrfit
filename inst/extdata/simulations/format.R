# Code used to simulate fit data
library(nmrfit)
library(RcppFaddeeva)


#------------------------------------------------------------------------
# Peak functions

# Lorenz
gen_lorenz <- function(x, a, b) {
  z <- (x - a)/b
  (1 + complex(im = z))/(1 + z^2)
}

#------------------------------------------------------------------------
# Peak generator

# Singlet
gen_singlet <- function(snr = 0, phase = 0, baseline = 0, f = gen_lorenz,
                        acqus = list(), procs = list()) {

  # Peak locations are constant
  n <- 2^11
  x <- seq(0, 1, length.out = n)
  y <- f(x, 0.5, 0.01)

  # Adding baseline
  x.baseline <- runif(1)
  y.baseline <- baseline*max(Re(y))*(x - x.baseline)^2
  y <- y + complex(re = y.baseline, im = y.baseline)

  # Adding phase
  y <- phase_spectrum(y, phase)

  # Adding noise
  s = snr*sd(Re(y))
  noise <- complex(re = rnorm(n, 0, s), im = rnorm(n, 0, s))
  y <- y + noise

  # Building acqus and procs lists
  d.acqus <- list(td = 2*n, sw.h = n, sfo1 = 10) 
  for ( i in names(acqus) ) d.acqus[i] <- acqus[i]

  d.procs <- list(si = n, lb = 0, gb = 0, ssb = 0, tm1 = 0, tm2 = 0)
  for ( i in names(procs) ) d.procs[i] <- procs[i]

  # Turning peak into NMRData
  d <- data.frame(direct.shift = x, intensity = y)
  out <- new("NMRData1D", processed = d, procs = d.procs, acqus = d.acqus)

  # Applying apodization
  out <- set_convolution(out)

  if ( length(unique(out@product)) > 1 ) {
    intensity <- out@processed$intensity
    convolution <- out@convolution
    intensity <- convolve(intensity, convolution, type = 'open')
    n.edge <- (length(convolution) - 1)/2
    index <- (1 + n.edge):(length(x) + n.edge)
    out@processed$intensity <- intensity[index]/max(Re(intensity[index]))
  }

  out
}

# Singlet
d <- gen_singlet(snr = 0.1, phase = 0, baseline = 0.2)#, procs = list(ssb = 1))
s <- nmrscaffold_1d('0.5 s', d, n.knots = 1, include.phase = TRUE)
f <- nmrfit_1d(s, include.convolution = FALSE)
plot(d)

# Doublet + singlet
y <- 0.5*f(x, 0.6, 0.01) + f(x, 0.45, 0.02) + f(x, 0.55, 0.02)

# Doublet + doublet 
y <- 0.5*f(x, 0.54, 0.01) + 0.5*f(x, 0.64, 0.01) + 
     f(x, 0.45, 0.02) + f(x, 0.55, 0.02)
#plot(x, y, type = 'l')


