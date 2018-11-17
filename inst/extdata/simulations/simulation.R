# Generuc function definitions for simulation
library(rnmrfit)
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

# Voigt
gen_voigt <- function(x, a, b, c) {
  z <- (x - a + complex(i = b))/(sqrt(2)*c)
  Faddeeva_w(z)/Faddeeva_w(complex(im = b*100)/(sqrt(2)*c*100))
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

