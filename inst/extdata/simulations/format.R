# Code used to simulate fit data
library(nmrfit)
library(RcppFaddeeva)

x <- seq(0, 1, length.out = 2^11)

#------------------------------------------------------------------------
# Lorenz

f <- function(x, a, b) {
  z <- (x - a)/b
  (1 + complex(im = z))/(1 + z^2)
}

# Singlet
y <- f(x, 0.5, 0.01)

# Doublet + singlet
y <- 0.5*f(x, 0.6, 0.01) + f(x, 0.45, 0.02) + f(x, 0.55, 0.02)

# Doublet + doublet 
y <- 0.5*f(x, 0.54, 0.01) + 0.5*f(x, 0.64, 0.01) + 
     f(x, 0.45, 0.02) + f(x, 0.55, 0.02)
plot(x, y, type = 'l')


