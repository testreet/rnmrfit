# Uses the nmrfit package to prepare set of ready to load examples.
library(nmrfit)

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

# Checking fit
d3 <- filter_1d(nmrdata3, 30.7, 33.3, round.up = TRUE)

scaffold <- nmrscaffold_1d(peaks, include.phase = FALSE, d3)
fit3 <- nmrfit_1d(scaffold)

print(fit3)
nmr.chol.apod <- nmrdata3

save(nmr.chol, nmr.chol.bad, nmr.chol.apod, file = 'nmr_chol.rda')

