% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMRData1D.R
\docType{methods}
\name{filter_1d}
\alias{filter_1d}
\alias{filter_1d,NMRData1D-method}
\title{Filter 1D chemical shift data by selecting specific regions}
\usage{
filter_1d(object, lower = NULL, upper = NULL, round.up = FALSE, ...)

\S4method{filter_1d}{NMRData1D}(object, lower, upper, round.up)
}
\arguments{
\item{object}{An NMRData1D object.}

\item{lower}{A lower bound for chemical shift (in the direct dimension).}

\item{upper}{An upper bound for chemical shift (in the direct dimension).}

\item{round.up}{True to round up the total number of points to a power of 2.
This is useful to ensure rapid convolution using the fft when
convolution is included in the fit.}
}
\value{
An NMRData1D object with filtered processed data.
}
\description{
Filters processed data to include only that which is contained between
a set of lower and upper bounds on chemical shift.
}
