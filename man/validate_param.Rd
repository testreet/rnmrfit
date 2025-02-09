% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{validate_param}
\alias{validate_param}
\title{Validate and extract acqus/procs parameters}
\usage{
validate_param(param.list, required.param, sublist = NULL)
}
\arguments{
\item{param.list}{A list possibly containing parameters.}

\item{required.param}{A vector of required parameters.}

\item{sublist}{Name of sublist that parameters may be stored in}
}
\value{
A list containing only the required parameters
}
\description{
Checks to see whether required acqus or procs parameters are present in 
main list or appriate sublist. Outputs the required parameters if
all are present and issues an error otherwise. This function is mainly
intended for internal use.
}
