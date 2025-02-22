% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bruker.R
\name{read_procs_1d}
\alias{read_procs_1d}
\title{Read 1D Bruker acqus parameters}
\usage{
read_procs_1d(path, number = NA, ...)
}
\arguments{
\item{path}{Character string that points to a scan directory or file.}

\item{number}{The processing file number to use if reading from
a directory. Defaults to smallest number in pdata folder.}

\item{...}{Arguments passed into process_jcamp().}
}
\value{
A list made up of nested lists with the processed acqus entries.
}
\description{
Reads the procs parameters from directory or file. Since the procs files 
are in JCAMP format, the actual parsing is just a thin wrapper around
read_jcamp(), process_jcamp(), and flatten_jcamp() that reads Bruker 
scan parameters and puts them in a flat list.
}
