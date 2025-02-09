% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bruker.R
\name{read_procs_file}
\alias{read_procs_file}
\title{Read Bruker procs parameters}
\usage{
read_procs_file(path, ...)
}
\arguments{
\item{path}{Character string that points to a procs file.}

\item{...}{Arguments passed into process_jcamp().}
}
\value{
A list made up of nested lists with the processed procs entries.
}
\description{
Reads the procs parameters from a single file.  Since the procs files are 
in JCAMP format, the actual parsing is just a thin wrapper around
read_jcamp(), process_jcamp(), and flatten_jcamp() that reads Bruker 
scan parameters and puts them in a flat list.
}
