% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing.R
\name{corr_group_delay}
\alias{corr_group_delay}
\title{Correct Bruker group delay artefact}
\usage{
corr_group_delay(signal, acqus.list)
}
\arguments{
\item{signal}{A vector of complex fid data.}

\item{acqus.list}{A list of acqus parameters that contains 'grpdly'
or 'dspfvs' and 'decim' entries. This list can be 
generated using read_acqus_1d() or through other means.
These parameters can also be nested within a list
item called 'acqus' if multiple dimensions are read
at once.}
}
\value{
Corrected fid data.
}
\description{
Corrects Bruker group delay artefact based on process of Westler and 
Abildgaard (1996): 
https://ucdb.googlecode.com/hg/application/ProSpectND/html/
        dmx_digital_filters.html
}
\details{
See the following for details about artefact:
http://nmr-analysis.blogspot.ca/2008/02/
       why-arent-bruker-fids-time-corrected.html
}
