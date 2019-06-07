# Definition of a super class structure for NMR data.



#------------------------------------------------------------------------------
#' Super class combining NMR data and scan parameters
#' 
#' Fourier transformed data is stored as a tibble (data.frame-like structure)
#' in the processed slot with the acquisition and processing parameters stored
#' as lists in separate slots. If the data is loaded from a JCAMP-DX file
#' without clear "acqus" and "procs" division, the parameters are stored in a
#' generic parameters list slot. Both parameter lists and chemical shift data
#' is clearly divided into "direct" or "indirect" components. If the data is
#' 1D, the "indirect" component is either omitted or NULL.
#' 
#' @slot processed A tibble containing chemical shift and intensity data.
#' @slot parameters A list of generic parameters, divided into "direct" and
#'                  "indirect" dimension components.
#' @slot acqus A list of parameters specific to acqus file, divided into
#'             "direct" and "indirect" dimension components.
#' @slot procs A list of parameters sepcific to procs file, divided into
#'             "direct" and "indirect" dimension components.
#' 
#' @name NMRData-class
#' @export
NMRData <- setClass("NMRData",
  representation(processed = "tbl",
                 parameters = "list",
                 acqus = "list",
                 procs = "list"),
  prototype(parameters = list(direct = list(), indirect = list()),
            acqus = list(direct = list(), indirect = list()),
            procs = list(direct = list(), indirect = list())))



#==============================================================================>
#  Validity check
#==============================================================================>



# A basic check of the components 
validNMRData <- function(object) {

  valid <- TRUE
  err <- c()

  processed <- object@processed
  acqus <- object@acqus
  procs <- object@procs

  # Dirty shortcut to avoid retyping the same thing
  add_err <- function(msg) {
    valid <- FALSE
    err <- c(err, msg)
  }

  #---------------------------------------
  # First the processed data

  # Processed should be a tibble
  msg <- '"processed" should be a tidyverse tibble'
  if (! "tbl" %in% class(processed) ) add_err(msg)

  # Intensity must be an ncmplx or ncmplx2 object
  msg <- '"intensity" columns must be "vctrs_cmplx1" or "vctrs_cmplx2" class'
  logic <- any(class(processed$intensity) %in% c('vctrs_cmplx1','vctrs_cmplx2'))
  if (! logic ) add_err(msg)

  if ( valid ) TRUE
  else err
}

setValidity("NMRData", validNMRData)



#==============================================================================>
#  Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
#' Access specific parameter
#'
#' Looks for named parameter in one of either acqus or procs before defaulting
#' to the generic parameters list.
#'
#' @param NMRData object.
#' @param parameter Name of parameter.
#' @param location One of either "acqus" or "procs".
#' @param dimension One of either "direct" or "indirect".
#' @param error TRUE to issue an error if parameter not found.
#'
#' @name get_parameter
#' @export
setGeneric("get_parameter", 
  function(object, ...) standardGeneric("get_parameter"))

#' @rdname get_parameter
#' @export
setMethod("get_parameter", "NMRData", 
  function(object, parameter, location = "acqus", 
           dimension = "direct", error = FALSE) {

    out <- slot(object, location)[[parameter]]
    if ( is.null(out) ) out <- object@parameters[[parameter]]

    err <- sprintf('Parameter "%s" not found', parameter)
    if ( error && is.null(out) ) stop(err) 

    out
  
  })



#------------------------------------------------------------------------------
#' Access Fourier-transformed data
#' @name processed
#' @export
setGeneric("processed", 
  function(object, ...) standardGeneric("processed"))

#' @rdname processed
#' @export
setMethod("processed", "NMRData", 
  function(object) object@processed)

#' Set Fourier-transformed data
#' @name processed-set
#' @export
setGeneric("processed<-", 
  function(object, value) standardGeneric("processed<-"))

#' @rdname processed-set
#' @export
setReplaceMethod("processed", "NMRData",
  function(object, value) {
    object@processed <- value
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
#' Access generic parameters
#' @name parameters
#' @export
setGeneric("parameters", 
  function(object, ...) standardGeneric("parameters"))

#' @rdname parameters
#' @export
setMethod("parameters", "NMRData", 
  function(object) object@parameters)

#' Set generic parameters
#' @name parameters-set
#' @export
setGeneric("parameters<-", 
  function(object, value) standardGeneric("parameters<-"))

#' @rdname parameters-set
#' @export
setReplaceMethod("parameters", "NMRData",
  function(object, value) {
   object@parameters <- value
   validObject(object)
   object 
  })



#------------------------------------------------------------------------------
#' Access acqus parameters
#' @name acqus
#' @export
setGeneric("acqus", 
  function(object, ...) standardGeneric("acqus"))

#' @rdname acqus
#' @export
setMethod("acqus", "NMRData", 
  function(object) object@acqus)

#' Set acqus parameters
#' @name acqus-set
#' @export
setGeneric("acqus<-", 
  function(object, value) standardGeneric("acqus<-"))

#' @rdname acqus-set
#' @export
setReplaceMethod("acqus", "NMRData",
  function(object, value) {
   object@acqus <- value
   validObject(object)
   object 
  })



#------------------------------------------------------------------------------
#' Access procs parameters
#' @name procs
#' @export
setGeneric("procs", 
  function(object, ...) standardGeneric("procs"))

#' @rdname procs
#' @export
setMethod("procs", "NMRData", 
  function(object) object@procs)

#' Set procs parameters
#' @name procs-set
#' @export
setGeneric("procs<-", 
  function(object, value) standardGeneric("procs<-"))

#' @rdname procs-set
#' @export
setReplaceMethod("procs", "NMRData",
  function(object, value) {
   object@procs <- value
   validObject(object)
   object 
  })



#==============================================================================>
#  Defining list and data.frame like behaviour
#==============================================================================>



#------------------------------------------------------------------------------
#' Convert NMRData object to list
#' @export
as.list.NMRData <- function(x) {
  list(processed = x@processed, parameters = x@parameters, 
       procs = x@procs, acqus = x@acqus)
}

setMethod("as.list", "NMRData", as.list.NMRData)



#------------------------------------------------------------------------------
#' Convert NMRData object to a tidyverse tibble
#' @export
as_tibble.NMRData <- function(x) {
  x@processed
}
