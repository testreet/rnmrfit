# Definition of a super class structure for NMR data.

#------------------------------------------------------------------------
#' Super class combining NMR data and scan parameters.
#'
#' See NMRData1D or NMRData2D.
#'
#' @slot processed A data.frame containing chemical shift and intensity data.
#' @slot acqus A list of acqus parameters.
#' @slot procs A list of procs parameters.
#'
#' @name NMRData-class
#' @export
NMRData <- setClass("NMRData",
                    slots = c(processed = "data.frame",
                              acqus = "list",
                              procs = "list"))



#========================================================================>
# Basic setter and getter functions
#========================================================================>


#------------------------------------------------------------------------
#' @templateVar slot processed 
#' @template NMRData_access
#' @name processed
#' @export
setGeneric("processed", 
           function(object, ...) standardGeneric("processed"))

#' @rdname processed
#' @export
setMethod("processed", "NMRData", 
          function(object) object@processed)

#' @templateVar slot processed 
#' @template NMRData_replacement
#' @name processed-set
#' @export
setGeneric("processed<-", 
           function(object, value) standardGeneric("processed<-"))

#' @rdname processed-set
#' @export
setReplaceMethod("processed", "NMRData",
                 function(object, value) {
                   object@parameters <- value
                   validObject(object)
                   object 
                 })

#------------------------------------------------------------------------
#' @templateVar slot procs 
#' @template NMRData_access
#' @name procs
#' @export
setGeneric("procs", 
           function(object, ...) standardGeneric("procs"))

#' @rdname procs
#' @export
setMethod("procs", "NMRData", 
          function(object) object@procs)

#' @templateVar slot procs 
#' @template NMRData_replacement
#' @name procs-set
#' @export
setGeneric("procs<-", 
           function(object, value) standardGeneric("procs<-"))

#' @rdname procs-set
#' @export
setReplaceMethod("procs", "NMRData",
                 function(object, value) {
                   object@parameters <- value
                   validObject(object)
                   object 
                 })

#------------------------------------------------------------------------
#' @templateVar slot acqus 
#' @template NMRData_access
#' @name acqus
#' @export
setGeneric("acqus", 
           function(object, ...) standardGeneric("acqus"))

#' @rdname acqus
#' @export
setMethod("acqus", "NMRData", 
          function(object) object@acqus)

#' @templateVar slot acqus 
#' @template NMRData_replacement
#' @name acqus-set
#' @export
setGeneric("acqus<-", 
           function(object, value) standardGeneric("acqus<-"))

#' @rdname acqus-set
#' @export
setReplaceMethod("acqus", "NMRData",
                 function(object, value) {
                   object@parameters <- value
                   validObject(object)
                   object 
                 })



#========================================================================>
# Defining list and data.frame like behaviour
#========================================================================>

#------------------------------------------------------------------------
#' Convert NMRData object to list
#'
#' Outputs all slot values as entries in a list.
#'
#' @param x NMRData object.
#'
#' @export
as.list.NMRData <- function(x) {
  list(processed = x@processed, procs = x@procs, acqus = x@acqus)
}

setMethod("as.list", "NMRData", as.list.NMRData)

#------------------------------------------------------------------------
#' Convert NMRData object to data.frame
#'
#' Outputs the processed data.frame, dropping acqus and procs lists.
#'
#' @param x NMRData object.
#'
#' @export
as.data.frame.NMRData <- function(x) {
  x@processed
}

setMethod("as.data.frame", "NMRData", as.data.frame.NMRData)


