# Definition of a super class structure for lineshape definition.



#========================================================================>
# Documentation entries
#========================================================================>



#' Inherited ellipses description
#' @param ... Additional arguments passed to inheriting methods.
#' @name methodEllipse
NULL



#========================================================================>
#  NMRScaffold -- peak description
#========================================================================>


#------------------------------------------------------------------------
#' Super class defining an NMR lineshape. 
#'
#' See NMRScaffold1D or NMRScaffold2D.
#'
#' @slot peaks A data.frame describing a series of singlets, with one row
#'             per peak. Specific peak parameters depend peak_type, but
#'             all peaks are characterized by a position (in ppm), height
#'             (in relative intensity units), and width (in ppm or Hz).
#' @slot baseline A vector of baseline B-spline coefficients, corresponding 
#'                to the B-spline knots.
#' @slot baseline_difference A vector of baseline B-spline coefficients, 
#'                           corresponding to the difference between real
#'                           and imaginary baselines. Although the two
#'                           are expected to be the same under most conditions,
#'                           this equality may break down in some cases.
#' @slot knots A vector of baseline B-spline knots, corresponding 
#'             to the B-spline coefficients.
#' @slot phase A vector of phase terms in radians.
#' @slot parameters A single vector combining peak, baseline, baseline 
#'                  difference, and phase terms. This slot is meant primarily 
#'                  for internal use and should generally be avoided.
#' @slot constraints A data.frame relating the position and width parameters
#'                   of the peaks, effectively combining singlets into
#'                   multiplets.
#' @slot normalized A logical variable indicating whether parameters have
#'                  been normalized with respect to a set of data. Position
#'                  and width are influenced by the chemical shift range of
#'                  the data, while peak height, baseline, and baseline 
#'                  difference are influenced by the maximum intensity of the 
#'                  data.
#' @slot peak_type The mathematical description of a singlet -- 
#'                 one of either 'lorenz', 'gauss', 'pvoigt', or 'voigt'.
#' @slot peak_units The units of peak width -- either 'ppm' or 'hz'. Although
#'                  Hz units are easier to interpret, the peak fitting
#'                  procedure uses ppm.
#' @slot lower_bounds An NMRScaffold object that sets lower feasible bounds
#'                    on all parameters during peak fitting. Can be set to NULL 
#'                    to leave unbounded
#' @slot upper_bounds An NMRScaffold object that sets upper feasible bounds
#'                    on all parameters during peak fitting. Can be set to NULL 
#'                    to leave unbounded
#' @slot nmrdata An optional NMRData object that serves as a reference for
#'               for normalization and peak_unit conversion. However, this 
#'               can also be provided to the individual functions as needed.
#'
#' @name NMRScaffold-class
#' @export
NMRScaffold <- setClass("NMRScaffold",
                        slots = c(peaks = 'data.frame',
                                  parameters = 'numeric',
                                  constraints = 'data.frame',
                                  normalized = 'logical',
                                  peak_type = 'character',
                                  peak_units = 'character',
                                  bounds = 'list',
                                  .sat_peaks = 'logical',
                                  .sat_constraints = 'logical',
                                  .j_constraints = 'logical'),
                        prototype = prototype(normalized = FALSE,
                                              peak_type = 'lorenz',
                                              peak_units = 'hz',
                                              bounds = list(lower = NULL,
                                                            upper = NULL)))



#========================================================================>
# Validation methods
#========================================================================>



#------------------------------------------------------------------------
#' Generic NRMScaffold validity test
#'
#' Combination of validity tests that are common to Scaffold1D and Scaffold2,
#' namely, everything except baseline and phase.
validNMRScaffold <- function(object) {

  valid <- TRUE
  msg <- c()

  peaks <- object@peaks
  constraints <- object@constraints
  peak.type <- object@peak_type

  # Peak type
  valid.types <- c('lorenz', 'gauss', 'pvoigt', 'voigt')
  if (! peak.type %in% valid.types ) {
    valid <- FALSE
    new.msg <- sprintf('"peak_type" must be one of %s', 
                          paste(valid.types, collapse = ', '))
    msg <- c(msg, new.msg)
  } 

  # Constraint ids must correspond to peaks 
  if ( nrow(constraints) > 0 ) {

    constraint.ids <- c(paste(constraints$id1, constraints$peak1),
                        paste(constraints$id2, constraints$peak2))
    peak.ids <- paste(peaks$id, peaks$peak)

    if (! all(constraint.ids %in% peak.ids) ) {
      valid <- FALSE
      msg <- c(msg,'All constraints must reference peaks by id and peak.')
    }
  }

  if (valid) TRUE
  else msg
}


setValidity("NMRScaffold", validNMRScaffold)

#------------------------------------------------------------------------
#' Validation function for baseline component of NMRScaffold
#'
#' Separate baseline validation functions due to the different structures
#' of 1D and 2D data.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @inheritParams methodEllipse
#'
#' @return A vector of integers or characters.
setGeneric("validNMRScaffoldBaseline", 
           function(object, ...) standardGeneric("validNMRScaffoldBaseline"))

#------------------------------------------------------------------------
#' Validation function for phase component of NMRScaffold
#'
#' Separate phase validation functions due to the different structures
#' of 1D and 2D data.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @inheritParams methodEllipse
#'
#' @return A vector of integers or characters.
setGeneric("validNMRScaffoldPhase", 
           function(object, ...) standardGeneric("validNMRScaffoldPhase"))



#========================================================================>
# Helper functions for determining column names (generics only)
#========================================================================>



#------------------------------------------------------------------------
#' Output columns of "peaks" data.frame that correspond to peak parameters
#'
#' Internal helper function used to generalize peak processing.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @param index TRUE to output column indeces, FALSE to output names
#' @param peak.type One of either 'lorenz', 'gauss', 'pvoigt', or 'voigt'
#' @inheritParams methodEllipse
#'
#' @return A vector of integers or characters.
#' @name data_columns
setGeneric(".data_columns", 
           function(object, index = FALSE, peak.type = NULL, ...) {
             standardGeneric(".data_columns")
           })

#------------------------------------------------------------------------
#' Output all columns of "peaks" data.frame
#'
#' Internal helper function used to generalize peak processing.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @param index TRUE to output column indeces, FALSE to output names
#' @param peak.type One of either 'lorenz', 'gauss', 'pvoigt', or 'voigt'
#' @inheritParams methodEllipse
#'
#' @return A vector of integers or characters.
#' @name all_columns
setGeneric(".all_columns", 
           function(object, index = FALSE, peak.type = NULL,  ...) {
             standardGeneric(".all_columns")
           })

#------------------------------------------------------------------------
#' Output peak position related columns of "peaks" data.frame
#'
#' Internal helper function used to generalize peak processing.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @param index TRUE to output column indeces, FALSE to output names
#' @param peak.type One of either 'lorenz', 'gauss', 'pvoigt', or 'voigt'
#' @inheritParams methodEllipse
#'
#' @return A vector of integers or characters.
#' @name position_columns
setGeneric(".position_columns", 
           function(object, index = FALSE, peak.type = NULL, ...) {
             standardGeneric(".position_columns")
           })

#------------------------------------------------------------------------
#' Output peak height related columns of "peaks" data.frame
#'
#' Internal helper function used to generalize peak processing.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @param index TRUE to output column indeces, FALSE to output names
#' @param peak.type One of either 'lorenz', 'gauss', 'pvoigt', or 'voigt'
#' @inheritParams methodEllipse
#'
#' @return A vector of integers or characters.
#' @name height_columns
setGeneric(".height_columns", 
           function(object, index = FALSE, peak.type = NULL , ...) {
             standardGeneric(".height_columns")
           })

#------------------------------------------------------------------------
#' Output peak width related columns of "peaks" data.frame
#'
#' Internal helper function used to generalize peak processing.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @param index TRUE to output column indeces, FALSE to output names
#' @param peak.type One of either 'lorenz', 'gauss', 'pvoigt', or 'voigt'
#' @inheritParams methodEllipse
#'
#' @return A vector of integers or characters.
#' @name width_columns
setGeneric(".width_columns", 
           function(object, index = FALSE, peak.type = NULL, ...) {
             standardGeneric(".width_columns")
           })



#========================================================================>
# Helper functions dealing with flattened parameters (generics only)
#========================================================================>



#------------------------------------------------------------------------
#' Generate a combined vector from individual parameters
#'
#' Combines a series of individual parameters into a single vector based
#' on the Scaffold being fit.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @param ... Individual parameters that depend on whether the object
#'            is 1D or 2D.
#'
#' @return A combined vector of parameters.
#' @name gen_parameters
setGeneric(".gen_parameters", 
           function(object, ...) standardGeneric(".gen_parameters"))

#------------------------------------------------------------------------
#' Copy individual parameters to combined vector
#'
#' Updates the combined vector of parameters based on the individual values
#' (prior to optimization).
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @param return.object TRUE to return an updated vector, TRUE to return
#'                      the combined vector itself.
#' @inheritParams methodEllipse
#'
#' @return Modified NMRScaffold object or numeric vector.
#' @name merge_parameters
setGeneric(".merge_parameters", 
           function(object, return.object = TRUE, ...) {
             standardGeneric(".merge_parameters")
           })

#------------------------------------------------------------------------
#' Copy combined vector to individual parameters
#'
#' Updates individual parameters based on the combined vector values
#' (after optimization).
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @inheritParams methodEllipse
#'
#' @return Modified NMRScaffold object.
#' @name spread_parameters
setGeneric(".spread_parameters", 
           function(object, ...) standardGeneric(".spread_parameters"))



#========================================================================>
# Internal functions relating to bounds
#========================================================================>



#------------------------------------------------------------------------
#' Propagate function to bounds
#'
#' Since lower and upper bounds should be members of the same class as
#' the object they are bounding, function applied on the parent object
#' can be directly applied to the bounds.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @param f Function to be applied.
#' @param ... Extra arguments to f.
#'
#' @return Modified NMRScaffold object.
#' @name propagate_to_bounds
setGeneric(".propagate_to_bounds", 
           function(object, ...) standardGeneric(".propagate_to_bounds"))

#' @rdname propagate_to_bounds
setMethod(".propagate_to_bounds", "NMRScaffold", 
  function(object, f, ...) {

    # Handling bounds if they exist
    lower <- object@bounds$lower
    upper <- object@bounds$upper

    if (! is.null(lower) ) {
      lower = f(lower, ...)
    }
    
    if (! is.null(upper) ) {
      upper = f(upper, ...)
    }

    object@bounds$lower <- lower
    object@bounds$upper <- upper

    object
  })

#------------------------------------------------------------------------
#' Drop bounds if they don't conform to parent object 
#'
#' If changes made to a parent NMRScaffold object make current bounds
#' no longer applicable, they must be dropped.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @inheritParams methodEllipse
#'
#' @return Modified NMRScaffold object.
#' @name drop_bounds
setGeneric(".drop_bounds", 
           function(object, ...) standardGeneric(".drop_bounds"))

#' @rdname drop_bounds
setMethod(".drop_bounds", "NMRScaffold", 
  function(object) {

    lower <- .is_conformant(object, object@bounds$lower)
    upper <- .is_conformant(object, object@bounds$upper)

    string <- c()
    if (! lower ) string <- c(string, 'lower')
    if (! upper ) string <- c(string, 'upper')
    string <- paste(string, collapse = ' and ')

    if (string != '' ) {
      valid <- FALSE
      msg <- sprintf(paste('The shape of %s bounds no longer',
                           'conforms to parameters, dropping.'), 
                     string)

      warning(msg, call. = FALSE)
    }

    if (! lower ) object@bounds$lower <- NULL
    if (! upper ) object@bounds$upper <- NULL

    object
  })

#------------------------------------------------------------------------
#' Check whether NMRScaffold objects conform
#'
#' For their parameters to be conformant, peaks must have the same rows and
#' columns, while baseline, baseline_difference, and phase must all be the  
#' same length.
#'
#' @param object1 NMRScaffold1D or NMRScaffold2D object.
#' @param object2 NMRScaffold1D or NMRScaffold2D object.
#'
#' @return TRUE if the two objects are conformant, FALSE otherwise.
#' @name is_conformant
setGeneric(".is_conformant", 
           function(object1, object2) standardGeneric(".is_conformant"))



#========================================================================>
# Misc internal functions
#========================================================================>



#------------------------------------------------------------------------
#' Check NMRData availability
#'
#' Validates NMRData input and issues an error if NMRData can not be found
#' as either a slot or an argument.
#'
#' @param object NMRScaffold1D or NMRScaffold2D object.
#' @param nmrdata NMRData1D, NMRData2D, or NULL
#' @inheritParams methodEllipse
#'
#' @return NMRData1D or NMRData2D object.
#' @name check_data
setGeneric(".check_data", 
           function(object, nmrdata, ...) standardGeneric(".check_data"))

#' @rdname check_data
setMethod(".check_data", c("NMRScaffold", "ANY"),
  function(object, nmrdata) {

    if ( is.null(nmrdata) ) {
      if ( is.null(object@nmrdata) ) {
        msg <- '"nmrdata" must be provided if the "nmrdata" slot is NULL'
        stop(msg, call. = FALSE)
      } else {
        nmrdata <- object@nmrdata
      }
    }

    if ( inherits(nmrdata, 'NMRData' )) {
      validObject(nmrdata)
    } else {
      msg <- '"nmrdata" must be a valid NMRData object.'
      stop(msg, call. = FALSE)
    }

    nmrdata
  })



#========================================================================>
# Basic setter and getter functions (to be inherited)
#========================================================================>

# Note, baseline generics are defined in NMRScaffold1D and NMRScaffold2D code

#------------------------------------------------------------------------
#' @templateVar slot peaks
#' @template NMRScaffold_access
#' @name peaks
#' @export
setGeneric("peaks", 
           function(object, ...) standardGeneric("peaks"))

#' @rdname peaks
#' @export
setMethod("peaks", "NMRScaffold", 
          function(object) object@peaks)

#' @templateVar slot peaks
#' @template NMRScaffold_replacement
#' @name peaks-set
#' @export
setGeneric("peaks<-", 
           function(object, value) standardGeneric("peaks<-"))

#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRScaffold",
                 function(object, value) {
                   object@peaks <- value
                   object <- .merge_parameters(object)
                   object <- .drop_bounds(object)
                   validObject(object)
                   object 
                 })

#------------------------------------------------------------------------
#' @templateVar slot phase
#' @template NMRScaffold_access
#' @name phase
#' @export
setGeneric("phase", function(object, ...) 
           standardGeneric("phase"))

#' @rdname phase
#' @export
setMethod("phase", "NMRScaffold", 
          function(object) object@phase)

#' @templateVar slot phase
#' @template NMRScaffold_replacement
#' @name phase-set
#' @export
setGeneric("phase<-", 
           function(object, value) standardGeneric("phase<-"))

#' @rdname phase-set
#' @export
setReplaceMethod("phase", "NMRScaffold",
                 function(object, value) {
                   if ( is.null(value) ) value <- numeric(0)
                   object@phase <- value
                   object <- .merge_parameters(object)
                   object <- .drop_bounds(object)
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @templateVar slot constraints
#' @template NMRScaffold_access
#' @name constraints
#' @export
setGeneric("constraints", function(object, ...) 
           standardGeneric("constraints"))

#' @rdname constraints
#' @export
setMethod("constraints", "NMRScaffold", 
          function(object) object@constraints)

#' @templateVar slot constraints
#' @template NMRScaffold_replacement
#' @name constraints-set
#' @export
setGeneric("constraints<-", 
           function(object, value) standardGeneric("constraints<-"))

#' @rdname constraints-set
#' @export
setReplaceMethod("constraints", "NMRScaffold",
  function(object, value) {
    if ( is.null(value) ) {
      value <- data.frame(id1 = character(0), id2 = character(0),
                          peak1 = character(0), peak2 = character(0),
                          difference = numeric(0), ratio = numeric(0))
    }
    object@constraints <- value
    validObject(object)
    object
    })

#------------------------------------------------------------------------
#' @templateVar slot bounds
#' @template NMRScaffold_access
#' @name bounds
#' @export
setGeneric("bounds", function(object, ...) 
           standardGeneric("bounds"))

#' @rdname bounds
#' @export
setMethod("bounds", "NMRScaffold", 
          function(object) object@bounds)

#' @templateVar slot bounds
#' @template NMRScaffold_replacement
#' @name bounds-set
#' @export
setGeneric("bounds<-", 
           function(object, value) standardGeneric("bounds<-"))

#' @rdname bounds-set
#' @export
setReplaceMethod("bounds", "NMRScaffold",
                 function(object, value) {
                   object@bounds <- value
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @templateVar slot lower_bounds
#' @template NMRScaffold_access
#' @name lower_bounds
#' @export
setGeneric("lower_bounds", function(object, ...) 
           standardGeneric("lower_bounds"))

#' @rdname lower_bounds
#' @export
setMethod("lower_bounds", "NMRScaffold", 
          function(object) object@bounds$lower)

#' @templateVar slot lower_bounds
#' @template NMRScaffold_replacement
#' @name lower_bounds-set
#' @export
setGeneric("lower_bounds<-", 
           function(object, value) standardGeneric("lower_bounds<-"))

#' @rdname lower_bounds-set
#' @export
setReplaceMethod("lower_bounds", "NMRScaffold",
                 function(object, value) {
                   object@bounds$lower <- value
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @templateVar slot upper_bounds
#' @template NMRScaffold_access
#' @name upper_bounds
#' @export
setGeneric("upper_bounds", function(object, ...) 
           standardGeneric("upper_bounds"))

#' @rdname upper_bounds
#' @export
setMethod("upper_bounds", "NMRScaffold", 
          function(object) object@bounds$upper)

#' @templateVar slot upper_bounds
#' @template NMRScaffold_replacement
#' @name upper_bounds-set
#' @export
setGeneric("upper_bounds<-", 
           function(object, value) standardGeneric("upper_bounds<-"))

#' @rdname upper_bounds-set
#' @export
setReplaceMethod("upper_bounds", "NMRScaffold",
                 function(object, value) {
                   object@bounds$upper <- value
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @templateVar slot nmrdata
#' @template NMRScaffold_access
#' @name nmrdata
#' @export
setGeneric("nmrdata", function(object, ...) 
           standardGeneric("nmrdata"))

#' @rdname nmrdata
#' @export
setMethod("nmrdata", "NMRScaffold", 
          function(object) object@nmrdata)

#' @templateVar slot nmrdata
#' @template NMRScaffold_replacement
#' @name nmrdata-set
#' @export
setGeneric("nmrdata<-", 
           function(object, value) standardGeneric("nmrdata<-"))

#' @rdname nmrdata-set
#' @export
setReplaceMethod("nmrdata", "NMRScaffold",
                 function(object, value) {
                   object@nmrdata <- value
                   validObject(object)
                   object
                 })



#========================================================================>
# Compound setter and getter functions (requiring internal calculations)
#========================================================================>



#------------------------------------------------------------------------
#' @templateVar slot peak_type
#' @template NMRScaffold_access
#' @name peak_type
#' @export
setGeneric("peak_type", function(object, ...) 
           standardGeneric("peak_type"))

#' @rdname peak_type
#' @export
setMethod("peak_type", "NMRScaffold", 
          function(object) object@peak_type)

#' Replace the "peak_type" slot of an NMRScaffold object 
#'
#' Generic method to replace the peak_type of NMRScaffold1D or NMRScaffold2D
#' object. This is a convenience function that makes some assumptions, 
#' see set_peak_type() for more details.
#'
#' @name peak_type-set
#' @export
setGeneric("peak_type<-", 
           function(object, value) standardGeneric("peak_type<-"))

#' @rdname peak_type-set
#' @export
setReplaceMethod("peak_type", "NMRScaffold",
                 function(object, value) {
                   object <- set_peak_type(object, value)
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' Convert an NMRScaffold1D or NMRScaffold2D object to a different peak.type
#'
#' As the various peak representations have different parameters, Converting
#' between different peak types is not always trivial. This function ensures
#' that conversions are reversible and path independent (switching directly
#' from lorenz to voigt will be the same as switching from lorenz to gauss
#' to voigt). Each conversion preserves the area of the peak and modifies
#' either the width or height. Switching from pvoigt to lorenz or gauss
#' preserves overall peak height, while switching from voigt to lorenz or
#' gauss preserves overall peak width. To go in the reverse direction,
#' it's necessary to specify what fraction of the original height or width
#' should remain attributed to the lorenz line shape.
#'
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @param peak.type Target peak.type -- one of lorenz, gauss, pvoigt, or voigt.
#' @param frac.lorenz If converting from lorenz or gauss to pvoigt or voigt,
#'                    the fraction of the output peak that is expected to be
#'                    lorenz. Ignored for the reverse conversion.
#' @inheritParams methodEllipse
#'
#' @return An NMRScaffold1D or NMRScaffold2D object with the new peak_type and
#'         modified parameters.
#'
#' @name set_peak_type
#' @export
setGeneric("set_peak_type", 
  function(object, peak.type, frac.lorenz = 0.9, ...) {
    standardGeneric("set_peak_type")
  })


#------------------------------------------------------------------------
#' @templateVar slot normalized
#' @template NMRScaffold_access
#' @name normalized
#' @export
setGeneric("normalized", function(object, ...) 
           standardGeneric("normalized"))

#' @rdname normalized
#' @export
setMethod("normalized", "NMRScaffold", 
          function(object) object@normalized)

#------------------------------------------------------------------------
#' Replace the "normalized" slot of an NMRScaffold object 
#'
#' Generic method to replace the normalized slot of NMRScaffold1D or 
#' NMRScaffold2D object. This is a convenience function, see set_normalized() 
#' for more details.
#'
#' @name normalized-set
#' @export
setGeneric("normalized<-", 
           function(object, value) standardGeneric("normalized<-"))

#' @rdname normalized-set
#' @export
setReplaceMethod("normalized", "NMRScaffold",
                 function(object, value) {
                   if ( is.null(object@nmrdata) ) {
                     msg <- '"nmrdata" is not set, use set_normalized()'
                     stop(msg)
                   }
                   object <- set_normalized(object, normalized = value)
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' Normalize or undo normalization of a NMRScaffold1D or NMRScaffold2D 
#'
#' The fit process benefits from scaling x and y data to a limited 0-1
#' range. Naturally, this requires all parameters be normalized to a 
#' set of data. By default, parameters are normalized to the internally
#' stored nmrdata slot, however, a different NMRData1D object can also
#' be provided.
#'
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @param nmrdata An optional NMRData object used for normalization.
#' @param normalized TRUE to normalize, FALSE to undo normalization.
#' @param include.bounds TRUE to perform the same operation on bounds.
#' @inheritParams methodEllipse
#'
#' @return A new NMRScaffold1D or NMRFit1D object with modified parameters.
#'
#' @name set_normalized
#' @export
setGeneric("set_normalized", 
  function(object, nmrdata = NULL, normalized = TRUE,  
           include.bounds = FALSE, ...) {
    standardGeneric("set_normalized")
  })

#------------------------------------------------------------------------
#' @templateVar slot peak_units
#' @template NMRScaffold_access
#' @name peak_units
#' @export
setGeneric("peak_units", function(object, ...) 
           standardGeneric("peak_units"))

#' @rdname peak_units
#' @export
setMethod("peak_units", "NMRScaffold", 
          function(object) object@peak_units)

#------------------------------------------------------------------------
#' Replace the "peak_units" slot of an NMRScaffold object 
#'
#' Generic method to replace the peak_units slot of NMRScaffold1D or 
#' NMRScaffold2D object. This is a convenience function, see set_peak_types() 
#' for more details.
#'
#' @name peak_units-set
#' @export
setGeneric("peak_units<-", 
           function(object, value) standardGeneric("peak_units<-"))

#' @rdname peak_units-set
#' @export
setReplaceMethod("peak_units", "NMRScaffold",
                 function(object, value) {
                   if ( is.null(object@nmrdata) ) {
                     msg <- '"nmrdata" is not set, use set_peak_units()'
                     stop(msg)
                   }
                   object <- set_peak_units(object, peak.units = value)
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' Switch NMRScaffold1D or NMRScaffold2D parameter units between Hz and ppm
#'
#' Although it's more practical to refer to peak width in terms of Hz, 
#' optimization requires that the width parameters be in ppm. This conversion
#' should generally be handled automatically, but it can be overriden manually
#' if desired.
#'
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @param peak.units One of either 'hz' or 'ppm'.
#' @param nmrdata An optional NMRData object needed for conversion parameters. 
#' @param include.bounds TRUE to perform the same operation on bounds.
#' @inheritParams methodEllipse
#'
#' @return An updated NMRScaffold1D or NMRScaffold2D object.
#'
#' @name set_peak_units
#' @export
setGeneric("set_peak_units", 
  function(object, nmrdata = NULL, peak.units = 'hz', 
           include.bounds = TRUE, ...) {
    standardGeneric("set_peak_units")
  })

#' @rdname set_peak_units
#' @export
setMethod("set_peak_units", "NMRScaffold",
  function(object, nmrdata = NULL, peak.units = 'hz', 
           include.bounds = TRUE) {

  # Ensuring that the nmrdata slot isn't misused quietly
  if (! class(nmrdata) %in% c('NMRData1D', 'NULL') ) {
    msg <- '"nmrdata" must NULL or class "NMRData1D'
    stop(msg)
  }

  # If target peak.units value matches current value, return
  if ( object@peak_units == peak.units ) return(object)

  valid.units <- c('hz', 'ppm')
  if (! peak.units %in% valid.units ) {
    msg <- sprintf('"peak_units" must be one of %s', 
                   paste(valid.units, collapse = ', '))
    stop(msg)
  }

  # Checking nmrdata 
  nmrdata <- .check_data(object, nmrdata)

  # Getting values 
  peaks <- object@peaks
  sfo1 <- nmrdata@acqus[['sfo1']]
  
  w.columns <- .width_columns(object, TRUE)

  if ( peak.units == 'ppm' ) {
    peaks[, w.columns] <- peaks[, w.columns]/sfo1
  } else {
    peaks[, w.columns] <- peaks[, w.columns]*sfo1
  }

  object@peaks <- peaks
  object@peak_units <- peak.units

  # Update parameter slot
  object <- .merge_parameters(object)

  # Updating bounds if necessary
  if ( include.bounds ) {
    .propagate_to_bounds(object, set_peak_units, nmrdata, peak.units)
  } else {
    object
  }
})



#========================================================================>
# Calculations based on current parameter values.
#========================================================================>



#------------------------------------------------------------------------
#' Generate lineshape function
#'
#' This is primarily an internal method that generates a lineshape
#' function that depends on the object input. 
#'
#' TO DO.
#'
#' Warning, a lineshape will be generated whether the parameters are
#' normalized or not, and the domain/range will vary as a result.
#'
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @inheritParams methodEllipse
#'
#' @return TO DO.
#'
#' @name f_lineshape
#' @export
setGeneric("f_lineshape", 
  function(object, ...) {
    standardGeneric("f_lineshape")
  })

#------------------------------------------------------------------------
#' Generate baseline function
#'
#' This is primarily an internal method that generates a function which
#' can be used to calculate the baseline at a desired chemical shift.
#' Any values outside the original fit domain will be set to 0. 
#'
#' TO DO.
#'
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @inheritParams methodEllipse
#'
#' @return TO DO.
#'
#' @name f_baseline
#' @export
setGeneric("f_baseline", 
  function(object,  ...) {
    standardGeneric("f_baseline")
  })


#------------------------------------------------------------------------
#' Calculate lineshape
#'
#' Calculate the lineshape associated with each defined peak.
#'
#' TO DO.
#'
#' @param object An NMRScaffold1D on NMRFit1D object.
#' @inheritParams methodEllipse
#'
#' @return To DO.
#'
#' @name calc_lineshape
#' @export
setGeneric("calc_lineshape", 
  function(object, ...) {
    standardGeneric("calc_lineshape")
  })

#------------------------------------------------------------------------
#' Calculate baseline
#'
#' Calculate baseline values at desired chemical shift.
#'
#' TO DO.
#'
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @inheritParams methodEllipse
#'
#' @return A vector of NMR intensity values corresponding to the baseline
#'         estimated in the fit.
#'
#' @name calc_baseline
#' @export
setGeneric("calc_baseline", 
  function(object, ...) {
    standardGeneric("calc_baseline")
  })

#------------------------------------------------------------------------
#' Calculate the areas of each peak
#'
#' Calculate the analytical area of each peak based on the peak_type.
#'
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @param type One of either 'analytical' or 'numerical'. Analytical
#'             calculation is the default, but numerical calculation using
#'             integrate() is left for debugging purposes.
#' @inheritParams methodEllipse
#'
#' @return A data.frame with three columns -- id, peak, and area.
#'
#' @name calc_area
#' @export
setGeneric("calc_area", 
  function(object, type = 'analytical', ...) {
    standardGeneric("calc_area")
  })



#========================================================================>
#  Constraint generation
#========================================================================>



#------------------------------------------------------------------------
#' Set conservative bounds on an NMRScaffold1D or NMRScaffold2D object
#'
#' Although any NMRScaffold object can act as a boundary on another
#' NMRScaffold object, this function provides a convenience method
#' for generating conservative normalized bounds. Essentially, most parameters
#' are effectively bound to an approximately 0-1 range following normalization.
#'
#' Warning, these bounds assume that the chemical shift domain is relatively
#' small (< 0.1 ppm)
#'
#' @param object An NMRScaffold1D or NMRscaffold2D object.
#' @param nmrdata An NMRData object.
#' @param position Fraction difference override applied to position parameters.
#' @param height Fraction difference override applied to height parameters.
#' @param width Fraction difference override applied to width parameters.
#' @inheritParams methodEllipse
#'
#' @return A new NMRScaffold1D or NMRScaffold2D object with modified parameters.
#'
#' @name set_conservative_bounds
#' @export
setGeneric("set_conservative_bounds", 
  function(object, nmrdata = NULL, position = NULL, 
           height = NULL, width = NULL, ...) {
    standardGeneric("set_conservative_bounds")
  })

#------------------------------------------------------------------------
#' Set relative bounds on an NMRScaffold1D or NMRScaffold2D object
#'
#' Although any NMRScaffold object can act as a boundary on another
#' NMRScaffold object, this function provides a convenience method
#' for generating bounds based on a relative difference from current
#' parameter values.
#'
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @param nmrdata An NMRData object.
#' @param overall Overall fraction difference applied to all parameters.
#' @param position Fraction difference override applied to position parameters.
#' @param height Fraction difference override applied to height parameters.
#' @param width Fraction difference override applied to width parameters.
#' @inheritParams methodEllipse
#'
#' @return A new NMRScaffold1D or NMRFit1D object with modified parameters.
#'
#' @export
setGeneric("set_relative_bounds", 
  function(object, nmrdata = NULL, overall = 0.1, position = NULL, 
           height = NULL, width = NULL, ...) {
    standardGeneric("set_relative_bounds")
  })
