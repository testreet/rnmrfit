# Definition of a class structure for 2D resonance data.

#' @import Rcpp
#' @useDynLib rnmrfit
NULL

#==============================================================================>
#  NMRResonance2D -- peak description 
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR resonance.
#' 
#' An NMRResonance2D object is nothing but a list of 2 NMRResonance1D objects,
#' one for the direct and one for the indirect dimension. All methods only
#' serve as wrappers for more conveniently working with the individual direct
#' and indirect components.
#' 
#' @slot resonances A list with "direct" and "indirect" elements, where both
#'                  elements are valid NMRResonance1D objects.
#' 
#' @name NMRResonance2D-class
#' @export
NMRResonance2D <- setClass("NMRResonance2D",
  slots = c(
    resonances = 'list'
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRResonance2D validity test
#'
validNMRResonance2D <- function(object) {

  resonances <- object@resonances
  direct <- resonances$direct
  indirect <- resonances$indirect

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Checking resonance names
  if (! any(names(resonances) %in% c('direct', 'indirect')) ) {
    valid <- FALSE
    new.err <- '"resonances" must contain "direct" and "indirect" elements.'
    err <- c(err, new.err)
  }

  #---------------------------------------
  # Both 1D elements must be valid NMRResonance1D objects
  logic1 <- ('NMRResonance1D' %in% class(direct)) && validObject(direct)
  logic2 <- ('NMRResonance1D' %in% class(indirect)) && validObject(indirect)

  if (! (logic1 && logic2) ) {
    valid <- FALSE
    new.err <- paste('"direct" and "indirect" components must be valid',
                     'NMRResonance1D objects.')
    err <- c(err, new.err)
  } else {

    # If both are valid objects check ids
    if ( direct@id != indirect@id ) {
      valid <- FALSE
      new.err <- '"direct" and "indirect" components must have the same id.'
      err <- c(err, new.err)
    }

  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRResonance2D", validNMRResonance2D)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRResonance2D object based on simplified peak list
#' 
#' Generates an NMRResonance2D object by combining direct and indirect peak
#' definitions. These peak definitions can be provided either with a string
#' that will be parsed by \code{parse_peaks_1d()} or a previously generated
#' NMRResonance1D object.
#' 
#' @param direct.peaks An NMRResonance1D object, numeric vector of singlet
#'                     chemical shifts, or a character string specifying
#'                     multiplets of the form "3 d 1.2". See ?parse_peaks_1d for
#'                     more information.
#' @param indirect.peaks An NMRResonance1D object, numeric vector of singlet
#'                       chemical shifts, or a character string specifying
#'                       multiplets of the form "3 d 1.2". See ?parse_peaks_1d
#'                       for more information.
#' @param direct.sf Sweep frequency (MHz) in the direct dimension -- needed to
#'                  convert coupling constants from Hz to ppm. In most cases, it
#'                  is recommended to set a single default value using
#'                  nmrsession_2d$sf$direct  = ..., but an override can be
#'                  provided here.
#' @param indirect.sf Sweep frequency (MHz) in the indirect dimension -- needed
#'                    to convert coupling constants from Hz to ppm. In most
#'                    cases, it is recommended to set a single default value
#'                    using nmrsession_2d$sf$indirect  = ..., but an override
#'                    can be provided here.
#' @param id A string specifying resonance name. If left empty, a name is
#'           automatically generated from the peaks argument.
#' @param width Initial estimate of peak width (in Hz). For Voigt lineshapes,
#'              this value is taken as the Lorentzian component, with the
#'              Gaussian component calculated from
#'              peak.width*frac.guass/(1-frac.gauss).
#' @param fraction.gauss Fraction of overall peak width that corresponds to a
#'                       Gaussian lineshape. A value of 0 corresponds to a
#'                       Lorentz peak whereas a value of 1 corresponds to a
#'                       Gaussian peak. Values in between 0 and 1 are modelled
#'                       as a Voigt lineshape but the specific value of
#'                       frac.gauss does not have a physical interpretation.
#' @param position.leeway A fraction specifying how tightly enforced the
#'                        coupling constraints on peak positions, should be.
#'                        E.g. coupling.leeway = 0 specifies that the j coupling
#'                        constant is exact, whereas couping.leeway = 0.1
#'                        specifies that the coupling constant may differ by +/-
#'                        10 percent.
#' @param width.leeway Similar to position.leeway but for peak widths.
#'                     Determines how strictly equal peak widths for all coupled
#'                     peaks are enforced.
#' @param area.leeway Similar to position.leeway but for peak areas. Determines
#'                    how strictly the coupling area ratios are enforced.
#' 
#' @return An NMRResonance2D object.
#' 
#' @export
nmrresonance_2d <- function(direct.peaks, indirect.peaks, 
                            direct.sf = nmrsession_2d('sf')$direct, 
                            indirect.sf = nmrsession_2d('sf')$indirect, 
                            id = NULL, width = 1, fraction.gauss = 0, 
                            position.leeway = 0, area.leeway = 0, 
                            width.leeway = 0) {

  #---------------------------------------
  # First, generating the direct/indirect components

  # If the direct.peaks component is already an NMRResonance 1D object,
  # then there is nothing to do, otherwise, pass input into nmrresonance_1d
  if ( 'NMRResonance1D' %in% class(direct.peaks) ) {
    direct <- direct.peaks
  } else {
    direct <- nmrresonance_1d(direct.peaks, sf = direct.sf, id = id,
                              width = width, fraction.gauss = fraction.gauss,
                              position.leeway = position.leeway,
                              area.leeway = area.leeway, 
                              width.leeway = width.leeway)
  }

  # Same for indirect
  if ( 'NMRResonance1D' %in% class(indirect.peaks) ) {
    indirect <- indirect.peaks
  } else {
    indirect <- nmrresonance_1d(indirect.peaks, sf = indirect.sf, id = id,
                                width = width, fraction.gauss = fraction.gauss,
                                position.leeway = position.leeway,
                                area.leeway = area.leeway, 
                                width.leeway = width.leeway)
  }

  #---------------------------------------
  # Aligning the ids

  # If an id is provided apply it directly to both components
  if (! is.null(id) ) {
    direct@id <- id
    indirect@id <- id
  }
  # If the id is not provided, and individual ids are different, combine them.
  else if ( direct@id != indirect.id ) {
    id <- paste(direct@id, indirect@id, sep = ' / ')
    direct@id <- id
    indirect@id <- id
  }

  #---------------------------------------
  # Generate object

  new('NMRResonance2D', resonances = list(direct = direct, indirect = indirect))
}



#==============================================================================>
#  Display function
#==============================================================================>



#------------------------------------------------------------------------------
#' Display NMRResonance2D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRResonance2D", 
  function(object) {

    direct <- object@resonances$direct
    indirect <- object@resonances$indirect

    cat('==================================\n')
    cat('An object of NMRResonance2D class:\n\n')


    cat('In the direct dimension:\n\n')
    cat('----------------------------------\n')
    show(direct)

    cat('In the indirect dimension:\n\n')
    cat('----------------------------------\n')
    show(indirect)

  })



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
# Direct

#---------------------------------------
#' Get direct dimension resonance object
#' 
#' Generic convenience method to access the direct dimension resonance
#' component of an NMRResonance2D or NMRSpecies2D object.
#' 
#' @param object An NMRResonance2D or NMRSpecies2D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name direct
#' @export
setGeneric("direct", 
  function(object, ...) standardGeneric("direct")
  )

#' @rdname direct
#' @export
setMethod("direct", "NMRResonance2D", 
  function(object) object@resonances$direct)

#---------------------------------------
#' Set direct dimension resonance object
#' 
#' Generic convenience method to set the direct dimension resonance component
#' of an NMRResonance2D or NMRSpecies2D object.
#' 
#' @param object An NMRResonance2D or NMRSpecies2D object.
#' @param value A valid NMRResonance1D object.
#' 
#' @name direct-set
#' @export
setGeneric("direct<-", 
  function(object, value) standardGeneric("direct<-"))

#' @rdname direct-set
#' @export
setReplaceMethod("direct", "NMRResonance2D",
  function(object, value) {
    object@resonance$direct <- value
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Indirect

#---------------------------------------
#' Get indirect dimension resonance object
#' 
#' Generic convenience method to access the indirect dimension resonance
#' component of an NMRResonance2D or NMRSpecies2D object.
#' 
#' @param object An NMRResonance2D or NMRSpecies2D object.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name indirect
#' @export
setGeneric("indirect", 
  function(object, ...) standardGeneric("indirect")
  )

#' @rdname indirect
#' @export
setMethod("indirect", "NMRResonance2D", 
  function(object) object@resonances$indirect)

#---------------------------------------
#' Set indirect dimension resonance object
#' 
#' Generic convenience method to set the indirect dimension resonance component
#' of an NMRResonance2D or NMRSpecies2D object.
#' 
#' @param object An NMRResonance2D or NMRSpecies2D object.
#' @param value A valid NMRResonance1D object.
#' 
#' @name indirect-set
#' @export
setGeneric("indirect<-", 
  function(object, value) standardGeneric("indirect<-"))

#' @rdname indirect-set
#' @export
setReplaceMethod("indirect", "NMRResonance2D",
  function(object, value) {
    object@resonance$indirect <- value
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Id

#' @rdname id
#' @export
setMethod("id", "NMRResonance2D", 
  function(object) object@resonances$direct@id)

#' @rdname id-set
#' @export
setReplaceMethod("id", "NMRResonance2D",
  function(object, value) {
    id <- as.character(value)
    object@resonances$direct@id <- id
    object@resonances$indirect@id <- id
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Peaks

#' @rdname peaks
#' @export
setMethod("peaks", "NMRResonance2D", 
  function(object, include.id = FALSE) {
    lapply(peaks, object@resonances, include.id = include.id)
  })

#' @export
setReplaceMethod("peaks", "NMRResonance2D",
  function(object, value) {
    err <- paste('There is no peak setter function for NMRResonance2D objects.',
                 '"direct" and "indirect" NMRResonance1D objects must be set',
                 'explicitly.')
    stop(err)
  })



#------------------------------------------------------------------------------
# Couplings

#---------------------------------------
#' Get object couplings 
#' 
#' Generic convenience method to access the coupling definitions of an
#' NMRResonance1D object. If used on an NMRSpecies1D or NMRFit1D object, the
#' function combines the couplings data frames of component resonances.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param include.id TRUE to return a column of resonance or species ids.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name couplings
#' @export
setGeneric("couplings", 
  function(object, include.id = FALSE, ...) standardGeneric("couplings")
  )

#' @rdname couplings
#' @export
setMethod("couplings", "NMRResonance1D", 
  function(object, include.id = FALSE) {
    couplings <- object@couplings
    if ( include.id && (nrow(couplings) > 0) ) {
      cbind(resonance.1 = object@id, resonance.2 = object@id, object@couplings)
    }
    else object@couplings
  })

#---------------------------------------
#' Set object couplings
#' 
#' Generic convenience method to set the coupling definitions of an
#' NMRResonance1D object.
#' 
#' @param object An NMRResonance1D object.
#' @param value A data frame with "id.1", "id.2", "position.difference", and
#'              "area.ratio" columns.
#' 
#' @name couplings-set
#' @export
setGeneric("couplings<-", 
  function(object, value) standardGeneric("couplings<-")
  )

#' @rdname couplings-set
#' @export
setReplaceMethod("couplings", "NMRResonance1D",
  function(object, value) {
    object@couplings <- value
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Bounds

#---------------------------------------
#' Get object bounds
#' 
#' Generic convenience method to access the bounds of an NMRResonance1D object.
#' If used on an NMRSpecies1D or NMRFit1D object, the function combines the
#' bounds data frames of component resonances.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param include.id TRUE to return a column of resonance or species ids.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name bounds
#' @export
setGeneric("bounds", 
  function(object, include.id = FALSE, ...) standardGeneric("bounds")
  )

#' @rdname bounds
#' @export
setMethod("bounds", "NMRResonance1D", 
  function(object, include.id = FALSE) {
    bounds <- .initialize_bounds(object)@bounds

    lower <- bounds$lower
    upper <- bounds$upper

    if ( include.id ) {
      if ( nrow(lower) > 0 ) {
        bounds$lower <- cbind(resonance = object@id, lower)
      }

      if ( nrow(upper) > 0 ) {
        bounds$upper <- cbind(resonance = object@id, upper)
      }
    }
    bounds
  })

#---------------------------------------
#' Set object bounds
#' 
#' Generic convenience method to set the bounds of an NMRResonance1D object.
#' 
#' @param object An NMRResonance1D object.
#' @param value A list with "lower" and "upper" elements, each containing a data
#'              frame with "peak", "position", "width", "height", and
#'              "fraction.gauss" columns.
#' 
#' @name bounds-set
#' @export
setGeneric("bounds<-", 
  function(object, value) standardGeneric("bounds<-"))

#' @rdname bounds-set
#' @export
setReplaceMethod("bounds", "NMRResonance1D",
  function(object, value) {
    object@bounds <- value
    validObject(object)
    object 
  })



#==============================================================================>
#  Initialization functions (generating parameter estimates based on data)
#==============================================================================>



#------------------------------------------------------------------------------
#' Initialize peak heights of an NMRResonance1D object
#' 
#' Generates peak height estimates based on spectral data. If applied to an
#' NMRSpecies1D or NMRFit1D object, the same initialization is propagated to
#' every component resonance. Since estimates cannot be provided for any peak
#' outside the given data range, any such peaks are ignored. Whether or not
#' warning/error messages are generated when that occurs is specified by the
#' exclusion.notification parameter. Similarly, exclusion.level provided
#' options to omit the whole resonance/species if a peak is found to be outside
#' the data bounds.
#' 
#' At this point, there is just one approach: take peak height as the intensity
#' of the data at the current position of the peak. There are plans to develop
#' more sophisticated approaches in the future.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param nmrdata An NMRData1D object.
#' @param exclusion.level A string specifying what to do when peaks are found to
#'                        fall outside of the data range: either 'species' to
#'                        exclude the whole species to which the offending peak
#'                        belongs, 'resonance' to exclude the resonance to which
#'                        the offending peak belongs, or 'peak' to exclude just
#'                        the peak itself.
#' @param exclusion.notification A function specifying how to report when peaks
#'                               are found to be outside the data range: 'none'
#'                               to ignore, 'message' to issue a message,
#'                               'warning' to issue a warning, and 'stop' to
#'                               issue an error.
#' @inheritParams methodEllipse
#' 
#' @return A new object with modified peak heights.
#' 
#' @name initialize_heights
#' @export
setGeneric("initialize_heights", 
  function(object, ...) {
    standardGeneric("initialize_heights")
  })

#' @rdname initialize_heights
#' @export
setMethod("initialize_heights", "NMRResonance1D",
  function(object, nmrdata, exclusion.level = nmrsession_1d$exclusion$level,
           exclusion.notification = nmrsession_1d$exclusion$notification) {

  # Checking nmrdata
  if ( class(nmrdata) != 'NMRData1D' ) {
    err <- '"nmrdata" must be a valid NMRData1D object.'
    stop(err)
  }
  else {
    validObject(nmrdata)
  }

  # Building an interpolating function betwewn chemical shift and intensity
  d <- processed(nmrdata)
  f <- approxfun(d$direct.shift, Re(d$intensity))

  # Excluding peaks that are outside the data frame
  peaks <- peaks(object) 
  logic <- (peaks$position > min(d$direct.shift)) & 
           (peaks$position < max(d$direct.shift))

  peaks <- peaks[logic, ]

  # Generating heights from interpolation
  peaks$height <- f(peaks$position)

  # Updating
  object.2 <- update_peaks(object, peaks, exclusion.level = exclusion.level,
                           exclusion.notification = exclusion.notification)
  peaks.2 <- peaks(object.2)
  
  all.columns <- colnames(peaks)
  data.columns <- c('position', 'width', 'height', 'fraction.gauss')
  id.columns <- all.columns[! all.columns %in% data.columns]

  peaks <- peaks(object)
  ids <- apply(peaks[, id.columns], 1, paste, collapse = '-')
  ids.2 <- apply(peaks.2[, id.columns], 1, paste, collapse = '-')

  peaks[ids %in% ids.2, ] <- peaks.2
  peaks(object) <- peaks

  object
})



#==============================================================================>
#  Bounds
#==============================================================================>



#------------------------------------------------------------------------------
#' Set general bounds of an NMRResonance1D object
#' 
#' This function provides a convenience method for generating bounds using a
#' simple set of lower and upper constraints on basic peak parameters such as
#' position, height, width, fraction.gauss. The term "general" refers to the
#' fact that the same bounds are applied to each and every peak, regardless of
#' current parameter values. These bounds can be normalized to a set of data
#' using the optional nmrdata argument. If applied to an NMRSpecies1D or
#' NMRFit1D object, the same bounds are propagated to every component
#' resonance.
#' 
#' In practice, general bounds are primarily useful for placing a hard
#' constraint on peak widths and preventing negative heights. Values of 0 for
#' widths and height can also cause issues during optimization, so simple
#' general bounds can be used to prevent errors.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param position A vector of two elements corresponding to a lower and upper
#'                 bound for peak position. If nmrdata is provided, 0
#'                 corresponds to the leftmost range of the data and 1 to the
#'                 rightmost. Otherwise, the units are in ppm.
#' @param height A vector of two elements corresponding to a lower and upper
#'               bound for peak height. If nmrdata is provided, 0 corresponds to
#'               the lowest value of spectral intensity and 1 to the largest.
#'               Otherwise, the units correspond to arbitrary spectral intensity
#'               values.
#' @param width A vector of two elements corresponding to a lower and upper
#'              bound for peak width in Hz. If nmrdata is provided, values are
#'              taken as fraction of the general data range. So 0.1 would
#'              correspond to a nominal peak width that covers a tenth of the
#'              general data range.
#' @param fraction.gauss A vector of two elements corresponding to a lower and
#'                       upper bound for the Gaussian fraction of the peak. This
#'                       can be set to c(0, 0) to force Lorentzian peaks. Any
#'                       values smaller than 0 will be treated as 0 and any
#'                       values greater than 1 will be treated as 1.
#' @param baseline Only applicable to NMRFit1D objects. A complex vector of two
#'                 elements corresponding to a lower and upper bound for both
#'                 real and imaginary baseline control points (which roughly
#'                 corresponds to a baseline value). If the vector has no
#'                 imaginary component, the imaginary baseline is left
#'                 unbounded. If nmrdata is provided, 0 corresponds to the
#'                 lowest value of spectral intensity and 1 to the largest.
#'                 Otherwise, the units correspond to arbitrary spectral
#'                 intensity values.
#' @param phase Only applicable to NMRFit1D objects. A vector of two elements
#'              corresponding to a lower and upper bound for the phase
#'              correction (in radians) at any point in the chemical shift
#'              range.
#' @param nmrdata An optional NMRData1D object that can serve as a reference
#'                point for the bounds.
#' @param widen FALSE to prevent new bounds from widening existing bounds.
#' @inheritParams methodEllipse
#' 
#' @return A new object with modified bounds.
#' 
#' @name set_general_bounds
#' @export
setGeneric("set_general_bounds", 
  function(object, ...) {
    standardGeneric("set_general_bounds")
  })

#' @rdname set_general_bounds
#' @export
setMethod("set_general_bounds", "NMRResonance1D",
  function(object, position = NULL, height = NULL, width = NULL,
           fraction.gauss = NULL, nmrdata = NULL, widen = FALSE) {
  
  # Initializing bounds
  object <- .initialize_bounds(object)
  lower <- object@bounds$lower
  upper <- object@bounds$upper

  #---------------------------------------
  # Scaling all bounds if nmrdata has been provided
  if (! is.null(nmrdata) ) {

    if ( class(nmrdata) != 'NMRData1D' ) {
      err <- '"nmrdata" must be a valid NMRData1D object.'
      stop(err)
    }
    else {
      validObject(nmrdata)
    }

    processed <- nmrdata@processed
    y.range <- max(Re(processed$intensity)) - min(Re(processed$intensity))
    x.range <- max(processed$direct.shift) - min(processed$direct.shift)

    position <- position * x.range + min(processed$direct.shift)
    height <- height * y.range

    sfo1 <- get_parameter(nmrdata, 'sfo1', 'procs')
    width <- width * (x.range[2] - x.range[1]) * sfo1
  }

  #---------------------------------------
  # Defining a bound check function
  .check_bounds <- function(bounds) {

    if ( length(bounds) != 2 ) {
      err <- paste("All bounds must be vectors of two elements consisting",
                   "of a lower and upper bound.")
      stop(err)
    }

    if ( bounds[1] > bounds[2] ) {
      err <- paste("Lower bound must be smaller than upper bound.",
                   "Proceeding with current constraints will result in a",
                   "fit error.")
      warning(err)
    }

  }

  #---------------------------------------
  # Creating a list of bounds to loop through each in term
  bounds = list(position = position, height = height, width = width,
                fraction.gauss = fraction.gauss)

  for ( parameter in names(bounds) ) {
    if ( length(bounds[[parameter]]) > 0 ) {
      .check_bounds(bounds[[parameter]])
      lower[[parameter]] <- bounds[[parameter]][1]
      upper[[parameter]] <- bounds[[parameter]][2]
    }
  }

  # Fraction gauss is a little different because it must be 0-1
  lower$fraction.gauss[lower$fraction.gauss < 0] <- 0
  upper$fraction.gauss[upper$fraction.gauss > 0] <- 1

  # Ensuring that parameters are only widened if desired
  columns <- c('position', 'width', 'height', 'fraction.gauss')

  new.lower <- unlist(lower[ , columns])
  old.lower <- unlist(object@bounds$lower[ , columns])

  new.upper <- unlist(upper[ , columns])
  old.upper <- unlist(object@bounds$upper[ , columns])
  
  if (! widen ) {
    new.lower <- ifelse(new.lower < old.lower, old.lower, new.lower)
    new.upper <- ifelse(new.upper > old.upper, old.upper, new.upper)
  }

  object@bounds$lower[ , columns] <- new.lower
  object@bounds$upper[ , columns] <- new.upper

  validObject(object)
  object
})



#------------------------------------------------------------------------------
#' Set offset bounds of an NMRResonance1D object
#' 
#' This function provides a convenience method for generating bounds using a
#' simple set of lower and upper constraints on basic peak parameters such as
#' position, height, width, fraction.gauss. The term "offset" refers to the
#' fact that bounds are applied as an offset to current values of the
#' parameters. These bounds can be expressed in absolute (e.g. -0.1 and +0.1
#' ppm) or relative (e.g. -1 percent and +1 percent) terms. If applied to an
#' NMRSpecies1D or NMRFit1D object, the same offset bounds are propogated to
#' every component resonance.
#' 
#' In practice, offset bounds are primarily useful for preventing peak
#' positions from drifting too much from initial guesses and for fine-tuning a
#' fit once an initial optimization is performed. It is not recommended to use
#' strict offset bounds based on rough initial parameter guesses.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param position A vector of two elements to be added to current peak
#'                 positions to generate a set of lower and upper bounds. If
#'                 relative is true, the values are treated as fractions to be
#'                 multipled by the current peak position before addition.
#'                 fraction of the current position.
#' @param height A vector of two elements to be added to current peak heights to
#'               generate a set of lower and upper bounds. If relative is true,
#'               the values are treated as fractions to be multipled by the
#'               current peak heights before addition.
#' @param width A vector of two elements to be added to current peak widths to
#'              generate a set of lower and upper bounds. If relative is true,
#'              the values are treated as fractions to be multipled by the
#'              current peak heights before addition.
#' @param relative TRUE to treat values as relative fractions, FALSE to apply
#'                 them directly.
#' @param widen FALSE to prevent new bounds from widening existing bounds.
#' @inheritParams methodEllipse
#' 
#' @return A new object with modified bounds.
#' 
#' @name set_offset_bounds
#' @export
setGeneric("set_offset_bounds", 
  function(object, ...) {
    standardGeneric("set_offset_bounds")
  })

#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRResonance1D",
  function(object, position = NULL, height = NULL, width = NULL, 
           relative = FALSE, widen = FALSE) {

  # Initializing bounds
  object <- .initialize_bounds(object)
  peaks <- object@peaks
  lower <- object@bounds$lower
  upper <- object@bounds$upper

  #---------------------------------------
  # Defining a bound check function
  .check_bounds <- function(bounds) {

    if ( length(bounds) != 2 ) {
      err <- paste("All bounds must be vectors of two elements consisting",
                   "of a lower and upper bound.")
      stop(err)
    }

    err2 <- "Proceeding with current constraints will result in a fit error."

    if ( bounds[1] > 0 ) {
      err <- paste("Lower offsets must be negative so that resulting bounds",
                   "include initial values.", err2)
      stop(err)
    }

    if ( bounds[2] < 0 ) {
      err <- paste("Upper offsets must be positive so that resulting bounds",
                   "include initial values.", err2)
      stop(err)
    }

    if ( bounds[1] > bounds[2] ) {
      err <- paste("Lower bound must be smaller than upper bound.", err2)
      warning(err)
    }

  }

  #---------------------------------------
  # Creating a list of bounds to loop through each in term
  bounds = list(position = position, height = height, width = width)

  for ( parameter in names(bounds) ) {
    if ( length(bounds[[parameter]]) > 0 ) {
      .check_bounds(bounds[[parameter]])

      lower.offset <- bounds[[parameter]][1]
      upper.offset <- bounds[[parameter]][2]
      
     if ( relative ) {
        lower.offset <- lower.offset*peaks[[parameter]]
        upper.offset <- upper.offset*peaks[[parameter]]
      } 

      lower[[parameter]] <- peaks[[parameter]] + lower.offset
      upper[[parameter]] <- peaks[[parameter]] + upper.offset
    }
  }

  # Ensuring that parameters are only widened if desired
  columns <- c('position', 'width', 'height', 'fraction.gauss')

  new.lower <- unlist(lower[ , columns])
  old.lower <- unlist(object@bounds$lower[ , columns])

  new.upper <- unlist(upper[ , columns])
  old.upper <- unlist(object@bounds$upper[ , columns])
  
  if (! widen ) {
    new.lower <- ifelse(new.lower < old.lower, old.lower, new.lower)
    new.upper <- ifelse(new.upper > old.upper, old.upper, new.upper)
  }

  object@bounds$lower[ , columns] <- new.lower
  object@bounds$upper[ , columns] <- new.upper

  validObject(object)
  object
})



#------------------------------------------------------------------------------
#' Set conservative bounds on an NMRResonance1D object
#' 
#' A convenience function that sets reasonable bounds on the fit. These bounds
#' are assumed to be widely applicable to most simple NMR data. Each set of
#' bounds can be turned on or off as necessary. A slightly better set of bounds
#' can be selected if a reference NMRData1D object is provided. If applied to
#' an NMRSpecies1D or NMRFit1D object, the same conservated bounds are
#' propogated to every component resonance.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param position Without reference data, position is limited to plus or minus
#'                 0.1 ppm. With reference data, the position of the peaks is
#'                 forced inside the domain of the data. FALSE to disable.
#' @param height Without reference data, height is set to strictly positive.
#'               With reference data, height is also limited to no more than 150
#'               percent of the maximum peak value. FALSE to disable.
#' @param width Without reference data, minimum peak width is set to almost, but
#'              not quite 0 Hz (1e-3 Hz) and a maximum peak width of 3 Hz. With
#'              reference data, peak width is prevented from being more than 20
#'              percent of the data range.  FALSE to disable.
#' @param baseline Without reference data, the baseline control points are left
#'                 unrestricted. With reference data, the real baseline control
#'                 points are limited to no more than 50 percent of the maximum
#'                 peak value. FALSE to disable.
#' @param phase With or without reference data, phase correction is limited to
#'              -pi/2 to pi/2.
#' @param nmrdata An optional NMRData1D object that can serve as a reference
#'                point for the bounds.
#' @param widen FALSE to prevent new bounds from widening existing bounds.
#' @inheritParams methodEllipse
#' 
#' @return A new object with modified bounds.
#' 
#' @name set_conservative_bounds
#' @export
setGeneric("set_conservative_bounds", 
  function(object, ...) {
    standardGeneric("set_conservative_bounds")
})

#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRResonance1D",
  function(object, position = TRUE,  height = TRUE, width = TRUE, 
           nmrdata = NULL, widen = FALSE) { 

  # First, do a single pass over general bounds with no reference
  if ( height )  gen.height <- c(0, Inf)
  else gen.height <- NULL

  if ( width ) gen.width <- c(0.003, 3)
  else gen.width <- NULL

  object <- set_general_bounds(object, height = gen.height, width = gen.width,
                               widen = widen)

  # Adding position offsets
  if ( position ) {
    object <- set_offset_bounds(object, position = c(-0.1, 0.1), widen = widen)
  }

  # If nmrdata is provided, add further constraints  
  if (! is.null(nmrdata) ) {
    
    if ( class(nmrdata) != 'NMRData1D' ) {
      err <- '"nmrdata" must be a valid NMRData1D object.'
      stop(err)
    } else {
      validObject(nmrdata)
    }

    if ( position )  gen.position <- c(0, 1)
    else gen.position <- NULL

    if ( height ) gen.height <- c(0, 1.5)
    else gen.height <- NULL

    if ( width ) gen.width <- c(0, 0.2)
    else gen.width <- NULL

    object <- set_general_bounds(object, position = gen.position, 
                                 height = gen.height, width = gen.width,
                                 nmrdata = nmrdata, widen = widen)
  }

  object
  })



#------------------------------------------------------------------------------
#' Set peak type of an NMRResonance1D object
#' 
#' The peak type of an NMRResonance1D object is governed by the fraction.gauss
#' parameter and can fluidly go from pure Lorentz to Gauss via the combined
#' Voigt lineshape. However, it can be computationally efficient to force a
#' specific peak type. This function provided a shortcut for doing so by
#' changing the fraction.gauss parameter as well as lower and upper bounds.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param peak.type One of either "lorentz", "voigt", "gauss", or "any" where
#'                  "any" clears existing bounds on the fraction.gauss
#'                  parameter.
#' @inheritParams methodEllipse
#' 
#' @return A new object with modified bounds and peak parameters.
#' 
#' @name set_peak_type
#' @export
setGeneric("set_peak_type", 
  function(object, ...) {
    standardGeneric("set_peak_type")
  })

#' @rdname set_general_bounds
#' @export
setMethod("set_peak_type", "NMRResonance1D",
  function(object, peak.type) {

    # Initializing bounds
    object <- .initialize_bounds(object)
    peaks <- object@peaks
    lower <- object@bounds$lower
    upper <- object@bounds$upper

    # Getting rid of empty spaces and capitals
    peak.type <- tolower(gsub('\\s', '', peak.type))
    peak.types <- c('lorentz', 'voigt', 'gauss', 'any')
    peak.type <- pmatch(peak.type, peak.types)

    if ( peak.type == 1 ) {
      lower$fraction.gauss <- 0
      upper$fraction.gauss <- 0
      peaks$fraction.gauss <- 0
    } else if ( peak.type == 2 ) {
      lower$fraction.gauss <- 1e-6
      upper$fraction.gauss <- 1 - 1e-6

      logic <- peaks$fraction.gauss < lower$fraction.gauss
      peaks[logic , 'fraction.gauss'] <- lower$fraction.gauss[logic] + 1e-6
      logic <- peaks$fraction.gauss > upper$fraction.gauss
      peaks[logic , 'fraction.gauss'] <- upper$fraction.gauss[logic] - 1e-6
    } else if ( peak.type == 3 ) {
      lower$fraction.gauss <- 1
      upper$fraction.gauss <- 1
      peaks$fraction.gauss <- 1
    } else if ( peak.type == 4 ) {
      lower$fraction.gauss <- 0
      upper$fraction.gauss <- 1
    } else {
      peak.types <- paste(peak.types, collapse = ', ')
      err <- sprintf('Peak type must be one of %s', peak.types)
      stop(err)
    }

    object@bounds <- list(lower = lower, upper = upper)
    object@peaks <- peaks

    object
  })



#========================================================================>
#  Lineshape and area calculations
#========================================================================>



#------------------------------------------------------------------------
#' Generate lineshape function
#' 
#' This is primarily an internal method that outputs a function (or a tbl_df
#' data frame of functions), where each function outputs spectral intensity
#' data given a vector input of chemical shifts.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmrsession_1d(sf = ...), but an override can be
#'           provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single function, FALSE to output a data frame of functions
#'                  that correspond to individual peaks.
#' @param include.id TRUE to include id as "resonance" column if outputting data
#'                   frame.
#' @param components 'r/i' to output both real and imaginary data, 'r' to output
#'                   only real and 'i' to output only imaginary.
#' @inheritParams methodEllipse
#' 
#' @return A function or tbl_df data frame of functions where each function
#'         outputs spectral intensity data given a vector input of chemical
#'         shifts. In the latter case, the functions are stored in a list column
#'         called f.
#' 
#' @name f_lineshape
#' @export
setGeneric("f_lineshape", 
  function(object, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           include.id = FALSE, components = 'r/i', ...) {
    standardGeneric("f_lineshape")
})

#' @rdname f_lineshape
#' @export
setMethod("f_lineshape", "NMRResonance1D",
  function(object, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           include.id = FALSE, components = 'r/i') {

    # Checking to make sure that sweep frequency is defined
    err <- '"sf" must be provided as input or set using nmrsession_1d()'
    if ( is.null(sf) ) stop(err)

    # Defining which components to return
    return.r <- grepl('r', tolower(components))
    return.i <- grepl('i', tolower(components))

    err <- '"components" must have at least one of either "r" or "i"'
    if ( return.r && return.i ) f_out <- function(y) {y}
    else if ( return.r ) f_out <- function(y) {Re(y)}
    else if ( return.i ) f_out <- function(y) {Im(y)}
    else stop(err)

    columns <- c('position', 'width', 'height', 'fraction.gauss')
    peaks <- peaks(object)
    parameters <- as.matrix(peaks[, columns])

    # Converting peak width to ppm
    parameters[, 2] <- parameters[, 2]/sf

    # If peaks are to be summed, just feed all parameters into the Rcpp function
    if ( sum.peaks ) {
      out <- function(x) {
        y <- .Call('_rnmrfit_lineshape_1d', PACKAGE = 'rnmrfit', x, parameters)
        f_out(y)
      }
    } 
    # Otherwise, generate a tbl_df data frame
    else {
      out <- as_tibble(peaks[, which(! colnames(peaks) %in% columns)])
      if ( include.id && (nrow(out) > 0) ) {
        if ( 'resonance' %in% colnames(out) ) {
          out <- cbind(species = object@id, out)
        }
        else {
          out <- cbind(resonance = object@id, out)
        }
      }

      parameters <- split(parameters, 1:nrow(parameters))
      
      # Generating a list of functions, each with their parameters enclosed
      functions <- lapply(parameters, function (p) {
        function(x) {
          p <- matrix(p, nrow = 1)
          y <- .Call('_rnmrfit_lineshape_1d', PACKAGE = 'rnmrfit', x, p)
          f_out(y)
        }
      })

      # Adding functions as a column
      out$f <- functions
    }

    out
    })



#------------------------------------------------------------------------
#' Calculate peak lineshape values
#' 
#' Calculated peak intensity values over a set of chemical shifts.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param direct.shift Vector of chemical shift data in ppm.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmrsession_1d(sf = ...), but an override can be
#'           provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single set of values, FALSE to output a data frame of values
#'                  that correspond to individual peaks.
#' @param sum.baseline TRUE to add baseline to every peak, if one is defined.
#'                     FALSE to exclude baseline. that correspond to individual
#'                     peaks.
#' @param include.id TRUE to include id as "resonance" column if outputting data
#'                   frame.
#' @param components 'r/i' to output both real and imaginary data, 'r' to output
#'                   only real and 'i' to output only imaginary.
#' @inheritParams methodEllipse
#' 
#' @return A vector of spectral intensity data or a data frame with columns
#'         "resonance" (optional), "peak", "direct.shift", and "intensity".
#' 
#' @name values
#' @export
setGeneric("values", 
  function(object, direct.shift, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           sum.baseline = FALSE, include.id = FALSE, components = 'r/i', ...) {
    standardGeneric("values")
})

#' @rdname values
#' @export
setMethod("values", "NMRResonance1D",
  function(object, direct.shift, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           sum.baseline = FALSE, include.id = FALSE, components = 'r/i') {

  # Generating baseline if necessaru
  if ( sum.baseline && (class(object) == 'NMRFit1D') ) {
    f <- f_baseline(object, components)
    baseline <- f(direct.shift)
  }
  else baseline <- rep(0, length(direct.shift))

  # Output depends on whether peaks are summed or not
  if ( sum.peaks ) {
    # Get function
    f <- f_lineshape(object, sf, sum.peaks, components)

    # And apply it to specified chemical shifts
    f(direct.shift) + baseline
  } 
  else {
    # Get data frame of functions
    d <- f_lineshape(object, sf, sum.peaks, include.id, components)

    # Defining function that generates necessary data frame
    f <- function(g) {
      data.frame(direct.shift = direct.shift, 
                 intensity = g[[1]](direct.shift) + baseline)
    }

    # And apply it for every peak
    group_by_if(d, function(x) {!is.list(x)}) %>% do( f(.$f) )
  }
  })



#------------------------------------------------------------------------
#' Calculate peak areas
#' 
#' Calculate total peak areas based on peak parameters.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmrsession_1d(sf = ...), but an override can be
#'           provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single area, FALSE to output a data frame of peak area
#'                  values.
#' @param include.id TRUE to include id as "resonance" column if outputting data
#'                   frame.
#' @inheritParams methodEllipse
#' 
#' @return A single overall area or a data frame of areas with columns
#'         "resonance" (optional), "peak", and "area".
#' 
#' @name areas
#' @export
setGeneric("areas", 
  function(object, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           include.id = FALSE, components = 'r/i', ...) {
    standardGeneric("areas")
})

#' @rdname areas 
#' @export
setMethod("areas", "NMRResonance1D",
  function(object, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           include.id = FALSE, components = 'r/i') {

  # Defining area function
  f <- function(position, width, height, fraction.gauss) {
    # If fraction is 0, treat as Lorentz
    if ( fraction.gauss == 0 ) {
      pi*width*height
    }
    # If fraction is 1, treat as Gauss
    else if ( fraction.gauss == 1) {
      sqrt(2*pi)*width*height
    }
    # Else, proceed as Voigt
    else {
      l.width <- width
      g.width <- width*fraction.gauss/(1 - fraction.gauss)
      Re(sqrt(2*pi)*g.width*height /
         Faddeeva_w(complex(im = l.width)/(sqrt(2)*g.width)))
    }
  }

  # Calculating areas
  areas <- peaks(object, include.id)

  areas <- areas %>%
    group_by_all() %>%
    summarize(area = f(position, width, height, fraction.gauss)) %>%
    ungroup() %>%
    select(-position, -width, -height, -fraction.gauss) %>%
    as.data.frame()

  # Sum if necessary
  if ( sum.peaks ) sum(areas$area)
  else areas
  })
