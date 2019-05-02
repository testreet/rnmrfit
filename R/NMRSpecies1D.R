# Definition of a class structure for collections of 1D resonances.



#==============================================================================>
#  NMRSpecies1D -- collection of NMRSpecies1D objects
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR species.
#' 
#' Essentially, this class is used to collect and define area relations between
#' multiple resonances as part of a single species.
#' 
#' @slot resonances A list of NMRSpecies1D objects.
#' @slot connections A data.frame relating the areas of the resonances.
#' @slot connections.leeway A value specifying how tightly enforced the
#'                          connection constraints on resonance areas should be.
#'                          E.g. a value of = 0 specifies that the area ratios
#'                          are exact, whereas 0.1 specifies that the
#'                          areas of the resonances may differ by +/- 10
#'                          percent from the specified ratios.
#' 
#' @name NMRSpecies1D-class
#' @export
NMRSpecies1D <- setClass("NMRSpecies1D",
  slots = c(
    resonances = 'list',
    connections = 'data.frame',
    connections.leeway = 'numeric',
  ),
  prototype = prototype(
    resonances = list(),
    connections = data.frame(),
    connections.leeway = 0
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRSpecies1D validity test
#'
validNMRSpecies1D <- function(object) {

  resonances <- object@resonances
  connections <- object@connections 

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Checking that all resonance list items are valid
  for ( resonance in resonances ) {
    logic1 <- class(resonance) != 'NMRResonance1D'
    logic2 <- ! validObject(resonance)
    if ( logic1 || logic2 ) {
      valid <- FALSE
      new.err <- paste('All elements of "resonances" list must be valid',
                       'NMRResonance1D objects')
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking connections 
  if ( nrow(connections) > 0 ) {

    valid.columns <- c('resonance.1', 'resonance.2', 'area.ratio')
    if (! identical(colnames(couplings), valid.columns) ) {
      valid <- FALSE
      new.err <- sprintf('"connections" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      err <- c(err, new.err)
    }

    valid.values <- names(resonances) 
    logic1 <- all(connections$resonance.1 %in% valid.values)
    logic2 <- all(connections$resonance.2 %in% valid.values)

    if (! (logic1 & logic2) ) {
      valid <- FALSE
      new.err <- paste('"connections" resonances must correspond to named',
                       'elements of "resonances" list.')
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRSpecies1D", validNMRSpecies1D)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRSpecies1D object
#' 
#' Generates an NMRSpecies1D object from a list of NMRResonance1D objects or a
#' list of character/numeric vectors that can be converted to NMRResonance1D
#' objects. See ?nmrresonance_1d for more details about this conversion.
#' 
#' @param resonances A list of NMRResonance1D objects or a list of
#'                   character/numeric vectors that can be converted to
#'                   NMRResonance1D objects. See ?nmrresonance_1d for more
#'                   details about this conversion. If list elements are named,
#'                   these names will be use to replace resonance names.
#' @param areas A vector of areas corresponding to the expected areas of the
#'              resonances. Set to NULL by default, signifying no fixed
#'              constraings.
#' @param name A string specifying species name. If left empty, a name is
#'             automatically generated from the resonance names.
#' @param connections.leeway A value specifying how tightly enforced the
#'                           connection constraints on resonance areas should
#'                           be. E.g. a value of = 0 specifies that the area
#'                           ratios are exact, whereas 0.1 specifies that the
#'                           areas of the resonances may differ by +/- 10
#'                           percent from the specified ratios.
#' 
#' @return An NMRSpecies1D object.
#' 
#' @export
nmrspecies_1d <- function(resonances, areas = NULL, name = NULL, 
                          connections.leeway = 0) {

  #---------------------------------------
  # Building peak list

  # Couplings are not added in every case
  add.couplings <- FALSE

  # If peaks is a character, parse coupling information
  if ( is.character(peaks) ) {
    coupling <- parse_peaks_1d(peaks)
    add.constraints <- TRUE

    # Initializing singlet at chemical shift
    if ( is.null(name) ) name <- peaks
    peaks <- data.frame(resonance = name, peak = 1, 
                        position = coupling$direct.shift,
                        width = width, height = 1, 
                        fraction.gauss = fraction.gauss)

    # If there is splitting to do, convert constants from Hz to ppm
    if ( any(coupling$number > 1) ) {

      # Checking to make sure that sweep frequency is defined
      err <- '"sf" must be provided as input or set using nmrsession_1d()'
      if ( is.null(sf) ) stop(err)

      # Converting coupling constant from Hz to ppm
      coupling$constant <- coupling$constant/sf

    }

    # Looping through the coupling to split the specified peaks
    for ( i in 1:length(coupling$number) ) {
      peaks <- split_peaks_1d(peaks, coupling$number[i], coupling$constant[i])
    }

    # Set flag to add couplings later
    add.couplings <- TRUE
  }
  # Otherwise, build peaks directly from singlets
  else {
    add.constraints <- FALSE
    middle <- (max(peaks) + min(peaks))/2 
    range <- paste(min(peaks), '..', max(peaks), sep = '')
    if ( is.null(name) ) name <- paste(middle, 'm', range)
    peaks <- data.frame(resonance = name, peak = 1:length(peaks), 
                        position = peaks, width = peak.width, height = 1, 
                        fraction.gauss = fraction.gauss)
  }

  #---------------------------------------
  # Adding coupling definitions

  # Starting with blanks first
  couplings <- data.frame()
  couplings.leeway = list(position = position.leeway, width = width.leeway,
                          area = area.leeway)

  nmrresonance = new('NMRSpecies1D', peaks = peaks, couplings = couplings, 
                                       couplings.leeway = couplings.leeway)

  # And then updating if necessary
  if ( add.couplings ) enforce_couplings_1d(nmrresonance)
  else nmrresonance

}



#==============================================================================>
#  Display function
#==============================================================================>



#------------------------------------------------------------------------------
#' Initialize empty set of bounds
#' 
#' Initalize lower and upper bounds with the correct dimensions, but set to
#' -Inf and +Inf respectively
#' 
#' @param object NMRSpecies1D object
#' @param overwrite TRUE to overwrite existing bounds (to reset them), FALSE to
#'                  quietly ignore any existing bounds.
#' @inheritParams methodEllipse
#' 
#' @return Modified NMRResonacen1D object.
#' @name initialize_bounds
setGeneric(".initialize_bounds", 
           function(object, overwrite = FALSE, ...) {
             standardGeneric(".initialize_bounds")
           })

#' @rdname initialize_bounds
setMethod(".initialize_bounds", "NMRSpecies1D", 
  function(object, overwrite = FALSE) {

    # Handling bounds if they exist
    bounds <- list(lower = object@bounds$lower, 
                   upper = object@bounds$upper)

    # Selecting default values
    values <- list(lower = -Inf, upper = +Inf)
    columns <- c('position', 'width', 'height', 'fraction.gauss')

    for ( name in names(bounds) ) {
      if ( is.null(bounds[[name]]) || overwrite ) {

        peaks <- object@peaks
        peaks[ , columns] <- values[[name]]

        bounds[[name]] <- peaks
      }
    }

    object@bounds$lower <- bounds$lower
    object@bounds$upper <- bounds$upper

    object
  })



#------------------------------------------------------------------------------
#' Display NMRSpecies1D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRSpecies1D", 
  function(object) {

    # Generating infinite bounds if empty
    object <- .initialize_bounds(object)

    peaks <- object@peaks
    bounds <- object@bounds
    couplings <- object@couplings

    cat('An object of NMRSpecies1D class\n\n')

    # Peaks
    cat('Peaks:\n\n')
    print(peaks)
    cat('\n')

    # Bounds
    columns <- c('position', 'width', 'height', 'fraction.gauss')
    lower <- unlist(bounds$lower[ , columns])
    upper <- unlist(bounds$upper[ , columns])
    
    range <- paste('(', lower, ', ', upper, ')', sep = '')
    peaks[ , columns] <- range

    cat('Bounds (lower, upper):\n\n')
    print(peaks)
    cat('\n')   

    # Couplings
    if ( nrow(couplings) > 0 ) {
      cat('Couplings:\n\n')
      print(couplings)
      cat('\n')
    }
    else {
      cat('No couplings defined.\n')
    }

  })



#==============================================================================>
# Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------------
# Peaks

#' @templateVar slot peaks
#' @template NMRSpecies1D_access
#' @name peaks
#' @export
setGeneric("peaks", 
  function(object, ...) standardGeneric("peaks"))

#' @rdname peaks
#' @export
setMethod("peaks", "NMRSpecies1D", 
  function(object) object@peaks)

#' @templateVar slot peaks
#' @template NMRSpecies1D_replacement
#' @name peaks-set
#' @export
setGeneric("peaks<-", 
  function(object, value) standardGeneric("peaks<-"))

#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRSpecies1D",
  function(object, value) {
    object@peaks <- value
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Couplings

#' @templateVar slot couplings
#' @template NMRSpecies1D_access
#' @name couplings
#' @export
setGeneric("couplings", 
  function(object, ...) standardGeneric("couplings"))

#' @rdname couplings
#' @export
setMethod("couplings", "NMRSpecies1D", 
  function(object) object@couplings)

#' @templateVar slot couplings
#' @template NMRSpecies1D_replacement
#' @name couplings-set
#' @export
setGeneric("couplings<-", 
  function(object, value) standardGeneric("couplings<-"))

#' @rdname couplings-set
#' @export
setReplaceMethod("couplings", "NMRSpecies1D",
  function(object, value) {
    object@couplings <- value
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Bounds

#' @templateVar slot bounds
#' @template NMRSpecies1D_access
#' @name bounds
#' @export
setGeneric("bounds", 
  function(object, ...) standardGeneric("bounds"))

#' @rdname bounds
#' @export
setMethod("bounds", "NMRSpecies1D", 
  function(object) object@bounds)

#' @templateVar slot bounds
#' @template NMRSpecies1D_replacement
#' @name bounds-set
#' @export
setGeneric("bounds<-", 
  function(object, value) standardGeneric("bounds<-"))

#' @rdname bounds-set
#' @export
setReplaceMethod("bounds", "NMRSpecies1D",
  function(object, value) {
    object@bounds <- value
    validObject(object)
    object 
  })



#==============================================================================>
#  Bounds
#==============================================================================>



#------------------------------------------------------------------------------
#' Set general bounds of an NMRSpecies1D object
#' 
#' This function provides a convenience method for generating bounds using a
#' simple set of lower and upper constraints on basic peak parameters such as
#' position, height, width, fraction.gauss. The term "general" refers to the
#' fact that the same bounds are applied to each and every peak, regardless of
#' current parameter values. These bounds can be normalized to a set of data
#' using the optional nmrdata argument.
#' 
#' In practice, general bounds are primarily useful for placing a hard
#' constraint on peak widths and preventing negative heights. Values of 0 for
#' widths and height can also cause issues during optimization, so simple
#' general bounds can be used to prevent errors.
#' 
#' @param object An NMRSpecies1D object.
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
#' @param nmrdata An optional NMRData1D object that can serve as a reference
#'                point for the bounds.
#' @param widen FALSE to prevent new bounds from widening existing bounds.
#' @inheritParams methodEllipse
#' 
#' @return A new NMRSpecies1D object with modified bounds.
#' 
#' @name set_general_bounds
#' @export
setGeneric("set_general_bounds", 
  function(object, position = NULL, height = NULL, width = NULL, 
           fraction.gauss = NULL, nmrdata = NULL, widen = FALSE, ...) {
    standardGeneric("set_general_bounds")
  })

#' @rdname set_general_bounds
#' @export
setMethod("set_general_bounds", "NMRSpecies1D",
  function(object, position = NULL, height = NULL, width = NULL, 
           nmrdata = NULL, widen = FALSE, ...) {

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
    y.range <- range(Re(processed$intensity))
    x.range <- range(processed$direct.shift)

    position <- position * x.range
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

  object
})



#------------------------------------------------------------------------------
#' Set offset bounds of an NMRSpecies1D object
#' 
#' This function provides a convenience method for generating bounds using a
#' simple set of lower and upper constraints on basic peak parameters such as
#' position, height, width, fraction.gauss. The term "offset" refers to the
#' fact that bounds are applied as an offset to current values of the
#' parameters. These bounds can be expressed in absolute (e.g. -0.1 and +0.1
#' ppm) or relative (e.g. -1 percent and +1 percent) terms.
#' 
#' In practice, offset bounds are primarily useful for preventing peak
#' positions from drifting too much from initial guesses and for fine-tuning a
#' fit once an initial optimization is performed. It is not recommended to use
#' strict offset bounds based on rough initial parameter guesses.
#' 
#' @param object An NMRSpecies1D object.
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
#' @return A new NMRSpecies1D object with modified bounds.
#' 
#' @name set_offset_bounds
#' @export
setGeneric("set_offset_bounds", 
  function(object, position = NULL, height = NULL, width = NULL, 
           relative = FALSE, widen = FALSE, ...) {
    standardGeneric("set_offset_bounds")
  })

#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRSpecies1D",
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

  object
})



#------------------------------------------------------------------------------
#' Set conservative bounds on an NMRSpecies1D object
#' 
#' A convenience function that sets reasonable bounds on the fit. These bounds
#' are assumed to be widely applicable to most simple NMR data. Each set of
#' bounds can be turned on or off as necessary. A slightly better set of bounds
#' can be selected if a reference NMRData1D object is provided.
#' 
#' @param object An NMRSpecies object.
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
#' @param nmrdata An optional NMRData1D object that can serve as a reference
#'                point for the bounds.
#' @param widen FALSE to prevent new bounds from widening existing bounds.
#' @inheritParams methodEllipse
#' 
#' @return A new NMRSpecies1D object with modified bounds.
#' 
#' @name set_conservative_bounds
#' @export
setGeneric("set_conservative_bounds", 
  function(object, position = TRUE, height = TRUE, width = TRUE,
           nmrdata = NULL, widen = FALSE, ...) {
    standardGeneric("set_conservative_bounds")
})

#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRSpecies1D",
  function(object, position = TRUE,  height = TRUE, width = TRUE, 
           nmrdata = NULL, widen = FALSE) { 

  # First, do a single pass over general bounds with no reference
  if ( height )  gen.height <- c(0, Inf)
  else gen.height <- NULL

  if ( width ) gen.width <- c(0.003, 3)
  else gen.width <- NULL

  object <- set_general_bounds(object, height = gen.height, width = gen.width)

  # Adding position offsets
  if ( position ) {
    object <- set_offset_bounds(object, position = c(-0.1, 0.1))
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
                                 nmrdata = nmrdata)
  }

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
#' @param object An NMRSpecies1D object.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmrsession_1d(sf = ...), but an override can be
#'           provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single function, FALSE to output a data frame of functions
#'                  that correspond to individual peaks.
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
           components = 'r/i', ...) {
    standardGeneric("f_lineshape")
})

#' @rdname f_lineshape
#' @export
setMethod("f_lineshape", "NMRSpecies1D",
  function(object, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           components = 'r/i') {

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
    parameters <- as.matrix(object@peaks[, columns])

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
      out <- as_tibble(object@peaks[, c('resonance', 'peak')])
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
#' @param object An NMRSpecies1D object.
#' @param direct.shift Vector of chemical shift data in ppm.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmrsession_1d(sf = ...), but an override can be
#'           provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single set of values, FALSE to output a data frame of values
#'                  that correspond to individual peaks.
#' @param components 'r/i' to output both real and imaginary data, 'r' to output
#'                   only real and 'i' to output only imaginary.
#' @inheritParams methodEllipse
#' 
#' @return A vector of spectral intensity data or a data frame with columns
#'         "resonance", "peak", "direct.shift", and "intensity".
#' 
#' @name values
#' @export
setGeneric("values", 
  function(object, direct.shift, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           components = 'r/i', ...) {
    standardGeneric("values")
})

#' @rdname values
#' @export
setMethod("values", "NMRSpecies1D",
  function(object, direct.shift, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           components = 'r/i') {

  # Output depends on whether peaks are summed or not
  if ( sum.peaks ) {
    # Get function
    f <- f_lineshape(object, sf, sum.peaks, components)

    # And apply it to specified chemical shifts
    f(direct.shift)
  } 
  else {
    # Get data frame of functions
    d <- f_lineshape(object, sf, sum.peaks, components)

    # Defining function that generates necessary data frame
    f <- function(g) {
      data.frame(direct.shift = direct.shift, intensity = g[[1]](direct.shift))
    }

    # And apply it for every peak
    d %>%
      group_by(resonance, peak) %>%
      do( f(.$f) )
  }
  })



#------------------------------------------------------------------------
#' Calculate peak areas
#' 
#' Calculate total peak areas based on peak parameters.
#' 
#' @param object An NMRSpecies1D object.
#' @param sf Sweep frequency (in MHz) -- needed to convert peak widths from Hz
#'           to ppm. In most cases, it is recommended to set a single default
#'           value using nmrsession_1d(sf = ...), but an override can be
#'           provided here.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single area, FALSE to output a data frame of peak area
#'                  values.
#' @inheritParams methodEllipse
#' 
#' @return A single overall area or a data frame of areas with columns
#'         "resonance", "peak", and "area".
#' 
#' @name areas
#' @export
setGeneric("areas", 
  function(object, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           components = 'r/i', ...) {
    standardGeneric("areas")
})

#' @rdname areas 
#' @export
setMethod("areas", "NMRSpecies1D",
  function(object, sf = nmrsession_1d('sf'), sum.peaks = TRUE, 
           components = 'r/i') {

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
  peaks <- object@peaks
  areas <- peaks %>%
    group_by(resonance, peak) %>%
    summarize(area = f(position, width, height, fraction.gauss)) %>%
    as.data.frame()

  # Sum if necessary
  if ( sum.peaks ) sum(areas$area)
  else areas
  })
