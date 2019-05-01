# Definition of a class structure for 1D resonance data.

#' @import Rcpp
#' @useDynLib rnmrfit
NULL

#==============================================================================>
#  NMRResonance1D -- peak description 
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR resonance.
#' 
#' Essentially, this class is used to define coupling relationships to group
#' individual peaks into resonances.
#' 
#' @slot peaks A data.frame describing a series of singlets, with one row per
#'             peak. All peaks are characterized by a position (in ppm), height
#'             (in relative intensity units), width (in ppm or Hz), and
#'             fraction.guass (in percent).
#' @slot couplings A data.frame relating the position and area parameters of
#'                 the peaks, effectively combining singlets into multiplets.
#' @slot couplings.leeway A list specifying how tightly enforced the coupling
#'                        constraints on peak positions, areas, and widths
#'                        should be. E.g. position = 0 specifies that the j
#'                        coupling constant is exact, whereas width = 0.1
#'                        specifies that the widths of individual peaks may
#'                        differ by +/- 10 percent.
#' @slot bounds A list of lower and upper bounds on the peak parameters where
#'              both bounds take the same shape as the peaks data.frame.
#' 
#' @name NMRResonance1D-class
#' @export
NMRResonance1D <- setClass("NMRResonance1D",
  slots = c(
    peaks = 'data.frame',
    couplings = 'data.frame',
    couplings.leeway = 'list',
    bounds = 'list'
  ),
  prototype = prototype(
    couplings = data.frame(),
    couplings.leeway = list(position = 0, width = 0, area = 0),
    bounds = list(lower = NULL, upper = NULL)
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRResonance1D validity test
#'
validNMRResonance1D <- function(object) {

  peaks <- object@peaks
  couplings <- object@couplings
  couplings.leeway <- object@couplings.leeway
  bounds <- object@bounds

  valid <- TRUE
  msg <- c()

  #---------------------------------------
  # Checking peak column names
  valid.columns <- c('id', 'peak', 'position', 
                     'width', 'height', 'fraction.gauss')

  if (! identical(colnames(peaks), valid.columns) ) {
    valid <- FALSE
    new.msg <- sprintf('"peaks" must have the following columns: %s',
                       paste(valid.columns, collapse = ', '))
    msg <- c(msg, new.msg)
  }

  #---------------------------------------
  # Checking couplings
  if ( nrow(couplings) > 0 ) {

    valid.columns <- c('id.1', 'id.2', 'peak.1', 'peak.2', 
                       'position.difference', 'area.ratio')
    if (! identical(colnames(couplings), valid.columns) ) {
      valid <- FALSE
      new.msg <- sprintf('"couplings" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      msg <- c(msg, new.msg)
    }

    valid.values <- paste(peaks$id, peaks$peak)
    logic1 <- all(paste(couplings$id.1, couplings$peak.1) %in% valid.values)
    logic2 <- all(paste(couplings$id.2, couplings$peak.2) %in% valid.values)

    if (! (logic1 & logic2) ) {
      valid <- FALSE
      new.msg <- '"couplings" ids and peaks must correspond to existing peaks.'
      msg <- c(msg, new.msg)
    }
  }

  #---------------------------------------
  # Checking that lower bounds match peaks
  if (! is.null(bounds$lower) ) {

    logic <- identical(colnames(bounds$lower), valid.columns)
    if (! logic ) {
      valid <- FALSE
      new.msg <- sprintf('"bounds$lower" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      msg <- c(msg, new.msg)
    }

    logic1 <- identical(bounds$lower$id, peaks$id)
    logic2 <- identical(bounds$lower$peak, peaks$peak)
    if (! (logic1 && logic2) ) {
      valid <- FALSE
      new.msg <- '"bounds$lower" id and peak columns must match "peaks"'
      msg <- c(msg, new.msg)
    }
  }

  #---------------------------------------
  # Checking that upper bounds match peaks
  if (! is.null(bounds$upper) ) {

    logic <- identical(colnames(bounds$upper), valid.columns)
    if (! logic ) {
      valid <- FALSE
      new.msg <- sprintf('"bounds$upper" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      msg <- c(msg, new.msg)
    }

    logic1 <- identical(bounds$upper$id, peaks$id)
    logic2 <- identical(bounds$upper$peak, peaks$peak)
    if (! (logic1 && logic2) ) {
      valid <- FALSE
      new.msg <- '"bounds$upper" id and peak columns must match "peaks"'
      msg <- c(msg, new.msg)
    }
  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else msg
}

# Add the extended validity testing
setValidity("NMRResonance1D", validNMRResonance1D)



#==============================================================================>
#  Helper functions for constructor
#==============================================================================>



#------------------------------------------------------------------------------
#' Parse peak coupling strings
#' 
#' Converts a coupling specification of the form '3.0 d 1.0' into the chemical
#' shift, number of peaks involved and the coupling between them. Currently,
#' the supported codes include s, d, t, q, pentet, sextet, septet, octet, and
#' nonet or any combination of the above. Essentially, the function splits the
#' coupling into a position (number), coupling description (text), and coupling
#' constant (number), parsing the latter two to extract the required codes and
#' numbers. The actual parsing is quite flexible so "1.0 dtpnt 1.5/2/1" will be
#' parsed the same as "1.0 d t pentet 1.5 2 1".
#' 
#' @param coupling.string A character vector with elements of the form '3.0 d
#'                        1.0' or '2 dt 1.0 1.2'
#' 
#' @return A list of three vectors -- chemical.shift, peak numbers, and coupling
#'         constants respectively.
#' 
#' @export
parse_peaks_1d <- function(coupling.string) {

  # Saving original string for later
  original.string <- coupling.string

  # Ensure lowercase
  coupling.string <- tolower(coupling.string)

  # Remove any and all brackets
  coupling.string <- str_replace_all(coupling.string, '[(){}]|\\[\\]', ' ')

  # Replace all possible separators with a single whitespace
  coupling.string <- str_replace_all(coupling.string, '[ ,;#$%_/]+', ' ')

  # Remove whitespace at beginning or end
  coupling.string <- str_trim(coupling.string, 'both')

  # First split directly after first number
  split <- str_split(coupling.string, '(?<=^[0-9.]{1,20})(?![0-9.])', 2)[[1]]
  direct.shift <- split[1]
  direct.shift <- as.numeric(direct.shift)
  coupling.string <- str_trim(split[2], 'both')

  # Then split directly before first remaining number
  split <- str_split(coupling.string, '(?<![0-9.])(?=[0-9])', 2)[[1]]
  codes <- str_trim(split[1], 'both')
  constants <- split[2]

  #---------------------------------------
  # Parsing codes

  # Saving original code text for later
  original <- codes

  # Try every iteration of long names down to three characters
  patterns <- c('pent', 'pnt', 'pentet', 'qui', 'qnt', 'quint', 'quintet',
                'sxt', 'sext', 'sextet', 'spt', 'sept', 'septet', 'hpt',
                'hept', 'heptet', 'oct', 'octet', 'non', 'nonet')
  numbers <- rep(5:9, c(7, 3, 6, 2, 2))

  # Arrange in order of longest to shortest
  index <- order(nchar(patterns), decreasing = TRUE)
  patterns <- patterns[index]
  numbers <- numbers[index]

  replacement <- as.character(numbers)
  names(replacement) <- patterns

  codes <- str_replace_all(codes, replacement)

  # Then parsing single character names
  replacement <- as.character(1:4)
  names(replacement) <- c('s', 'd', 't', 'q')

  codes <- str_replace_all(codes, replacement)

  # Finally, remove any spaces and split by single numbers
  codes <- str_replace_all(codes, '\\s', '')
  codes <- suppressWarnings(as.numeric(str_split(codes, '')[[1]]))

  msg <- 'The coupling definition "%s" could not be parsed.'
  if ( any( is.na(codes) ) ) stop(sprintf(msg, original))

  # If there is only one singlet, no need to parse constants
  if ( identical(codes, 1) ) {
    return(list(direct.shift = direct.shift, numbers = codes, constants = NA))
  }

  #---------------------------------------
  # Parsing constants

  # Saving original constants text for later
  original <- constants

  # Split by single spaces (previous conversions should have removed issues)
  constants <- suppressWarnings(as.numeric(str_split(constants, ' ')[[1]]))

  msg <- 'The coupling constant definition "%s" could not be parsed.'
  if ( any( is.na(constants) ) ) stop(sprintf(msg, original))

  # Making sure that the number of codes matches constants
  msg <- '"%s" was not parsed with an equal number of couplings and constants.'
  logic <- sum(codes != 1) != sum(! is.na(constants))
  if ( logic ) stop(sprintf(msg, original.string))

  # Filling in NA values in case there are singlets
  for (i in which(codes == 1)) {
    if (! is.na(constants[i]) ) {
      constants <- append(constants, NA, i-1)
    }
  }

  list(direct.shift = direct.shift, numbers = codes, constants = constants)
}



#------------------------------------------------------------------------
#' Split peaks data.frame according to specified coupling.
#' 
#' @param peaks A peaks data.frame from NMRResonance1D object.
#' @param number The number of output peaks per input peak.
#' @param constant The coupling constant of the split.
#' 
#' @return A modified data.frame suitable for NMRResonance1D object.
#' 
#' @export
split_peaks_1d <- function(peaks, number, constant) {

  # Singlets do not require splitting
  if ( number == 1 ) return(peaks)

  # Calculating offsets
  offsets <- constant*(number - 1) * seq(-0.5, 0.5, length = number)

  # Calculating height ratios based on pascal's triangle
  ratios <- choose(number - 1, 0:(number - 1))
  heights <- ratios/sum(ratios)

  # Replicating heights and offsets to match data frame
  n <- nrow(peaks)
  offsets <- rep(offsets, n)
  heights <- rep(heights, n)

  # Replicating each row to match number of output peaks
  peaks <- peaks[rep(1:n, each = number), ]
  peaks$height <- peaks$height*heights
  peaks$position <- peaks$position + offsets

  # Resetting peak number and row names
  peaks <- peaks[order(peaks$position), ]
  peaks$peak <- 1:nrow(peaks)
  rownames(peaks) <- peaks$peak

  peaks
}



#------------------------------------------------------------------------
#' Enforce coupling relations between peaks.
#' 
#' Adds a set of constraints between specified peaks that fixes the differences
#' between peak positions and the ratios between peak areas to whatever the
#' current positions and peak heights are. These constraints can be relaxed by
#' setting leeway parameters that convert hard equality constraints to soft
#' inequality constraints around a fraction of the fixed values.
#' 
#' @param nmrresonance An NMRResonance1D object to be modified.
#' @param peaks Vector of peak numbers to enforce coupling. The default is to
#'              select all peaks.
#' 
#' @return An NMRResonance1D object with a new set of coupling constraints.
#' 
#' @export
enforce_couplings_1d <- function(nmrresonance, peaks = NULL) {

  # Checking validity
  if ( class(nmrresonance) != 'NMRResonance1D' ) {
    msg <- '"nmrresonance" must be a valid NMRResonance1D object.'
    stop(msg)
  }
  else {
    validObject(nmrresonance)
  }

  # Filtering peaks
  d <- nmrresonance@peaks
  if ( is.null(peaks) ) peaks <- d$peak
  d <- d[d$peak %in% peaks, ]

  # Calculating
  index.1 <- 1:(nrow(d) - 1)
  index.2 <- 2:nrow(d)

  differences <- d$position[index.2] - d$position[index.1]
  ratios <- d$height[index.2]/d$height[index.1]
  couplings <- data.frame(id.1 = d$id[index.1], id.2 = d$id[index.2], 
                          peak.1 = d$peak[index.1], peak.2 = d$peak[index.2],
                          position.difference = differences,
                          area.ratio = ratios)

  nmrresonance@couplings <- rbind(nmrresonance@couplings, couplings)
  nmrresonance
}




#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRResonance1D object based on simplified peak list
#' 
#' Generates an NMRResonance1D object by converting multiplet definitions into
#' a set of singlets related by constraints on their position and area. Note
#' that peak height parameters are arbitrary at this point. Initial guesses for
#' peak height can be generated during the fit process or overriden manually.
#' 
#' @param peaks A numeric vector of singlet chemical shifts or a character
#'              string specifying multiplets of the form "3 d 1.2". See
#'              ?parse_peaks_1d for more information.
#' @param id A string specifying resonance name. If left empty, a name is
#'           automatically generated from the peaks argument.
#' @param peak.width Initial estimate of peak width (in Hz). For Voigt
#'                   lineshapes, this value is taken as the Lorentzian
#'                   component, with the Gaussian component calculated from
#'                   peak.width*frac.guass/(1-frac.gauss).
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
#' @param area.leeway Similar to position.leeway but for peak areas. 
#'                      Determines how strictly the coupling area ratios are
#'                      enforced.
#' 
#' @return An NMRResonance1D object.
#' 
#' @export
nmrresonance_1d <- function(peaks, id = NULL, peak.width = 1, 
                            fraction.gauss = 0, position.leeway = 0,
                            area.leeway = 0, width.leeway = 0) {

  #---------------------------------------
  # Building peak list

  # Couplings are not added in every case
  add.couplings <- FALSE

  # If peaks is a character, parse coupling information
  if ( is.character(peaks) ) {
    coupling <- parse_peaks_1d(peaks)
    add.constraints <- TRUE

    # Initializing singlet at chemical shift
    if ( is.null(id) ) id <- peaks
    peaks <- data.frame(id = id, peak = 1, position = coupling$direct.shift,
                        width = peak.width, height = 1, 
                        fraction.gauss = fraction.gauss)

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
    if ( is.null(id) ) id <- paste(middle, 'm', range)
    peaks <- data.frame(id = id, peak = 1:length(peaks), position = peaks,
                        width = peak.width, height = 1, 
                        fraction.gauss = fraction.gauss)
  }

  #---------------------------------------
  # Adding coupling definitions

  # Starting with blanks first
  couplings <- data.frame()
  couplings.leeway = list(position = position.leeway, width = width.leeway,
                          area = area.leeway)

  nmrresonance = new('NMRResonance1D', peaks = peaks, couplings = couplings, 
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
#' @param object NMRResonance1D object
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
setMethod(".initialize_bounds", "NMRResonance1D", 
  function(object, overwrite = FALSE) {

    # Handling bounds if they exist
    bounds <- list(lower = object@bounds$lower, 
                   upper = object@bounds$upper)

    # Selecting default values
    values <- list(lower = -Inf, upper = +Inf)
    data.columns <- c('position', 'width', 'height', 'fraction.gauss')

    for ( name in names(bounds) ) {
      if ( is.null(bounds[[name]]) || overwrite ) {

        peaks <- object@peaks
        peaks[ , data.columns] <- values[[name]]

        bounds[[name]] <- peaks
      }
    }

    object@bounds$lower <- bounds$lower
    object@bounds$upper <- bounds$upper

    object
  })



#------------------------------------------------------------------------------
#' Display NMRResonance1D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRResonance1D", 
  function(object) {

    # Generating infinite bounds if empty
    object <- .initialize_bounds(object)

    peaks <- object@peaks
    bounds <- object@bounds
    couplings <- object@couplings

    cat('An object of NMRResonance1D class\n\n')

    # Peaks
    cat('Peaks:\n\n')
    print(peaks)
    cat('\n')

    # Bounds
    data.columns <- c('position', 'width', 'height', 'fraction.gauss')
    lower <- unlist(bounds$lower[ , data.columns])
    upper <- unlist(bounds$upper[ , data.columns])
    
    range <- paste('(', lower, ', ', upper, ')', sep = '')
    peaks[ , data.columns] <- range

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
#' @template NMRResonance1D_access
#' @name peaks
#' @export
setGeneric("peaks", 
  function(object, ...) standardGeneric("peaks"))

#' @rdname peaks
#' @export
setMethod("peaks", "NMRResonance1D", 
  function(object) object@peaks)

#' @templateVar slot peaks
#' @template NMRResonance1D_replacement
#' @name peaks-set
#' @export
setGeneric("peaks<-", 
  function(object, value) standardGeneric("peaks<-"))

#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRResonance1D",
  function(object, value) {
    object@peaks <- value
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Couplings

#' @templateVar slot couplings
#' @template NMRResonance1D_access
#' @name couplings
#' @export
setGeneric("couplings", 
  function(object, ...) standardGeneric("couplings"))

#' @rdname couplings
#' @export
setMethod("couplings", "NMRResonance1D", 
  function(object) object@couplings)

#' @templateVar slot couplings
#' @template NMRResonance1D_replacement
#' @name couplings-set
#' @export
setGeneric("couplings<-", 
  function(object, value) standardGeneric("couplings<-"))

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

#' @templateVar slot bounds
#' @template NMRResonance1D_access
#' @name bounds
#' @export
setGeneric("bounds", 
  function(object, ...) standardGeneric("bounds"))

#' @rdname bounds
#' @export
setMethod("bounds", "NMRResonance1D", 
  function(object) object@bounds)

#' @templateVar slot bounds
#' @template NMRResonance1D_replacement
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
#' using the optional nmrdata argument.
#' 
#' In practice, general bounds are primarily useful for placing a hard
#' constraint on peak widths and preventing negative heights. Values of 0 for
#' widths and height can also cause issues during optimization, so simple
#' general bounds can be used to prevent errors.
#' 
#' @param object An NMRResonance1D object.
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
#' @return A new NMRResonance1D object with modified bounds.
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
setMethod("set_general_bounds", "NMRResonance1D",
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
      msg <- '"nmrdata" must be a valid NMRData1D object.'
      stop(msg)
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
      msg <- paste("All bounds must be vectors of two elements consisting",
                   "of a lower and upper bound.")
      stop(msg)
    }

    if ( bounds[1] > bounds[2] ) {
      msg <- paste("Lower bound must be smaller than upper bound.",
                   "Proceeding with current constraints will result in a",
                   "fit error.")
      warning(msg)
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
#' Set offset bounds of an NMRResonance1D object
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
#' @param object An NMRResonance1D object.
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
#' @return A new NMRResonance1D object with modified bounds.
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
      msg <- paste("All bounds must be vectors of two elements consisting",
                   "of a lower and upper bound.")
      stop(msg)
    }

    msg2 <- "Proceeding with current constraints will result in a fit error."

    if ( bounds[1] > 0 ) {
      msg <- paste("Lower offsets must be negative so that resulting bounds",
                   "include initial values.", msg2)
      stop(msg)
    }

    if ( bounds[2] < 0 ) {
      msg <- paste("Upper offsets must be positive so that resulting bounds",
                   "include initial values.", msg2)
      stop(msg)
    }

    if ( bounds[1] > bounds[2] ) {
      msg <- paste("Lower bound must be smaller than upper bound.", msg2)
      warning(msg)
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
#' Set conservative bounds on an NMRResonance1D object
#' 
#' A convenience function that sets reasonable bounds on the fit. These bounds
#' are assumed to be widely applicable to most simple NMR data. Each set of
#' bounds can be turned on or off as necessary. A slightly better set of bounds
#' can be selected if a reference NMRData1D object is provided.
#' 
#' @param object An NMRResonance object.
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
#' @return A new NMRResonance1D object with modified bounds.
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
setMethod("set_conservative_bounds", "NMRResonance1D",
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
      msg <- '"nmrdata" must be a valid NMRData1D object.'
      stop(msg)
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
#' data frame of functions), where the resulting function(s) output spectral
#' intensity data given a vector input of chemical shifts.
#' 
#' @param object An NMRResonance1D object.
#' @param sum.peaks TRUE to add all individual peaks together and output a
#'                  single function, FALSE to output a data frame of functions
#'                  that correspond to individual peaks.
#' @param components 'r/i' to output both real and imaginary data, 'r' to output
#'                   only real and 'i' to output only imaginary.
#' @inheritParams methodEllipse
#' 
#' @return A function or tbl_df data frame of functions where each function
#'         outputs spectral intensity data given a vector input of chemical
#'         shifts.
#' 
#' @name f_lineshape
#' @export
setGeneric("f_lineshape", 
  function(object, sum.peaks = TRUE, components = 'r/i', ...) {
    standardGeneric("f_lineshape")
})

#' @rdname f_lineshape
#' @export
setMethod("f_lineshape", "NMRResonance1D",
  function(object, sum.peaks = TRUE, components = 'r/i') {

  # Converting object peak_units to ppm
  object <- set_peak_units(object, peak.units = 'ppm')

  # Performing calculation based on peak.type
  peak.type <- object@peak_type

  # Extracting values
  peaks <- as.matrix(object@peaks[, -(1:2)])

  # Define a function for a single set of peakseters
  if ( peak.type == 'lorenz' ) {
    f <- function(i, x) {
      p <- peaks[i, 1]
      h <- peaks[i, 2]
      w <- peaks[i, 3]
      z <- (x - p)/w
      
      h*complex(re = 1, im = z)/(1 + z^2)
    }
  } else if ( peak.type == 'gauss' ) {
    f <- function(i, x) {
      p <- peaks[i, 1]
      h <- peaks[i, 2]
      w <- peaks[i, 3]
      z <- (x - p)/(sqrt(2)*w)
      
      h*Faddeeva_w(z)
    }
  } else if ( peak.type == 'pvoigt' ) {
    f <- function(i, x) {
      p <- peaks[i, 1]
      lh <- peaks[i, 2]
      gh <- peaks[i, 3]
      w <- peaks[i, 4]
      z <- (x - p)/w
      
      lh*complex(re = 1, im = z)/(1 + z^2) + gh*Faddeeva_w(z)
    }
  } else if ( peak.type == 'voigt' ) {
    f <- function(i, x) {
      p <- peaks[i, 1]
      h <- peaks[i, 2]
      lw <- peaks[i, 3]
      gw <- peaks[i, 4]
      z <- (x - p + complex(im = lw))/(sqrt(2)*gw)
      
      h*Faddeeva_w(z)/Faddeeva_w(complex(im = lw)/(sqrt(2)*gw))
    }
  }

  # Then wrap the internal function with the ability to apply
  # over multiple peaks.
  f_outer <- function(x, i, include.convolution = TRUE) {

    # If convolution is required, the weights may need to be adjusted to the
    # specified x values.
    nmrdata <- object@nmrdata
    logic <- include.convolution && 
             (! is.null(nmrdata) ) &&
             (length(nmrdata@convolution) > 0)

    if ( logic ) {

      convolution <- nmrdata@convolution
      processed <- nmrdata@processed

      dx <- abs(mean(diff(processed$direct.shift)))
      new.dx <-  sort(unique(diff(x)))

      if ( (length(new.dx) > 1) && any(abs(diff(new.dx)) > 1e-10) ) {
        msg <- paste('Convolution vectors can only be considered for evenly',
                     'sampled data.')
        stop(msg, call. = TRUE)
      } else {
        new.dx <- abs(mean(new.dx))
      }

      new.n <- round(length(convolution)*new.dx/dx)
      convolution <- spline(1:length(convolution), convolution, n = new.n)$y

      f_inner <- function(i, x) {
        out <- convolve(f(i, x), convolution, type = 'open')
        n.edge <- (length(convolution) - 1)/2
        index <- (1 + n.edge):(length(x) + n.edge)
        out[index]
      }
    } else {
      f_inner <- f
    }

    out.list <- lapply(i, f_inner, x = x)

    n.peaks <- length(i)
    n.points <- length(x)
    
    intensity <- unlist(out.list)
    direct.shift <- rep(x, n.peaks)

    cbind(object@peaks[rep(i, each = n.points), 1:2],
          direct.shift = direct.shift,
          intensity = intensity)
  }

  # Return the function
  f_outer
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
#' @rdname calc_lineshape
#' @export
setMethod("calc_lineshape", "NMRScaffold1D",
  function(object, direct.shift = NULL, include.convolution = TRUE) {

  # If direct.shift is NULL check for available data
  if ( is.null(direct.shift) ) {
    if (! is.null(object@nmrdata) ) {
      direct.shift <- object@nmrdata@processed$direct.shift
    } else {
      msg <- '"direct.shift" must be provided if "nmrdata" slot is NULL.'
      stop(msg)
    }
  }

  # Generating function
  f <- f_lineshape(object)

  # Calculate all lineshapes
  n <- nrow(object@peaks)

  # Apply it
  f(direct.shift, 1:n, include.convolution)
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

#------------------------------------------------------------------------
#' @rdname calc_area
#' @export
setMethod("calc_area", "NMRScaffold1D",
  function(object, type = 'analytical') {

  # Performing calculation based on peak.type
  peak.type <- object@peak_type
  d <- object@peaks

  if ( type == 'analytical' ) {

    if ( peak.type == 'lorenz' ) {
      d <- mutate(d, area = pi*width*height)
    } else if ( peak.type == 'gauss' ) {
      d <- mutate(d, area = sqrt(2*pi)*width*height)
    } else if ( peak.type == 'pvoigt' ) {
      d <- mutate(d, area = width*(pi*l.height + sqrt(pi)*g.height))
    } else if ( peak.type == 'voigt' ) {
      d <- mutate(d, area = Re(sqrt(2*pi)*g.width*height/
                            Faddeeva_w(complex(im = l.width)/(sqrt(2)*g.width))))
    }

  } else if ( type == 'numerical' ) {

    # Generate lineshape function
    f <- f_lineshape(object, FALSE)

    # Determine how many rows there are
    n <- nrow(d)

    f_area <- function(i) {
      integrate(function(x) Re(f(x)), -Inf, Inf, i = i)$value
    }

    areas <- lapply(1:n, f_area)
    d$area <- unlist(areas)
  }

  # Dropping parameter columns and returning
  select(d, id, peak, area)
})




