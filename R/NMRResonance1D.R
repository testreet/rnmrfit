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
#' @slot id A name to be used for the object in output data.
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
    id = 'character',
    peaks = 'data.frame',
    couplings = 'data.frame',
    couplings.leeway = 'list',
    bounds = 'list'
  ),
  prototype = prototype(
    id = 'resonance',
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

  id <- object@id
  peaks <- object@peaks
  couplings <- object@couplings
  couplings.leeway <- object@couplings.leeway
  bounds <- object@bounds

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Checking name
  if ( length(id) != 1 ) {
    valid <- FALSE
    new.err <- '"name" must be a character vector of length 1.'
    err <- c(err, new.err)
  }

  #---------------------------------------
  # Checking peak column names
  valid.columns <- c('peak', 'position', 'width', 'height', 'fraction.gauss')

  if (! identical(colnames(peaks), valid.columns) ) {
    valid <- FALSE
    new.err <- sprintf('"peaks" must have the following columns: %s',
                       paste(valid.columns, collapse = ', '))
    err <- c(err, new.err)
  }

  #---------------------------------------
  # Checking that lower bounds match peaks
  if (! is.null(bounds$lower) ) {

    logic <- identical(colnames(bounds$lower), valid.columns)
    if (! logic ) {
      valid <- FALSE
      new.err <- sprintf('"bounds$lower" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      err <- c(err, new.err)
    }

    logic1 <- identical(bounds$lower$resonance, peaks$resonance)
    logic2 <- identical(bounds$lower$peak, peaks$peak)
    if (! (logic1 && logic2) ) {
      valid <- FALSE
      new.err <- '"bounds$lower" resonance and peak columns must match "peaks"'
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking that upper bounds match peaks
  if (! is.null(bounds$upper) ) {

    logic <- identical(colnames(bounds$upper), valid.columns)
    if (! logic ) {
      valid <- FALSE
      new.err <- sprintf('"bounds$upper" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      err <- c(err, new.err)
    }

    logic1 <- identical(bounds$upper$resonance, peaks$resonance)
    logic2 <- identical(bounds$upper$peak, peaks$peak)
    if (! (logic1 && logic2) ) {
      valid <- FALSE
      new.err <- '"bounds$upper" resonance and peak columns must match "peaks"'
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking couplings
  if ( nrow(couplings) > 0 ) {

    valid.columns <- c('peak.1', 'peak.2', 'position.difference', 'area.ratio')
    if (! identical(colnames(couplings), valid.columns) ) {
      valid <- FALSE
      new.err <- sprintf('"couplings" must have the following columns: %s',
                         paste(valid.columns, collapse = ', '))
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
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

  err <- 'The coupling definition "%s" could not be parsed.'
  if ( any( is.na(codes) ) ) stop(sprintf(err, original))

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

  err <- 'The coupling constant definition "%s" could not be parsed.'
  if ( any( is.na(constants) ) ) stop(sprintf(err, original))

  # Making sure that the number of codes matches constants
  err <- '"%s" was not parsed with an equal number of couplings and constants.'
  logic <- sum(codes != 1) != sum(! is.na(constants))
  if ( logic ) stop(sprintf(err, original.string))

  # Filling in NA values in case there are singlets
  for (i in which(codes == 1)) {
    if (! is.na(constants[i]) ) {
      constants <- append(constants, NA, i-1)
    }
  }

  list(direct.shift = direct.shift, numbers = codes, constants = constants)
}



#------------------------------------------------------------------------------
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



#------------------------------------------------------------------------------
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
    err <- '"nmrresonance" must be a valid NMRResonance1D object.'
    stop(err)
  }
  else {
    validObject(nmrresonance)
  }

  # Filtering peaks
  d <- nmrresonance@peaks
  if ( is.null(peaks) ) peaks <- d$peak
  d <- d[d$peak %in% peaks, ]

  # There must be more than one peak to add couplings
  if ( nrow(d) <= 1 ) return(nmrresonance)

  # Calculating
  index.1 <- 1:(nrow(d) - 1)
  index.2 <- 2:nrow(d)

  differences <- d$position[index.2] - d$position[index.1]
  ratios <- d$height[index.2]/d$height[index.1]
  couplings <- data.frame(peak.1 = d$peak[index.1], peak.2 = d$peak[index.2],
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
#' @param sf Sweep frequency (in MHz) -- needed to convert coupling constants
#'           from Hz to ppm. In most cases, it is recommended to set a single
#'           default value using nmrsession_1d(sf = ...), but an override can be
#'           provided here.
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
#' @return An NMRResonance1D object.
#' 
#' @export
nmrresonance_1d <- function(peaks, sf = nmrsession_1d('sf'), id = NULL, 
                            width = 1, fraction.gauss = 0, 
                            position.leeway = 0, area.leeway = 0, 
                            width.leeway = 0) {

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
    peaks <- data.frame(peak = 1, position = coupling$direct.shift,
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
    if ( is.null(id) ) id <- paste(middle, 'm', range)
    peaks <- data.frame(peak = 1:length(peaks), 
                        position = peaks, width = width, height = 1, 
                        fraction.gauss = fraction.gauss)
  }

  #---------------------------------------
  # Adding coupling definitions

  # Starting with blanks first
  couplings <- data.frame()
  couplings.leeway = list(position = position.leeway, width = width.leeway,
                          area = area.leeway)

  nmrresonance = new('NMRResonance1D', id = id, peaks = peaks, 
                                       couplings = couplings, 
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
#' Display NMRResonance1D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRResonance1D", 
  function(object) {

    # Generating infinite bounds if empty
    object <- .initialize_bounds(object)

    id <- object@id
    peaks <- object@peaks
    bounds <- object@bounds
    couplings <- object@couplings

    cat('An object of NMRResonance1D class\n\n')

    # ID 
    cat('Id: ')
    cat(id)
    cat('\n\n')

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
# Id

#---------------------------------------
#' Get object id
#' 
#' Generic convenience method to access the id of an NMRResonance1D or
#' NMRSpecies1D object.
#' 
#' @param object An NMRResonance1D or NMRSpecies1D object.
#' @param ... Additional arguments passed to inheriting methods.
#'
#' @name id
#' @export
setGeneric("id", 
  function(object, ...) standardGeneric("id")
  )

#' @rdname id
#' @export
setMethod("id", "NMRResonance1D", 
  function(object) object@id)

#---------------------------------------
#' Set object id
#' 
#' Generic convenience method to set the id of an NMRResonance1D or
#' NMRSpecies1D object.
#' 
#' @param object An NMRResonance1D or NMRSpecies1D object.
#' @param value New id (converted to character).
#'
#' @name id-set
#' @export
setGeneric("id<-", 
  function(object, value) standardGeneric("id<-")
  )

#' @rdname id-set
#' @export
setReplaceMethod("id", "NMRResonance1D",
  function(object, value) {
    object@id <- as.character(value)
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Peaks

#---------------------------------------
#' Get object peaks
#' 
#' Generic convenience method to access the peak definitions of an
#' NMRResonance1D object. If used on an NMRSpecies1D or NMRFit1D object, the
#' function combines the peaks data frames of component resonances.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param include.id TRUE to return a column of resonance or species ids.
#' @param ... Additional arguments passed to inheriting methods.
#' 
#' @name peaks
#' @export
setGeneric("peaks", 
  function(object, include.id = FALSE, ...) standardGeneric("peaks")
  )

#' @rdname peaks
#' @export
setMethod("peaks", "NMRResonance1D", 
  function(object, include.id = FALSE) {
    if ( include.id ) cbind(resonance = object@id, object@peaks)
    else object@peaks
  })

#---------------------------------------
#' Set object peaks
#' 
#' Generic convenience method to set the peak definitions of an NMRResonance1D
#' NMRSpecies1D, or NMRFit1D object. This is primarily intended as an internal
#' method, so use with caution. Changing peak definitions after an object has
#' been defined may have unpredictable consequences.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param value A data frame with "position", "width", "height", and
#'              "fraction.gauss" columns. Peaks may be defined by one to three
#'              columns of "peak", "resonance", and "species" depending on the
#'              nature of the original object.
#' 
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
    if ( include.id ) {
      bounds$lower <- cbind(resonance = object@id, bounds$lower)
      bounds$upper <- cbind(resonance = object@id, bounds$upper)
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
#' Update NMRResonance1D peak parameters
#' 
#' This is an internal function whose role is to update existing peak
#' parameters, while accounting for exclusion criteria and generating relevant
#' errors/messages.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param peaks A data frame with "position", "width", "height", and
#'              "fraction.gauss" columns. Peaks may be defined by one to three
#'              columns of "peak", "resonance", and "species" depending on the
#'              nature of the original object.
#' @param exclusion.level A string specifying what to do when peaks are found to
#'                        fall outside of the data range: either 'species' to
#'                        exclude the whole species to which the offending peak
#'                        belongs, 'resonance' to exclude the resonance to which
#'                        the offending peak belongs, or 'peak' to exclude just
#'                        the peak itself.
#' @param exclusion.action A function specifying what to do when peaks are found
#'                         to be outside the data range. The default action is
#'                         to generate a warning message. The function must take
#'                         two inputs, a string with a notification message and
#'                         a data.frame of all the excluded peaks.
#' @inheritParams methodEllipse
#' 
#' @return A new object with modified peak parameters.
#' 
#' @name update_peaks
setGeneric("update_peaks", 
  function(object, ...) {
    standardGeneric("update_peaks")
  })

#' @rdname update_peaks
#' @export
setMethod("update_peaks", "NMRResonance1D",
  function(object, peaks, exclusion.level = nmrsession_1d$exclusion$level,
           exclusion.action = nmrsession_1d$exclusion$action) {

  # Check that columns match before continuing
  current.peaks <- peaks(object)
  err <- '"peaks" columns must match those of current peaks data.frame.'
  if (! all(colnames(peaks) %in% colnames(current.peaks))) stop(err)

  # Isolating ids
  data.columns <- c('position', 'width', 'height', 'fraction.gauss')
  all.columns <- colnames(current.peaks)
  id.columns <- all.columns[! all.columns %in% data.columns]

  # First, check if there are any generic ids missing
  current.ids <- apply(current.peaks[, id.columns], 1, paste, collapse = '-')
  new.ids <- apply(peaks[, id.columns], 1, paste, collapse = '-')
  logic <- ! current.ids %in% new.ids

  if ( any(logic) ) {

    msg <- 'Some peaks were found outside data range and were not updated.'
    
    # There are three conditions under which the logic vector will expand
    # to exclude more terms: exclusion level is species and species are present,
    # exclusion level is species but only resonances are present, 
    # exclusion level is resonance and resonances are present, but in this third
    # case, the print message needs to change to reflect species

    if ( exclusion.level == 'species' ) {

      if ( 'species' %in% id.columns ) {
        species <- unique(current.peaks$species)
        logic <- logic | ( current.peaks$species %in% species ) 
        msg <- sprintf('%s\nThe following species were excluded entirely: %s',
                       msg, paste(species, collapse = ', '))
      }
      else if ( 'resonance' %in% id.columns ) {
        resonances <- unique(current.peaks$resonance)
        logic <- logic | ( current.peaks$resonances %in% resonances ) 
        msg <- sprintf('%s\nThe following resonances were excluded entirely: %s',
                       msg, paste(resonances, collapse = ', '))
      }
    }
    else if ( exclusion.level == 'resonance' ) {

      # In this case, the ids are built up iteratively
      if ( 'resonance' %in% id.columns ) {
        ids <- current.peaks$resonance

        if ( 'species' %in% id.columns ) {
          ids <- paste(current.peaks$species, ids, paste = '-')
        }

        pattern <- paste('^', ids, sep = '')
        logic <- logic | str_detect(current.ids, pattern)
        msg <- sprintf('%s\nThe following resonances were excluded entirely: %s',
                       msg, paste(unique(ids[logic]), collapse = ', '))
      }
    }
  }
  else {

    ids <- current.peaks$peak

    if ( 'resonance' %in% id.columns ) {
      ids <- paste(current.peaks$resonance, ids, paste = '-')
    }

    if ( 'species' %in% id.columns ) {
      ids <- paste(current.peaks$species, ids, paste = '-')
    }

    msg <- sprintf('%s\nThe following peaks were excluded: %s',
                     msg, paste(unique(ids[logic]), collapse = ', '))
  }

})

#------------------------------------------------------------------------------
#' Initialize peak heights of an NMRResonance1D object
#' 
#' Generates peak height estimates based on spectral data. If applied to an
#' NMRSpecies1D or NMRFit1D object, the same initialization is propagated to
#' every component resonance.
#' 
#' At this point, there is just one approach: take peak height as the intensity
#' of the data at the current position of the peak. There are plans to develop
#' more sophisticated approaches in the future.
#' 
#' @param object An NMRResonance1D, NMRSpecies1D, or NMRFit1D object.
#' @param nmrdata An NMRData1D object.
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
  function(object, nmrdata) {

  # Checking nmrdata
  if ( class(nmrdata) != 'NMRData1D' ) {
    err <- '"nmrdata" must be a valid NMRData1D object.'
    stop(err)
  }
  else {
    validObject(nmrdata)
  }

  # Building a simple interpolating function betwen chemical shift and intensity
  d <- processed(nmrdata)
  f <- approxfun(d$direct.shift, d$intensity)

  # Checking that all peak positions are inside the provided data
  p <- object@peaks
  logic <- (p$position > min(d$direct.shift)) & 
           (p$position < max(d$direct.shift))


  # Generating heights from interpolation

  d$height <- f(d$position)


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
      if ( include.id ) {
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
