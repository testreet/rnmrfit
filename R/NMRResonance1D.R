# Definition of a class structure for 1D resonance data.



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
#' @slot couplings A data.frame relating the position and height parameters of
#'                 the peaks, effectively combining singlets into multiplets.
#' @slot couplings.leeway A list specifying how tightly enforced the coupling
#'                        constraints on peak positions, heights, and widths
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
    couplings.leeway = list(position = 0, width = 0, height = 0),
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
                       'position.difference', 'height.ratio')
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
#' @param height.leeway Similar to position.leeway but for peak heights.
#'                      Determines how strictly the coupling height ratios are
#'                      enforced.
#' 
#' @return An NMRScaffold1D object.
#' 
#' @export
nmrresonance_1d <- function(peaks, id = NULL, peak.width = 1, 
                            fraction.gauss = 0, position.leeway = 0,
                            height.leeway = 0, width.leeway = 0) {

  #---------------------------------------
  # Building peak list

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

  # Generate rough estimates for height
  #d <- nmrdata@processed
  #logic <- which_approx(d$direct.shift, peaks$position)
  #height <- Re(d$intensity)[logic]
  #height <- ifelse(height < 0, 0.1*max(Re(d$intensity)), height)

  couplings <- data.frame()
  couplings.leeway = list(position = position.leeway, width = width.leeway,
                          height = height.leeway)

  new('NMRResonance1D', peaks = peaks, couplings = couplings, 
                        couplings.leeway = couplings.leeway)

}



#========================================================================>
# Display function
#========================================================================>



#------------------------------------------------------------------------
#' Display NMRScaffold1D object
#'
#' Display a quick summary of scaffold parameters.
#'
#' @export
setMethod("show", "NMRScaffold1D", 
  function(object) {

    peaks <- object@peaks
    baseline <- object@baseline
    baseline.diff <- object@baseline_difference
    knots <- object@knots
    phase <- object@phase
    peak.type <- object@peak_type
    peak.units <- object@peak_units
    normalized <- object@normalized
    constraints <- object@constraints
    nmrdata <- object@nmrdata

    cat('An object of NMRScaffold1D class\n\n')
    cat(sprintf('Peak type: %s\n', peak.type))  
    cat(sprintf('Peak width units: %s\n', peak.units))  
    cat(sprintf('Normalized parameters: %s\n\n', normalized))  

    cat('Peaks:\n')
    print(peaks)
    cat('\n')

    if ( length(baseline) > 0 ) {
      cat('Baseline:\n')
      print(baseline)

      if ( length(baseline.diff) > 0 ) {
        cat('Baseline difference:\n')
        print(baseline.diff)
      }

      cat('Knots:\n')
      print(knots)
      cat('\n')
    }

    if ( length(phase) > 0 ) {
      cat('Phase:\n')
      print(phase)
      cat('\n')
    }

    if ( nrow(constraints) > 0 ) {
      cat('Constraints:\n')
      print(constraints)
      cat('\n')
    }

    if (! is.null(nmrdata) ) {
      cat('Data:\n')
      print(summary(nmrdata))
      cat('\n')
    }
  })



#========================================================================>
# Helper functions for determining column names 
#========================================================================>



#------------------------------------------------------------------------
#' @rdname data_columns
setMethod(".data_columns", "NMRScaffold1D",
  function(object, index = FALSE, peak.type = NULL) {
    if ( is.null(peak.type) ) peak.type <- object@peak_type
    column.names <- list('lorenz' = c('position', 'height', 'width'),
                         'gauss' = c('position', 'height', 'width'),
                         'pvoigt' = c('position', 'l.height', 
                                      'g.height', 'width'),
                         'voigt' = c('position', 'height', 
                                     'l.width', 'g.width'))
    column.indexes <- list('lorenz' = 3:5, 'gauss' = 3:5,
                           'pvoigt' = 3:6, 'voigt' = 3:6)

    if ( index ) column.indexes[[peak.type]]
    else column.names[[peak.type]]
  })

#------------------------------------------------------------------------
#' @rdname all_columns
setMethod(".all_columns", "NMRScaffold1D",
  function(object, index = FALSE, peak.type = NULL) {
    if ( is.null(peak.type) ) peak.type <- object@peak_type
    if ( index ) c(1:2, .data_columns(object, TRUE, peak.type))
    else c('id', 'peak', .data_columns(object, FALSE, peak.type))
  })

#------------------------------------------------------------------------
#' @rdname position_columns
setMethod(".position_columns", "NMRScaffold1D",
  function(object, index = FALSE, peak.type = NULL) {
    if ( is.null(peak.type) ) peak.type <- object@peak_type
    if ( index ) 3
    else 'position'
  })

#------------------------------------------------------------------------
#' @rdname height_columns
setMethod(".height_columns", "NMRScaffold1D",
  function(object, index = FALSE, peak.type = NULL) {
    if ( is.null(peak.type) ) peak.type <- object@peak_type
    column.names <- list('lorenz' = 'height', 'gauss' = 'height',
                         'pvoigt' = c('l.height', 'g.height'), 
                         'voigt' = 'height')
    column.indexes <- list('lorenz' = 4, 'gauss' = 4,
                           'pvoigt' = 4:5, 'voigt' = 4)
    if ( index ) column.indexes[[peak.type]]
    else column.names[[peak.type]]
  })

#------------------------------------------------------------------------
#' @rdname width_columns
setMethod(".width_columns", "NMRScaffold1D",
  function(object, index = FALSE, peak.type = NULL) {
    if ( is.null(peak.type) ) peak.type <- object@peak_type
    column.names <- list('lorenz' = 'width', 'gauss' = 'width',
                         'pvoigt' = 'width', 'voigt' = c('l.width', 'g.width'))
    column.indexes <- list('lorenz' = 5, 'gauss' = 5,
                           'pvoigt' = 6, 'voigt' = 5:6)
    if ( index ) column.indexes[[peak.type]]
    else column.names[[peak.type]]
  })



#========================================================================>
# Internal functions dealing with flattened parameters
#========================================================================>



#------------------------------------------------------------------------
#' @rdname gen_parameters
setMethod(".gen_parameters", "NMRScaffold1D",
  function(object, peaks, baseline, baseline.diff, phase) {

    peak.names <- paste(rep(colnames(peaks)[-(1:2)], nrow(peaks)),
                        rep(1:nrow(peaks), each = ncol(peaks) - 2),   
                        sep = '.')

    # The rest of the parameters are handled together
    combined <-  list(baseline, baseline.diff, phase)
    combined.names <- c('baseline', 'baseline.diff', 'phase')
    combined.lengths <- vapply(combined, length, 0)

    # Condensing and returning
    parameters <- c(as.numeric(t(peaks[,-(1:2)])), unlist(combined))
    names(parameters) <- c(peak.names, rep(combined.names, combined.lengths)) 

    parameters
  })

#------------------------------------------------------------------------
#' @rdname merge_parameters
setMethod(".merge_parameters", "NMRScaffold1D",
  function(object, return.object) {

    # Getting values
    peaks <- object@peaks
    baseline <- object@baseline
    baseline.diff <- object@baseline_difference
    phase <- object@phase

    parameters <- .gen_parameters(object, peaks, baseline, baseline.diff, phase)

    if ( return.object ) {
      object@parameters <- parameters
      object
    } else {
      parameters
    }
  })

#------------------------------------------------------------------------
#' @rdname spread_parameters
setMethod(".spread_parameters", "NMRScaffold1D",
  function(object) {

    # Getting values
    peaks <- object@peaks
    baseline <- object@baseline
    baseline.diff <- object@baseline_difference
    phase <- object@phase

    parameters <- object@parameters

    # Determining vector breakdown
    nc <- as.integer(nrow(peaks)*(ncol(peaks)-2))
    nb1 <- as.integer(length(baseline))
    nb2 <- as.integer(length(baseline.diff))
    np <- as.integer(length(phase))

    # Splitting up the p results into curve, baseline, and phase components
    fit.data <- t(matrix(parameters[1:nc], nrow = (ncol(peaks) - 2)))
    peaks[, -(1:2)] <- as.data.frame(fit.data)
    object@peaks <-peaks 
    parameters <- parameters[-(1:nc)]

    # Extracting baseline
    if (nb1 > 0) {
      object@baseline <- parameters[1:nb1]
      parameters <- parameters[-(1:nb1)]
    }

    # Extracting  baseline difference
    if ( nb2 > 0 ) {
      object@baseline_difference <- parameters[1:nb2]
      parameters <- parameters[-(1:nb2)]
    }

    # Extracting phase
    if ( np > 0 ) {
      object@phase <- parameters[1:np]
    }

   object 
  })



#========================================================================>
# Internal functions relating to bounds
#========================================================================>



#' @rdname propagate_to_bounds
setMethod(".propagate_to_bounds", "NMRScaffold1D", 
  function(object, f, ...) callNextMethod() )

#' @rdname drop_bounds
setMethod(".drop_bounds", "NMRScaffold1D", 
  function(object, f, ...) callNextMethod())

#' @rdname is_conformant
setMethod(".is_conformant", "NMRScaffold1D", 
  function(object1, object2) {

  # If object2 is NULL, then it's conformant
  if ( is.null(object2) ) return(TRUE)

  # If object1 and object2 aren't the same class, they're not conformant
  if ( class(object1) != class(object2) ) return(FALSE)

  # Running through all the parameters
  if ( nrow(object1@peaks) != nrow(object2@peaks) ||
       ncol(object1@peaks) != ncol(object2@peaks) ||
       length(object1@baseline) != length(object2@baseline) ||
       length(object1@baseline_difference) != 
       length(object2@baseline_difference) ||
       length(object1@phase) != length(object2@phase) ||
       object1@peak_type != object2@peak_type ) FALSE
  else TRUE
  })



#========================================================================>
# Misc internal functions
#========================================================================>



#' @rdname check_data
setMethod(".check_data", c("NMRScaffold1D", "ANY"),
  function(object, nmrdata) {
    callNextMethod()
  })



#========================================================================>
# Basic setter and getter functions
#========================================================================>



#------------------------------------------------------------------------
#' @rdname peaks
#' @export
setMethod("peaks", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRScaffold1D",
                 function(object, value) callNextMethod())

#------------------------------------------------------------------------
#' @templateVar slot baseline
#' @template NMRScaffold_access
#' @name baseline
#' @export
setGeneric("baseline", function(object, ...) 
           standardGeneric("baseline"))

#' @rdname baseline
#' @export
setMethod("baseline", "NMRScaffold1D", 
          function(object) object@baseline)

#' @templateVar slot baseline
#' @template NMRScaffold_replacement
#' @name baseline-set
#' @export
setGeneric("baseline<-", 
           function(object, value) standardGeneric("baseline<-"))

#' @rdname baseline-set
#' @export
setReplaceMethod("baseline", "NMRScaffold1D",
                 function(object, value) {
                   if ( is.null(value) ) value <- numeric(0)
                   object@baseline <- value
                   object <- .merge_parameters(object)
                   object <- .drop_bounds(object)
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @templateVar slot baseline_difference 
#' @template NMRScaffold_access
#' @name baseline_difference 
#' @export
setGeneric("baseline_difference", function(object, ...) 
           standardGeneric("baseline_difference"))

#' @rdname baseline_difference 
#' @export
setMethod("baseline_difference", "NMRScaffold1D", 
          function(object) object@baseline_difference)

#' @templateVar slot baseline_difference 
#' @template NMRScaffold_replacement
#' @name baseline_difference-set
#' @export
setGeneric("baseline_difference<-", 
           function(object, value) standardGeneric("baseline_difference<-"))

#' @rdname baseline_difference-set
#' @export
setReplaceMethod("baseline_difference", "NMRScaffold1D",
                 function(object, value) {
                   if ( is.null(value) ) value <- numeric(0)
                   object@baseline_difference <- value
                   object <- .merge_parameters(object)
                   object <- .drop_bounds(object)
                   validObject(object)
                   object
                 })


#------------------------------------------------------------------------
#' @templateVar slot knots 
#' @template NMRScaffold_access
#' @name baseline_knots 
#' @export
setGeneric("baseline_knots", function(object, ...) 
           standardGeneric("baseline_knots"))

#' @rdname baseline_knots 
#' @export
setMethod("baseline_knots", "NMRScaffold1D", 
          function(object) object@knots)

#' @templateVar slot knots
#' @template NMRScaffold_replacement
#' @name baseline_knots-set
#' @export
setGeneric("baseline_knots<-", 
           function(object, value) standardGeneric("baseline_knots<-"))

#' @rdname baseline_knots-set
#' @export
setReplaceMethod("baseline_knots", "NMRScaffold1D",
                 function(object, value) {
                   if ( is.null(value) ) value <- numeric(0)
                   
                   if ( length(value) != length(object@knots) ) {
                     msg <- paste('Baseline parameters no longer correspond to',
                                  'knot number, setting baseline to zero.')
                     warning(msg, call. = FALSE)
                     object@baseline <- rep(0, length(object@baseline)) 
                     object <- .merge_parameters(object)
                     object <- .drop_bounds(object)
                   }

                   object@knots <- value
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @rdname phase
#' @export
setMethod("phase", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname phase-set
#' @export
setReplaceMethod("phase", "NMRScaffold1D",
                 function(object, value) callNextMethod())

#------------------------------------------------------------------------
#' @rdname constraints
#' @export
setMethod("constraints", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname constraints-set
#' @export
setReplaceMethod("constraints", "NMRScaffold1D",
                 function(object, value) callNextMethod())

#------------------------------------------------------------------------
#' @rdname bounds
#' @export
setMethod("bounds", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname bounds-set
#' @export
setReplaceMethod("bounds", "NMRScaffold1D",
                 function(object, value) callNextMethod())

#------------------------------------------------------------------------
#' @rdname lower_bounds
#' @export
setMethod("lower_bounds", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname lower_bounds-set
#' @export
setReplaceMethod("lower_bounds", "NMRScaffold1D",
                 function(object, value) callNextMethod())

#------------------------------------------------------------------------
#' @rdname upper_bounds
#' @export
setMethod("upper_bounds", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname upper_bounds-set
#' @export
setReplaceMethod("upper_bounds", "NMRScaffold1D",
                 function(object, value) callNextMethod())

#------------------------------------------------------------------------
#' @rdname nmrdata
#' @export
setMethod("nmrdata", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname nmrdata-set
#' @export
setReplaceMethod("nmrdata", "NMRScaffold1D",
                 function(object, value) callNextMethod())



#========================================================================>
# Compound setter and getter functions (requiring internal calculations)
#========================================================================>



#------------------------------------------------------------------------
#' @rdname peak_type
#' @export
setMethod("peak_type", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname peak_type-set
#' @export
setReplaceMethod("peak_type", "NMRScaffold1D",
                 function(object, value) callNextMethod())

#------------------------------------------------------------------------
#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRScaffold1D", 
  function(object, peak.type, frac.lorenz = 0.9) {

  valid.types <- c('lorenz', 'gauss', 'pvoigt', 'voigt')
  if (! peak.type %in% valid.types ) {
    msg <-sprintf('"peak_type" must be one of %s', 
                  paste(valid.types, collapse = ', '))
    stop(msg, call. = FALSE)
  }

  current.type <- object@peak_type
  
  # Do nothing if the scaffold is already associated with correct peak.type
  if ( current.type == peak.type ) {
    return(object)
  }

  # All of the following operations will center on the peaks data.frame
  d <- object@peaks
  d$area <- calc_area(object)$area

  # Final columns will be the same for next three cases
  columns <- c(.all_columns(object, peak.type = 'lorenz'), 'area')


  # Conversion to lorenz if necessary
  if ( current.type == 'gauss' ) {

    d$width <- d$area/(pi*d$height)
    d <- d[, columns]

  } else if ( current.type == 'pvoigt' ) {

    frac.lorenz <- d$l.height/(d$l.height + d$g.height)
    d$height <- d$l.height + d$g.height
    d$width <- d$area/(pi*d$height)
    d <- d[, columns]

  } else if ( current.type == 'voigt' ) {

    frac.lorenz <- d$l.width/(d$l.width + d$g.width)
    d$width <- d$l.width + d$g.width
    d$height <- d$area/(pi*d$width)
    d <- d[, columns]

  }

  # Conversion from lorenz to desired peak
  if ( peak.type == 'gauss' ) {

    d$width <- d$area/(sqrt(2*pi)*d$height)
    d <- d[, .all_columns(object, peak.type = 'gauss')]

  } else if ( peak.type == 'pvoigt' ) {

    d$l.height <- frac.lorenz*d$height
    d$g.height <- (1 - frac.lorenz)*d$height
    new.width <- d$area/(pi*d$l.height + sqrt(pi)*d$g.height)
    d$width <- ifelse(is.finite(new.width), new.width, d$width)
    d <- d[, .all_columns(object, peak.type = 'pvoigt')]

  } else if ( peak.type == 'voigt' ) {

    d$l.width = frac.lorenz*d$width
    d$g.width = (1 - frac.lorenz)*d$width
    
    new.height = d$area*Re(Faddeeva_w(complex(im=d$l.width)/(sqrt(2)*d$g.width)))/
                       (sqrt(2*pi)*d$g.width)
    d$height <- ifelse(is.finite(new.height), new.height, d$height)
    d <- d[, .all_columns(object, peak.type = 'voigt')]

  }

  # If our target was lorenz, still need to drop area 
  if ( peak.type == 'lorenz' ) {
    d <- d[ , colnames(d) != 'area']
  }

  # Insert new peak parameters
  object@peak_type <- peak.type
  object@peaks <- d

  # Update parameter slot
  object <- .merge_parameters(object)

  # Bounds must have the same peak_type
  .propagate_to_bounds(object, set_peak_type, peak.type, frac.lorenz)
})

#------------------------------------------------------------------------
#' @rdname normalized
#' @export
setMethod("normalized", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname normalized-set
#' @export
setReplaceMethod("normalized", "NMRScaffold1D",
                 function(object, value) callNextMethod())

#------------------------------------------------------------------------
#' @rdname set_normalized
#' @export
setMethod("set_normalized", "NMRScaffold1D", 
  function(object, nmrdata = NULL, normalized = TRUE, 
           include.bounds = FALSE) {

  # Updating bounds, first, if necessary
  if ( include.bounds ) {
    object <- .propagate_to_bounds(object, set_normalized, 
                                   nmrdata = nmrdata, normalized = normalized)
  } 

  # If target normalization state matches current state, return
  if ( object@normalized == normalized ) return(object)

  # Checking nmrdata
  nmrdata <- .check_data(object, nmrdata)

  # Extracting the processed portion of nmrdata
  processed <- nmrdata@processed

  peaks <- object@peaks
  baseline <- object@baseline
  baseline.diff <- object@baseline_difference
  knots <- object@knots
  constraints <- object@constraints
  peak.type <- object@peak_type

  # Calculate scaling factors
  x.offset <- min(processed$direct.shift)
  x.scale <- max(processed$direct.shift) - x.offset

  y.scale <- max(abs(Re(processed$intensity)))

  p.columns <- .position_columns(object, TRUE)
  h.columns <- .height_columns(object, TRUE)
  w.columns <- .width_columns(object, TRUE)

  if ( normalized ) {
    # Scaling based on peak type
    peaks[, p.columns] <- (peaks[, p.columns] - x.offset)/x.scale
    peaks[, h.columns] <- peaks[, h.columns]/y.scale
    peaks[, w.columns] <- peaks[, w.columns]/x.scale

    # Scaling couplin relationships
    constraints$difference <- constraints$difference/x.scale

    # Baseline
    baseline <- baseline/y.scale
    baseline.diff <- baseline.diff/y.scale
    knots <- (knots - x.offset)/x.scale
  } else {
    # Scaling based on peak type
    peaks[, p.columns] <- (peaks[, p.columns]*x.scale) + x.offset
    peaks[, h.columns] <- peaks[, h.columns]*y.scale
    peaks[, w.columns] <- peaks[, w.columns]*x.scale

    # Scaling coupling relationships
    constraints$difference <- constraints$difference*x.scale

    # Baseline
    baseline <- baseline*y.scale
    baseline.diff <- baseline.diff*y.scale
    knots <- knots*x.scale + x.offset
  }

  # Setting new values
  object@peaks <- peaks
  object@baseline <- baseline
  object@baseline_difference <- baseline.diff
  object@knots <- knots
  object@constraints <- constraints
  object@normalized <- normalized 

  # Update parameter slot
  object <- .merge_parameters(object)
})

#------------------------------------------------------------------------
#' @rdname peak_units
#' @export
setMethod("peak_units", "NMRScaffold1D", 
          function(object) callNextMethod())

#' @rdname peak_units-set
#' @export
setReplaceMethod("peak_units", "NMRScaffold1D",
                 function(object, value) callNextMethod())

#------------------------------------------------------------------------
#' @rdname set_peak_units
#' @export
setMethod("set_peak_units", "NMRScaffold1D",
  function(object, nmrdata = NULL, peak.units = 'hz', 
           include.bounds = TRUE) callNextMethod())



#========================================================================>
# Calculations based on current parameter values.
#========================================================================>



#------------------------------------------------------------------------
#' @rdname f_lineshape
#' @export
setMethod("f_lineshape", "NMRScaffold1D",
  function(object) {

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
#' @rdname f_baseline
#' @export
setMethod("f_baseline", "NMRScaffold1D",
  function(object) {

  # Converting object peak_units to ppm
  object <- set_peak_units(object, peak.units = 'ppm')

  # Extracting parameters
  baseline <- object@baseline
  baseline.diff <- object@baseline_difference
  knots <- object@knots

  degree1 <- length(baseline) - length(knots) - 1
  degree2 <- length(baseline.diff) - length(knots) - 1
  if ( length(knots) == 0 ) knots <- NULL

  f <- function(x, include.convolution = TRUE) {

    if ( degree1 >= 1 ) {
      basis <- bs(x, degree = degree1, knots = knots, intercept = TRUE)
      y <- c(basis %*% matrix(baseline, ncol = 1))
    } else {
      y <- rep(0, length(x))
    }

    # Tack on baseline difference
    if ( degree2 >= 0 ) {
      if ( degree2 != degree1 ) {
        basis <- bs(x, degree = degree2, knots = knots, intercept = TRUE) 
      }
      y.diff <- c(basis %*% matrix(baseline.diff, ncol = 1))
    } else {
      y.diff <- rep(0, length(x))
    }

    complex(re = y, im = y + y.diff )
  }

  f
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
#' @rdname calc_baseline 
#' @export
setMethod("calc_baseline", "NMRScaffold1D",
  function(object, direct.shift = NULL) {

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
  f <- f_baseline(object)

  # Apply it
  f(direct.shift)
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



#========================================================================>
#  Constraint generation
#========================================================================>



#------------------------------------------------------------------------
# Checks simple bounds to make sure they are valid
.check_bounds <- function(constraints) {

  if ( length(constraints) != 2 ) {
    msg <- paste("All constraints must be vectors of two elements consisting",
                 "of a lower and upper bound.")
    stop(msg)
  }

  if ( constraints[1] > constraints[2] ) {
    msg <- paste("Lower bound must be greater than upper bound.",
                 "Proceeding with current constraints will result in a",
                 "fit error.")
    warning(msg)
  }

}

#------------------------------------------------------------------------
# Updates bounds with an option to prevent widening
.update_bounds <- function(old.lower, new.lower,
                           old.upper, new.upper, widen = FALSE) {

  if ( length(new.lower) == 1) new.lower <- rep(new.lower, length(old.lower))
  if ( length(new.upper) == 1) new.upper <- rep(new.upper, length(old.upper))

  logic.lower <- is.finite(new.lower)
  logic.upper <- is.finite(new.upper)

  if (! widen) {
    logic.lower <- logic.lower & (new.lower > old.lower)
    logic.upper <- logic.upper & (new.upper < old.upper)
  }

  old.lower[logic.lower] <- new.lower[logic.lower]
  old.upper[logic.upper] <- new.upper[logic.upper]

  list(lower = old.lower, upper = old.upper)
}

#------------------------------------------------------------------------
#' @rdname initialize_bounds
setMethod(".initialize_bounds", "NMRScaffold1D", 
  function(object, overwrite = FALSE) {

    # Handling bounds if they exist
    bounds <- list(lower = object@bounds$lower, 
                   upper = object@bounds$upper)

    # Selecting default values
    values <- list(lower = -Inf, upper = +Inf)

    for ( name in names(bounds) ) {
      if ( is.null(bounds[[name]]) || overwrite ) {
        bound <- object

        peaks <- bound@peaks
        d.columns <- .data_columns(object, TRUE)
        peaks[ , d.columns] <- values[[name]]
        bound@peaks <- peaks

        # Since baseline and phase elements must be numeric, ensure that NA
        # values take numeric form
        bound@baseline <- rep(values[[name]], length(bound@baseline))
        bound@baseline_difference <- rep(values[[name]], 
                                         length(object@baseline_difference))
        bound@phase <- rep(values[[name]], length(object@phase))

        bound@constraints <- data.frame(
          id1 = character(), id2 = character(), 
          peak1 = character(), peak2 = character(),
          difference = numeric(), ratio = numeric())

        bounds[[name]] <- .merge_parameters(bound, TRUE)
      }
    }

    object@bounds$lower <- bounds$lower
    object@bounds$upper <- bounds$upper

    object
  })

#------------------------------------------------------------------------
#' @rdname set_absolute_bounds
#' @export
setMethod("set_absolute_bounds", "NMRScaffold1D",
  function(object, position = NULL, height = NULL, width = NULL, 
           baseline = NULL, phase = NULL, 
           normalized = FALSE, peak.units = 'hz', widen = FALSE) {

  # Initializing bounds
  object <- .initialize_bounds(object)
  lower <- object@bounds$lower
  upper <- object@bounds$upper

  # Ensuring proper normalization and peak units
  lower <- set_normalized(lower, normalized = normalized)
  upper <- set_normalized(upper, normalized = normalized)

  lower <- set_peak_units(lower, peak.units = peak.units)
  upper <- set_peak_units(upper, peak.units = peak.units)

  # Getting current values
  l.peaks <- lower@peaks
  u.peaks <- upper@peaks

  l.baseline <- lower@baseline
  u.baseline <- upper@baseline

  l.baseline.diff <- lower@baseline_difference
  u.baseline.diff <- upper@baseline_difference

  l.phase <- lower@phase
  u.phase <- upper@phase

  p.columns <- .position_columns(object, TRUE)
  h.columns <- .height_columns(object, TRUE)
  w.columns <- .width_columns(object, TRUE)

  # Position
  if (! is.null(position) ) {
    .check_bounds(position)
    new.bounds <- .update_bounds(l.peaks[, p.columns], position[1],
                                u.peaks[, p.columns], position[2], widen)
    l.peaks[, p.columns] <- new.bounds$lower 
    u.peaks[, p.columns] <- new.bounds$upper
  }

  # Height
  if (! is.null(height) ) {
    .check_bounds(height)
    # If there are multiple width columns, flatten them
    n <- length(h.columns)
    new.bounds <- .update_bounds(unlist(l.peaks[, h.columns]), height[1],
                                 unlist(u.peaks[, h.columns]), height[2], widen)
    l.peaks[, h.columns] <- matrix(new.bounds$lower, ncol = n)
    u.peaks[, h.columns] <- matrix(new.bounds$upper, ncol = n)
  }

  # Width
  if (! is.null(width) ) {
    .check_bounds(width)
    # If there are multiple width columns, flatten them
    n <- length(w.columns)
    new.bounds <- .update_bounds(unlist(l.peaks[, w.columns]), width[1],
                                 unlist(u.peaks[, w.columns]), width[2], widen)
    l.peaks[, w.columns] <- matrix(new.bounds$lower, ncol = n)
    u.peaks[, w.columns] <- matrix(new.bounds$upper, ncol = n)
  }

  # Baseline -- both real and imaginary components treated in the same way
  if (! is.null(baseline) ) {
    .check_bounds(baseline)
    new.bounds <- .update_bounds(l.baseline, baseline[1],
                                u.baseline, baseline[2], widen)
    l.baseline <- new.bounds$lower
    u.baseline <- new.bounds$upper

    new.bounds <- .update_bounds(l.baseline.diff, baseline[1],
                                u.baseline.diff, baseline[2], widen)
    l.baseline.diff <- new.bounds$lower
    u.baseline.diff <- new.bounds$upper
  }

  # Phase
  if (! is.null(phase) ) {
    .check_bounds(phase)
    new.bounds <- .update_bounds(l.phase, phase[1],
                                u.phase, phase[2], widen)
    l.phase <- new.bounds$lower
    u.phase <- new.bounds$upper
  }

  # Combine parameters
  l.parameters <- .gen_parameters(lower, l.peaks, l.baseline, 
                                          l.baseline.diff, l.phase) 
  u.parameters <- .gen_parameters(upper, u.peaks, u.baseline, 
                                          u.baseline.diff, u.phase) 

  # Generating and setting scaffold objects
  lower <- new("NMRScaffold1D", lower, peaks = l.peaks, baseline = l.baseline,
               baseline_difference = l.baseline.diff, phase = l.phase, 
               parameters = l.parameters)

  upper <- new("NMRScaffold1D", upper, peaks = u.peaks, baseline = u.baseline,
               baseline_difference = u.baseline.diff, phase = u.phase, 
               parameters = u.parameters)

  # Setting the boundaries
  object@bounds$lower <- lower
  object@bounds$upper <- upper

  # Propagating normalization and peak units
  object <- set_normalized(object, normalized = object@normalized, 
                           include.bounds = TRUE)
  object <- set_peak_units(object, peak.units = object@peak_units,
                           include.bounds = TRUE)

  object
})

#------------------------------------------------------------------------
#' @rdname set_relative_bounds
#' @export
setMethod("set_relative_bounds", "NMRScaffold1D",
  function(object, position = NULL, height = NULL, width = NULL, 
           baseline = NULL, phase = NULL, 
           normalized = TRUE, widen = FALSE) {

  # Initializing reference values for fractions 
  reference <- set_normalized(object, normalized = normalized)

  r.peaks <- reference@peaks
  r.baseline <- reference@baseline
  r.baseline.diff <- reference@baseline_difference
  r.phase <- reference@phase

  p.columns <- .position_columns(reference, TRUE)
  h.columns <- .height_columns(reference, TRUE)
  w.columns <- .width_columns(reference, TRUE)

  # Initializing bounds
  reference <- .initialize_bounds(reference)
  lower <- reference@bounds$lower
  upper <- reference@bounds$upper

  # Ensuring proper normalization
  lower <- set_normalized(lower, normalized = normalized)
  upper <- set_normalized(upper, normalized = normalized)

  # Getting current values
  l.peaks <- lower@peaks
  u.peaks <- upper@peaks

  l.baseline <- lower@baseline
  u.baseline <- upper@baseline

  l.baseline.diff <- lower@baseline_difference
  u.baseline.diff <- upper@baseline_difference

  l.phase <- lower@phase
  u.phase <- upper@phase

  # Position
  if (! is.null(position) ) {
    .check_bounds(position)
    l.position <- r.peaks[, p.columns]*position[1]
    u.position <- r.peaks[, p.columns]*position[2]
    new.bounds <- .update_bounds(l.peaks[, p.columns], l.position,
                                 u.peaks[, p.columns], u.position, widen)
    l.peaks[, p.columns] <- new.bounds$lower 
    u.peaks[, p.columns] <- new.bounds$upper
  }

  # Height
  if (! is.null(height) ) {
    .check_bounds(height)
    # If there are multiple height columns, flatten them
    n <- length(w.columns)
    l.height <- unlist(r.peaks[, h.columns])*height[1]
    u.height <- unlist(r.peaks[, h.columns])*height[2]
    new.bounds <- .update_bounds(unlist(l.peaks[, h.columns]), l.height,
                                 unlist(u.peaks[, h.columns]), u.height, widen)
    l.peaks[, h.columns] <- matrix(new.bounds$lower, ncol = n)
    u.peaks[, h.columns] <- matrix(new.bounds$upper, ncol = n)
  }

  # Width
  if (! is.null(width) ) {
    .check_bounds(width)
    # If there are multiple width columns, flatten them
    n <- length(w.columns)
    l.width <- unlist(r.peaks[, w.columns])*width[1]
    u.width <- unlist(r.peaks[, w.columns])*width[2]
    new.bounds <- .update_bounds(unlist(l.peaks[, w.columns]), l.width,
                                 unlist(u.peaks[, w.columns]), u.width, widen)
    l.peaks[, w.columns] <- matrix(new.bounds$lower, ncol = n)
    u.peaks[, w.columns] <- matrix(new.bounds$upper, ncol = n)
  }

  # Baseline -- both real and imaginary components treated in the same way
  if (! is.null(baseline) ) {
    .check_bounds(baseline)
    new.bounds <- .update_bounds(l.baseline, r.baseline*baseline[1],
                                u.baseline, r.baseline*baseline[2], 
                                widen)
    l.baseline <- new.bounds$lower
    u.baseline <- new.bounds$upper

    new.bounds <- .update_bounds(l.baseline.diff, r.baseline.diff*baseline[1],
                                u.baseline.diff, r.baseline.diff*baseline[2], 
                                widen)
    l.baseline.diff <- new.bounds$lower
    u.baseline.diff <- new.bounds$upper
  }

  # Phase
  if (! is.null(phase) ) {
    .check_bounds(phase)
    new.bounds <- .update_bounds(l.phase, r.phase*phase[1],
                                u.phase, r.phase*phase[2], widen)
    l.phase <- new.bounds$lower
    u.phase <- new.bounds$upper
  }

  # Combine parameters
  l.parameters <- .gen_parameters(lower, l.peaks, l.baseline, 
                                         l.baseline.diff, l.phase) 
  u.parameters <- .gen_parameters(upper, u.peaks, u.baseline, 
                                         u.baseline.diff, u.phase) 

  # Generating and setting scaffold objects
  lower <- new("NMRScaffold1D", lower, peaks = l.peaks, baseline = l.baseline,
               baseline_difference = l.baseline.diff, phase = l.phase, 
               parameters = l.parameters)

  upper <- new("NMRScaffold1D", upper, peaks = u.peaks, baseline = u.baseline,
               baseline_difference = u.baseline.diff, phase = u.phase, 
               parameters = u.parameters)

  # Setting the boundaries
  object@bounds$lower <- lower
  object@bounds$upper <- upper

  # Propagating normalization and peak units
  object <- set_normalized(object, normalized = object@normalized, 
                           include.bounds = TRUE)
  object <- set_peak_units(object, peak.units = object@peak_units,
                           include.bounds = TRUE)

  object
  })



#------------------------------------------------------------------------
#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRScaffold1D",
  function(object, position = TRUE,  height = TRUE, width = TRUE, 
           baseline = TRUE, phase = TRUE, widen = FALSE) { 

  # First, do a single pass over absolute normalized bounds
  if ( position ) {
    abs.position <- c(0, 1)
  } else {
    abs.position <- NULL
  }

  if ( height ) {
    abs.height <- c(0, 1.5)
  } else {
    abs.height <- NULL
  }

  if ( baseline ) {
    abs.baseline <- c(-0.25, 0.25)
  } else {
    abs.baseline <- NULL
  }

  if ( phase ) {
    abs.phase <- c(-pi/2, pi/2)
  } else {
    abs.phase <- NULL
  }

  object <- set_absolute_bounds(object, position = abs.position, 
                                height = abs.height, baseline = abs.baseline,
                                phase = abs.phase, normalized = TRUE,
                                widen = widen)

  # Width constraints are not normalized
  if ( width ) {
    object <- set_absolute_bounds(object, width = c(0.003, 3), 
                                  normalized = FALSE, peak.units = 'hz',
                                  widen = widen)
  }

  # Finally, the only relative constraints are on position
  if ( position ) {
    object <- set_relative_bounds(object, position = c(0.8, 1.2), 
                                  normalized = TRUE, widen = widen)
  }

  object
  })
