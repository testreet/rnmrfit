# Definition of a class structure for 1D peak scaffold data.



#========================================================================>
# Documentation entries
#========================================================================>



# Ensuring that NMRScaffold is loaded
#' @include NMRScaffold.R
NULL

#' Inherited ellipses description
#' @param ... Additional arguments passed to inheriting methods.
#' @name methodEllipse
NULL



#========================================================================>
#  NMRScaffold1D -- peak description 
#========================================================================>



# Defining a union between NMRData1D and NULL to use in slot
NMRData1DorNULL <- setClassUnion('NMRData1DorNULL',
                                 c('NMRData1D', 'NULL'))

#------------------------------------------------------------------------
#' A class representing a set of NMR peaks and their relationships. 
#'
#' TO DO.
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
#' @slot convolution TO DO.
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
#' @slot lower_bounds An NMRScaffold1D object that sets lower feasible bounds
#'                    on all parameters during peak fitting. Can be set to NULL 
#'                    to leave unbounded
#' @slot upper_bounds An NMRScaffold1D object that sets upper feasible bounds
#'                    on all parameters during peak fitting. Can be set to NULL 
#'                    to leave unbounded
#' @slot nmrdata An optional NMRData1D object that serves as a reference for
#'               for normalization and peak_unit conversion. However, this 
#'               can also be provided to the individual functions as needed.
#'
#' @name NMRScaffold1D-class
#' @export
NMRScaffold1D <- setClass("NMRScaffold1D",
                          contains = "NMRScaffold",
                          slots = c(baseline = 'numeric',
                                    baseline_difference = 'numeric',
                                    knots = 'numeric',
                                    phase = 'numeric',
                                    convolution = 'numeric',
                                    nmrdata = 'NMRData1DorNULL'),
                          prototype = prototype(normalized = FALSE,
                                                peak_type = 'lorenz',
                                                peak_units = 'hz',
                                                bounds = list(lower = NULL,
                                                              upper = NULL),
                                                nmrdata = NULL))



#========================================================================>
# Validation methods
#========================================================================>

# Defining specific cases for baseline and phase

#------------------------------------------------------------------------
#' NMRScaffold1D validity test
#' 
#' Appends baseline, phase, and parameter checks on top of the checks
#' handled by the NMRScaffold super class
validNMRScaffold1D <- function(object) {

  peaks <- object@peaks
  baseline <- object@baseline
  baseline.diff <- object@baseline_difference
  knots <- object@knots
  phase <- object@phase 
  parameters <- object@parameters
  constraints <- object@constraints
  bounds <- object@bounds

  valid <- TRUE
  msg <- c()

  # Checking peak column names
  valid.columns <- .all_columns(object)
  if (! identical(colnames(peaks), valid.columns) ) {
    valid <- FALSE
    new.msg <- sprintf('"peaks" must have the following columns: %s',
                          paste(valid.columns, collapse = ', '))
    msg <- c(msg, new.msg)
  }
    
  # The parameter vector must be longer than the knot vector
  if ( (length(baseline) > 0) && (length(baseline) <= length(knots)) ) {
    valid <- FALSE
    new.msg <- paste('"knots" vector must be shorter than "baseline" vector:',
               'n.baseline = n.knots + baseline.degree + 1')
    msg <- c(msg, new.msg)
  }

  # Phase length
  if ( (length(phase) > 0) && (length(phase) != 1) ) {
    valid <- FALSE
    new.msg <- '"phase", if it exists, must be of length 1.'
    msg <- c(msg, new.msg)
  }

  # Making sure that individual parameters are equivalent to combined vector
  current.parameters <- .merge_parameters(object, FALSE)
  if (! identical(parameters, current.parameters) ) {
    valid <- FALSE
    new.msg <- paste('The combined "parameters" slot does not match values', 
                        'of the individual slots.')
    msg <- c(msg, new.msg)
  }

  # Checking constraint columns
  valid.columns <- c('id1', 'id2', 'peak1', 'peak2', 
                     'difference', 'ratio')

  if (! identical(colnames(constraints), valid.columns)) {
    valid <- FALSE
    new.msg <- sprintf('"constraints" must have the following columns: %s',
                          paste(valid.columns, collapse = ', '))
    msg <- c(msg, new.msg)
  } 

  # Checking bounds
  valid.bounds <- c('lower', 'upper')
  if ( (length(bounds) > 0) && (any(!names(bounds) %in% valid.bounds)) ) {
    valid <- FALSE
    msg <- c(msg, '"bounds" can only have "lower" and "upper" elements.')
  }

  lower <- .is_conformant(object, bounds$lower)
  upper <- .is_conformant(object, bounds$upper)

  string <- c()
  if (! lower ) string <- c(string, 'lower')
  if (! upper ) string <- c(string, 'upper')
  string <- paste(string, collapse = ' and ')

  if (string != '' ) {
    valid <- FALSE
    new.msg <- sprintf('The shape of %s bounds does not conform to parameters', 
                       new.msg)
    msg <- c(msg, string)
  }

  if (valid) TRUE
  else msg
}

# Add the extended validity testing
setValidity("NMRScaffold1D", validNMRScaffold1D)



#========================================================================>
# Constructor for initialization
#========================================================================>



#------------------------------------------------------------------------
#' Generate an NMRScaffold1D object based on simplified peak list
#'
#' Generates an NMRScaffold1D object by converting multiplet definitions into
#' a set of singlets related by constraints on their position and area.
#  An NMRData1D object is required to come up with rough estimates for peak 
#' height and convert coupling constants from Hz to ppm. If no data is 
#' provided, peak heights are assumed to be 0. This is a bad initial 
#' guess and is unlikely to result in a good fit.
#'
#' Multiplets are specified in the from of "3 d 1.2", which indicates
#' a doublet at 3 ppm, with a coupling constant of 1.2 Hz. Supported  
#' multiplet codes include s, d, t, q, quint, sext, sept, oct, dd, dt, 
#' td, and tt. The latter 4 expect two coupling values. The exact 
#' syntax of the specification is quite flexible, as the function  
#' removes all brackets and most punctuation and converts common spacers 
#' like - and _ into whitespace. Satellite coupling constants (in Hz) can be
#' added using "j" in the specification, so "3 s j 10" would indicate
#' a singlet at 3 ppm and 1 pair of satellites with a coupling constant of
#' 10 Hz. It's also possible to offset the center of the satellite
#' doublet from the center of the main isotopic peak using a colon. In the
#' specification of "3 s j 10:0.01", the center of the satellite doublet is
#' offset 0.01 ppm upfield. 
#'
#' @param peak.list A vector or list of character strings specifying
#'                  multiplets of the form "3 d 1.2". See details section.
#' @param peak.type One of lorenz, gauss, pvoigt, or voigt.
#' @param nmrdata An NMRData1D object.
#' @param peak.width Initial estimate of peak width (in Hz)
#' @param baseline.degree Degree of B-spline polynomials to use for baseline.
#'                        As the degree can be 0, use NULL to remove
#'                        baseline entirely (default value).
#' @param n.knots The number of equally spaced interior knots used for 
#'                B-spline baseline. The knots can be changed after the
#'                scaffold in initialized.
#' @param include.difference Include a baseline difference term to account
#'                           for slight difference between real and imaginary
#'                           components. By default, this difference is
#'                           initialized to be the same as the baseline
#'                           term. This can be changed later.
#' @param include.phase Include a phasing term in the fit.
#' @param store FALSE to not store the nmrdata object. Although storing
#'              an nmrdata object is likely to be more convenient,
#'              it may result in unnecessary data copies.
#'
#' @return An NMRScaffold1D object.
#'
#' @export
nmrscaffold_1d <- function(peak.list, nmrdata, peak.type = 'lorenz', 
                           peak.width = 1, baseline.degree = 3, n.knots = 3,
                           include.difference = TRUE, include.phase = TRUE, 
                           store = TRUE) {

  # Checking nmrdata
  if ( class(nmrdata) == 'NMRData1D' ) {
    validObject(nmrdata)
  } else {
    msg <- '"nmrdata" must be a valid NMRData1D object.'
    stop(msg)
  }

  # Forcing include.difference to FALSE if there is no baseline
  if ( is.null(baseline.degree) && ( include.difference ) ) {
    msg <- paste('Setting "include.difference" to TRUE without a baseline',
                 'can not be processed, ignoring')
    warning(msg, call. = FALSE)
    include.difference <- FALSE
  }

  # Parsing peak.list
  parsed <- parse_coupling(peak.list)
  ids <- parsed$ids
  chemical.shift <- parsed$chemical.shift
  peaks <- parsed$peaks
  coupling <- parsed$coupling

  # Converting Hertz coupling into ppm
  sfo1 <- nmrdata@acqus[['sfo1']]
  coupling <- lapply(coupling, function(x) x/sfo1)

  # Initializing output
  peaks.list = list()
  constraints.list = list()

  # Parsing singlets (they have no constraints)
  is.s <- unlist(lapply(coupling, function(x) all(is.na(x))) )
  for (i in which(is.s)) {
    peaks.frame <- data.frame(id = ids[i], peak = 1, 
                              position = chemical.shift[i])
    peaks.list[[as.character(ids[i])]] <- peaks.frame
  }

  # Parsing multiplets
  for (i in which(!is.s)) {
    singlets <- split_coupling(peaks[[i]], coupling[[i]])
    singlets[, 1] <- singlets[, 1] + chemical.shift[i]
    n <- nrow(singlets)
    
    # First the peaks 
    peaks.frame <- data.frame(id = ids[i], peak = 1:n, 
                              position = singlets[, 1])
    peaks.list[[as.character(ids[i])]] <- peaks.frame

    # Then the constraints
    differences <- singlets[2:n, 1] - singlets[1:(n-1), 1]
    ratios = singlets[2:n, 2]/singlets[1:(n-1), 2]
    
    constraints.frame <- data.frame(id1 = ids[i], id2 = ids[i], 
                                    peak1 = 1:(n-1), peak2 = 2:n,
                                    difference = differences, ratio = ratios)

    constraints.list[[as.character(ids[i])]] <- constraints.frame
  }

  # Stitching together lists of data.frames
  peaks <- bind_rows(peaks.list)
  constraints <- bind_rows(constraints.list)

  n.peaks <- nrow(peaks)
  n.constraints <- nrow(constraints)

  # Initializing constraint and peak logic
  satellite.peaks <- list(rep(FALSE, n.peaks))
  satellite.constraints <- list(rep(FALSE, n.constraints))
  coupling.constraints <- list(rep(TRUE, n.constraints))

  # Parsing satellite data
  sat.ids <- parsed$sat.ids
  sat.offsets <- parsed$sat.offsets
  sat.numbers <- parsed$sat.numbers
  sat.coupling <- parsed$sat.coupling

  # Looping through sat.ids if there are any
  for (id in sat.ids) {

    id.factor <- factor(id, levels = levels(sat.ids))

    np <- nrow(peaks.list[[id]])
    last.peak <- np

    # Adding a doublet per peak per satellite
    for (j in 1:sat.numbers[[id]]) {
      coupling <- sat.coupling[[id]][j]
      offsets <- sat.offsets[[id]][j]

      left.position <- peaks.list[[id]]$position - 
                       coupling/2/sfo1 -
                       offsets
      right.position <- peaks.list[[id]]$position + 
                        coupling/2/sfo1 -
                        offsets

      # Interleaving and adding peak numbers
      peak.numbers <- seq(last.peak + 1, 2*np + last.peak)
      peak.positions <- c(rbind(left.position, right.position))

      peaks.frame <- data.frame(id = id.factor,
                                peak = peak.numbers,
                                position = peak.positions)

      peaks.list[[paste(id, 'j', j, sep = '-')]] <- peaks.frame
      last.peak <- last.peak + 2*np

      # Generating coupling constraints between satellites
      constraints.frame1 <- data.frame(id1 = id.factor, id2 = id.factor,
                                       peak1 = peak.numbers[seq(1, 2*np, 2)],
                                       peak2 = peak.numbers[seq(2, 2*np, 2)],
                                       difference = coupling/sfo1,
                                       ratio = 1)

      # Generating coupling constraints between satellites and main peaks
      constraints.frame2 <- data.frame(id1 = id.factor, id2 = id.factor,
                                       peak1 = peak.numbers[seq(1, 2*np, 2)],
                                       peak2 = 1:np,
                                       difference = coupling/2/sfo1 - 
                                                    offsets,
                                       ratio = 1)

      constraints.frame <- rbind(constraints.frame1, constraints.frame2)
      constraints.list[[paste(id, 'j', j, sep = '-')]] <- constraints.frame    

      # Tacking on constraint and peak logic
      satellite.peaks <- c(satellite.peaks, rep(TRUE, 2*np))
      satellite.constraints <- c(satellite.constraints, rep(TRUE, 2*np))
      coupling.logic <- rep(c(TRUE, FALSE), each = np)
      coupling.constraints <- c(coupling.constraints, coupling.logic)
    }
  }

  # Recombining following the addition of satellites
  peaks <- bind_rows(peaks.list)
  constraints <- bind_rows(constraints.list)

  satellite.peaks <- unlist(satellite.peaks) 
  satellite.constraints <- unlist(satellite.constraints)
  coupling.constraints <- unlist(coupling.constraints)

  # If there are no constraints, proper columns must still be enforced
  if ( nrow(constraints) == 0 ) {
    constraints <- data.frame(id1 = character(0), id2 = character(0),
                              peak1 = character(0), peak2 = character(0),
                              difference = numeric(0), ratio = numeric(0))
  }

  # Generate rough estimates for height
  d <- nmrdata@processed
  logic <- which_approx(d$direct.shift, peaks$position)
  height <- Re(d$intensity)[logic]
  height <- ifelse(height < 0, 0.1*max(Re(d$intensity)), height)

  # Adding on estimates for height and width
  peaks$height <- height
  peaks$width <- peak.width

  # Calculating scaling factors for knots
  x.offset <- min(d$direct.shift)
  x.scale <- max(d$direct.shift) - x.offset

  # Generating the scaffold object for lorenz peaks
  if ( n.knots > 0 ) {
    knots <- seq(0, 1, length.out = (n.knots + 2))[-c(1, n.knots+2)]
    knots <- knots*x.scale + x.offset
  } else {
    knots <- numeric(0)
  }

  if ( is.null(baseline.degree) ) { 
    baseline <- numeric(0)
  } else {
    baseline <- rep(0, n.knots + baseline.degree + 1)
  }

  if ( include.difference ) baseline.diff <- baseline
  else baseline.diff <- numeric(0)

  if ( include.phase ) phase <- 0
  else phase <- numeric(0)

  # Combining with a small hack to use the right version of .gen_parameters)
  object <- new('NMRScaffold1D')
  parameters <- .gen_parameters(object, peaks, baseline, baseline.diff, phase)

  # Dropping data from storage if desired
  if (! store ) nmrdata <- NULL

  scaffold <- new('NMRScaffold1D', peaks = peaks,
                  constraints = constraints, baseline = baseline,
                  baseline_difference = baseline.diff, knots = knots, 
                  phase = phase, nmrdata = nmrdata, parameters = parameters, 
                  .sat_peaks = satellite.peaks,
                  .sat_constraints = satellite.constraints,
                  .j_constraints = coupling.constraints)

  # Converting scaffold based on desired peak type
  set_peak_type(scaffold, peak.type)
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
    d$width <- d$area/(pi*d$l.height + sqrt(pi)*d$g.height)
    d <- d[, .all_columns(object, peak.type = 'pvoigt')]

  } else if ( peak.type == 'voigt' ) {

    d$l.width = frac.lorenz*d$width
    d$g.width = (1 - frac.lorenz)*d$width
    d$height = d$area*Re(Faddeeva_w(complex(im=d$l.width)/(sqrt(2)*d$g.width)))/
                      (sqrt(2*pi)*d$g.width)
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
  logic.upper <- is.finite(old.lower)

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
    new.bounds <- .update_bounds(l.peaks[, h.columns], height[1],
                                u.peaks[, h.columns], height[2], widen)
    l.peaks[, h.columns] <- new.bounds$lower 
    u.peaks[, h.columns] <- new.bounds$upper
  }

  # Width
  if (! is.null(width) ) {
    .check_bounds(width)
    new.bounds <- .update_bounds(l.peaks[, w.columns], width[1],
                                u.peaks[, w.columns], width[2], widen)
    l.peaks[, w.columns] <- new.bounds$lower 
    u.peaks[, w.columns] <- new.bounds$upper
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
    l.height <- r.peaks[, h.columns]*height[1]
    u.height <- r.peaks[, h.columns]*height[2]
    new.bounds <- .update_bounds(l.peaks[, h.columns], l.height,
                                u.peaks[, h.columns], u.height, widen)
    l.peaks[, h.columns] <- new.bounds$lower 
    u.peaks[, h.columns] <- new.bounds$upper
  }

  # Width
  if (! is.null(width) ) {
    .check_bounds(width)
    l.width <- r.peaks[, w.columns]*width[1]
    u.width <- r.peaks[, w.columns]*width[2]
    new.bounds <- .update_bounds(l.peaks[, w.columns], l.width,
                                u.peaks[, w.columns], u.width, widen)
    l.peaks[, w.columns] <- new.bounds$lower 
    u.peaks[, w.columns] <- new.bounds$upper
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
  print(lower_bounds(object))
  object <- set_normalized(object, normalized = object@normalized, 
                           include.bounds = TRUE)
  print(lower_bounds(object))
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
    abs.height <- c(0, 0.5)
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
    object <- set_absolute_bounds(object, width = c(0.3, 3), 
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
