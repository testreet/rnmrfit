# Definition of a class structure for peak fitting



#========================================================================>
#  Documentation entries
#========================================================================>



# Ensuring that NMRScaffold1D is loaded
#' @include NMRScaffold1D.R
NULL



#========================================================================>
#  NMRFit1D -- a combination of scaffold and data
#========================================================================>



#------------------------------------------------------------------------
#' A class representing a set of NMR peaks as well as a set of data.
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
#' @slot status Numeric status code of optimization.
#' @slot status_message Message issued by optimization.
#' @slot time A rough estimate of optimization time (in seconds).
#'
#' @name NMRFit1D-class
#' @export
NMRFit1D <- setClass('NMRFit1D',
                     slots = c(nmrdata = 'NMRData1D',
                               status = 'integer',
                               status_message = 'character',
                               time = 'numeric'),
                     contains = 'NMRScaffold1D')

#========================================================================>
#  Defining constraint functions for fitting
#========================================================================>

#------------------------------------------------------------------------
# Generate a constraint function on overall parameters based on lower bounds
.f_lower <- function(object) {

  # Get variables
  parameters <- object@parameters
  bounds <- object@bounds$lower@parameters

  mat = -diag(length(parameters))

  function(p) {
    constraints <- -p + bounds
    list(constraints = constraints, jacobian = mat)
  }
  
}

#------------------------------------------------------------------------
# Generate a constraint function on overall parameters based on upper bounds
.f_upper <- function(object) {

  # Get variables
  parameters <- object@parameters
  bounds <- object@bounds$upper@parameters

  mat = diag(length(parameters))

  function(p) {
    constraints <- p - bounds
    list(constraints = constraints, jacobian = mat)
  }
  
}

#------------------------------------------------------------------------
# Generate a constraint function on peak positions
.f_position <- function(object, leeway, logic) {

  # If there are no valid constraints, output NULL
  if (! any(logic) ) return(NULL)
  
  # Get variables
  peak.type <- object@peak_type
  peaks <- object@peaks
  parameters <- object@parameters
  constraints <- object@constraints[logic, ]

  n.peaks <- nrow(peaks)
  n.cols <- ncol(peaks) - 2
  n.parameters <- length(parameters)
  n.constraints <- nrow(constraints)

  # Generate a list of all positions
  p.index <- seq(.position_columns(object, TRUE) - 2, n.parameters, n.cols)
  j.rows <- 1:n.constraints

  # Subset the positions based on desired peaks
  logic.1 <- match(paste(constraints$id1, constraints$peak1),
                   paste(peaks$id, peaks$peak))
  peaks.1 <- p.index[logic.1]

  logic.2 <- match(paste(constraints$id2, constraints$peak2),
                   paste(peaks$id, peaks$peak))
  peaks.2 <- p.index[logic.2]

  # The constraints are based based on peak differences
  differences <- constraints$difference

  # Generate jacobian matrix (it will remain constant)
  mat <- matrix(0, nrow = n.constraints, ncol = n.parameters)
  mat[cbind(j.rows, peaks.2)] <- 1
  mat[cbind(j.rows, peaks.1)] <- -1

  # If no leeway, generate equality constraints, otherwise inequality
  if ( leeway == 0 ) {
    function(p) {
      constraints <- p[peaks.2] - p[peaks.1] - differences

      # For positions, the jacobian is constant
      list(constraints = constraints, jacobian = mat)
    }
  } else {
    function(p) {
      constraints <- c( p[peaks.2] - p[peaks.1] - (1 + leeway)*differences,
                       -p[peaks.2] + p[peaks.1] + (1 - leeway)*differences)

      # For positions, the jacobian is constant
      list(constraints = constraints, jacobian = rbind(mat, -mat))
    }
  }

}

#------------------------------------------------------------------------
# Generate a constraint function on peak position order
.f_order <- function(object, min.distance, ordered.peaks) {

  # If there are no ordered peaks, output NULL
  if ( length(ordered.peaks) == 0 ) return(NULL)
  
  # Get variables
  peak.type <- object@peak_type
  peaks <- object@peaks
  parameters <- object@parameters

  n.peaks <- nrow(peaks)
  n.cols <- ncol(peaks) - 2
  n.parameters <- length(parameters)

  # Generate a list of all positions
  p.index <- seq(.position_columns(object, TRUE) - 2, n.parameters, n.cols)

  # Subset the positions based on desired peaks
  logic <- match(ordered.peaks, peaks$id)
  n.constraints <- length(logic)

  # Ensuring correct order
  positions <- unlist(peaks[, .position_columns(object, TRUE)])[logic]
  peak.order <- order(positions)

  peaks.1 <- p.index[peak.order][1:(n.constraints-1)]
  peaks.2 <- p.index[peak.order][2:n.constraints]

  # The constraints are set from provided tolerance
  differences <- min.distance 

  # Generate jacobian 
  mat <- matrix(0, nrow = n.constraints-1, ncol = n.parameters)
  rows <- 1:(n.constraints - 1)
  mat[cbind(rows, peaks.2)] <- -1
  mat[cbind(rows, peaks.1)] <- 1

  function(p) {
    constraints <- -p[peaks.2] + p[peaks.1] + differences

    list(constraints = constraints, jacobian = mat)
  }

}

#------------------------------------------------------------------------
# Generate a constraint function on peak heights 
.f_height <- function(object, leeway, logic) {

  # If there are no valid constraints, output NULL
  if (! any(logic) ) return(NULL)
  
  # Get variables
  peak.type <- object@peak_type
  peaks <- object@peaks
  parameters <- object@parameters
  constraints <- object@constraints[logic, ]

  n.peaks <- nrow(peaks)
  n.cols <- ncol(peaks) - 2
  n.parameters <- length(parameters)
  n.constraints <- nrow(constraints)

  # Generate a list of all positions
  p.index <- seq(.height_columns(object, TRUE)[1] - 2, n.parameters, n.cols)
  j.rows <- 1:n.constraints

  # Subset the positions based on desired peaks
  logic.1 <- match(paste(constraints$id1, constraints$peak1),
                   paste(peaks$id, peaks$peak))
  peaks.1 <- p.index[logic.1]

  logic.2 <- match(paste(constraints$id2, constraints$peak2),
                   paste(peaks$id, peaks$peak))
  peaks.2 <- p.index[logic.2]

  # The constraints are based based on peak differences
  ratios <- constraints$ratio

  # Generate jacobian matrix (initialized here, but will not remain constant)
  mat <- matrix(0, nrow = n.constraints, ncol = n.parameters)

  # If the peak type is pvoigt, there are two heights to consider
  if ( peak.type == 'pvoigt' ) {
    j.rows <- 1:(n.constraints*2)
    peaks.1 <- c(peaks.1, peaks.1 + 1)
    peaks.2 <- c(peaks.2, peaks.2 + 1)
    mat <- rbind(mat, mat)
    ratios <- c(ratios, ratios)
  }

  # If no leeway, generate equality constraints, otherwise inequality
  if ( leeway == 0 ) {
    function(p) {
      constraints <- p[peaks.2]/p[peaks.1] - ratios

      mat[cbind(j.rows, peaks.2)] <- 1/p[peaks.1]
      mat[cbind(j.rows, peaks.1)] <- -p[peaks.2]/p[peaks.1]^2

      list(constraints = constraints, jacobian = mat)
    }
  } else {
    function(p) {
      constraints <- c( p[peaks.2]/p[peaks.1] - (1 + leeway)*ratios,
                       -p[peaks.2]/p[peaks.1] + (1 - leeway)*ratios)

      mat[cbind(j.rows, peaks.2)] <- 1/p[peaks.1]
      mat[cbind(j.rows, peaks.1)] <- -p[peaks.2]/p[peaks.1]^2

      list(constraints = constraints, jacobian = rbind(mat, -mat))
    }
  }

}

#------------------------------------------------------------------------
# Generate a constraint function on peak widths 
.f_width <- function(object, leeway, logic) {

  # If there are no valid constraints, output NULL
  if (! any(logic) ) return(NULL)
  
  # Get variables
  peak.type <- object@peak_type
  peaks <- object@peaks
  parameters <- object@parameters
  constraints <- object@constraints[logic, ]

  n.peaks <- nrow(peaks)
  n.cols <- ncol(peaks) - 2
  n.parameters <- length(parameters)
  n.constraints <- nrow(constraints)

  # Generate a list of all positions
  p.index <- seq(.width_columns(object, TRUE)[1] - 2, n.parameters, n.cols)
  j.rows <- 1:n.constraints

  # Subset the positions based on desired peaks
  logic.1 <- match(paste(constraints$id1, constraints$peak1),
                   paste(peaks$id, peaks$peak))
  peaks.1 <- p.index[logic.1]

  logic.2 <- match(paste(constraints$id2, constraints$peak2),
                   paste(peaks$id, peaks$peak))
  peaks.2 <- p.index[logic.2]

  # The constraints are based based on peak differences
  ratios <- constraints$ratio

  # Generate jacobian matrix (initialized here, but will not remain constant)
  mat <- matrix(0, nrow = n.constraints, ncol = n.parameters)

  # If the peak type is voigt, there are two widths to consider
  if ( peak.type == 'voigt' ) {
    j.rows <- 1:(n.constraints*2)
    peaks.1 <- c(peaks.1, peaks.1 + 1)
    peaks.2 <- c(peaks.2, peaks.2 + 1)
    mat <- rbind(mat, mat)
    ratios <- c(ratios, ratios)
  }

  # If no leeway, generate equality constraints, otherwise inequality
  if ( leeway == 0 ) {
    function(p) {
      constraints <- p[peaks.2]/p[peaks.1] - ratios

      mat[cbind(j.rows, peaks.2)] <- 1/p[peaks.1]
      mat[cbind(j.rows, peaks.1)] <- -p[peaks.2]/p[peaks.1]^2

      list(constraints = constraints, jacobian = mat)
    }
  } else {
    function(p) {
      constraints <- c( p[peaks.2]/p[peaks.1] - (1 + leeway)*ratios,
                       -p[peaks.2]/p[peaks.1] + (1 - leeway)*ratios)

      mat[cbind(j.rows, peaks.2)] <- 1/p[peaks.1]
      mat[cbind(j.rows, peaks.1)] <- -p[peaks.2]/p[peaks.1]^2

      list(constraints = constraints, jacobian = rbind(mat, -mat))
    }
  }

}

#------------------------------------------------------------------------
# The constraint function merely goes through a list of provided functions
# and combines the overall results
.f_constraint <- function(f.list) {
  constraints <- list()
  jacobian <- list()

  function(p) {
    for (i in 1:length(f.list)) {
      out <- f.list[[i]](p)
      constraints[[i]] <- out$constraints
      jacobian[[i]] <- out$jacobian
    }

    list(constraints = unlist(constraints), jacobian = do.call(rbind, jacobian))
  }
}

#------------------------------------------------------------------------
#' Fit an NMRScaffold1d object to a set of data
#'
#' To do.
#'
#' @param object An NMRScaffold1D object.
#' @param nmrdata An NMRData1D object. Optional if the NMRScaffold1D object
#'                has the nmrdata slot set. 
#' @param normalized TRUE to use a normalized fit. If the parameters or bounds
#'                   are not already normalized, they will be normalized to
#'                   the data. If FALSE, previously normalized data will be
#'                   rescaled before fit.
#' @param bounds TRUE to use bounds (setting conservative bounds if none exist),
#'               FALSE to ignore any attached bounds.
#' @param coupling.leeway Fraction of possible error in specified coupling
#'                        constants. Switching to a non-zero value results
#'                        in the use of inequality constraints.
#' @param area.leeway Fraction of possible error in specified area ratios of
#'                    coupled peaks. Switching to a non-zero value results
#'                    in the use of inequality constraints.
#' @param satellite.leeway Fraction of possible error in the symmetry of
#'                         satellite positions as a function of the distance
#'                         between a satellite and the main peak. Switching
#'                         to a non-zero value results in the use of 
#'                         inequality constraints.
#' @param ordered.peaks A vector of peak ids (names) whose order relative to
#'                      each other is not allowed to change. This option can
#'                      be useful for preventing an arbitrary collection of 
#'                      singlets used to fit a more complex peak from drifting
#'                      past or on top of each other, leading to ill-defined 
#'                      optimization.
#' @param include.convolution TRUE to include convolution vector of nmrdata
#'                            into fit (if present), FALSE to ignore it.
#' @param components String specifying the complex data componets to use for 
#'                   the fit (one of either "r/i", "r", or "i"). 
#' @param opts NLOPT option list used to overried default values.
#'
#' @return An NMRFit1D object.
#'
#' @export
nmrfit_1d <- function(object, nmrdata = NULL, normalized = TRUE, 
                      bounds = TRUE, coupling.leeway = 0, area.leeway = 0, 
                      satellite.leeway = 0, ordered.peaks = character(),
                      include.convolution = TRUE, components = "r/i", 
                      opts = NULL) {

  # Checking object 
  if ( class(object) %in% c('NMRScaffold1D', 'NMRFit1D') ) {
    validObject(object)
  } else {
    msg <- '"object" must be a valid NMRScaffold1D or NMRFit1D object.'
    stop(msg)
  }

  # Setting components vector
  if ( components == "r/i" ) components <- c(1, 1)
  else if ( components == "r" ) components <- c(1, 0)
  else if ( components == "i" ) components <- c(0, 1)
  else stop('"components" must be one of "r/i", "r", or "i"')

  # Checking nmrdata
  nmrdata <- .check_data(object, nmrdata)

  # Generating bounds if they don't exist
  lower <- object@bounds$lower
  upper <- object@bounds$upper

  if ( bounds ) {
    if ( is.null(lower) || is.null(upper) ) {
      bounded <- set_conservative_bounds(object, nmrdata)
    }

    if ( is.null(lower) ) lower <- bounded@bounds$lower
    if ( is.null(upper) ) upper <- bounded@bounds$upper

    object@bounds$lower <- lower
    object@bounds$upper <- upper
  }

  # Normalizing data
  processed <- nmrdata@processed

  x <- as.numeric(processed$direct.shift)
  y <- as.complex(processed$intensity)

  if ( normalized ) {
    x.offset <- min(x)
    x.scale <- max(x) - x.offset
    x <- (x - x.offset)/x.scale

    y.scale <- max(abs(Re(y)))
    y <- y/y.scale
  }

  # Normalizing parameters
  object <- set_normalized(object, nmrdata, normalized, TRUE)

  # Converting width to ppm
  object <- set_peak_units(object, nmrdata, 'ppm', TRUE)

  # Getting values
  peaks <- object@peaks
  baseline <- object@baseline
  baseline.diff <- object@baseline_difference
  knots <- object@knots
  phase <- object@phase

  peak.type <- object@peak_type

  # Flattening parameters
  p0 <- object@parameters

  # Defining general constraint functions. Since equality/inequality constraints
  # can vary, a list of functions is prepared and the actual constraint
  # function merely goes through each one and combines the output.
  j.constraints <- object@.j_constraints
  sat.constraints <- object@.sat_constraints

  eq <- list()
  ineq <- list()

  # First, adding coupled peak position differences
  if ( coupling.leeway == 0 ) {
    eq <- c(eq, .f_position(object, coupling.leeway, j.constraints))
  } else {
    ineq <- c(ineq, .f_position(object, coupling.leeway, j.constraints))
  }

  # Second, adding coupled peak height and width ratios
  if ( area.leeway == 0 ) {
    eq <- c(eq, .f_height(object, area.leeway, j.constraints))
    eq <- c(eq, .f_width(object, area.leeway, j.constraints))
  } else {
    ineq <- c(ineq, .f_height(object, area.leeway, j.constraints))
    ineq <- c(ineq, .f_width(object, area.leeway, j.constraints))
  }

  # Third, adding satellite position constraints  
  logic <- sat.constraints & (! j.constraints)
  if ( satellite.leeway == 0 ) {
    eq <- c(eq, .f_position(object, satellite.leeway, logic))
  } else {
    ineq <- c(ineq,.f_position(object, satellite.leeway, logic))
    #ineq <- c(ineq,.f_width(object, satellite.leeway, logic))
  }
  eq <- c(eq, .f_width(object, 0, logic))

  # Finally, forcing order on specified peaks
  if ( length(ordered.peaks) > 0 ) {
    ineq <- c(ineq, .f_order(object, x[2] - x[1], ordered.peaks))
  }

  # Adding simple lower and upper bounds if they exist
  if ( bounds ) {
    lower <- object@bounds$lower
    upper <- object@bounds$upper

    if (! is.null(object@bounds$lower) ) ineq <- c(ineq, .f_lower(object))
    if (! is.null(object@bounds$upper) ) ineq <- c(ineq, .f_upper(object))
  }

  # Generating overall constraint functions
  if ( length(eq) > 0 ) f_equality <- .f_constraint(eq)
  else f_equality <- NULL

  if ( length(ineq) > 0 ) f_inequality <- .f_constraint(ineq)
  else f_inequality <- NULL

  # Setting up default options
  d.opts <- list(algorithm = 'NLOPT_LD_SLSQP', 
                 xtol_rel = 1e-8,
                 maxeval = 200)

  # Overwriting defaults with opts argument
  if (! is.null(opts) ) {
    d.opts[names(opts)] = opts
  }

  # Defining optimization function
  n <- as.integer(nrow(processed))
  nc <- as.integer(nrow(peaks)*(ncol(peaks)-2))
  nb1 <- as.integer(length(baseline))
  nb2 <- as.integer(length(baseline.diff))
  nk <- as.integer(length(knots))
  np <- as.integer(length(phase))

  gradient <- as.numeric(rep(0, nc + nb1 + nb2 + np))
  evaluation <- as.numeric(0)
  components = as.integer(components)

  # Initializing fit environment for R fit
  if ( include.convolution && (length(nmrdata@convolution) > 0) ) {
    convolution <- nmrdata@convolution
  } else {
    convolution <- NULL
  }

  fit.env$init(x = x, y = y, nc = nc, nb1 = nb1, nb2 = nb2, knots = knots, 
               np = np, lineshape = peak.type, convolution = convolution, 
               components = NULL)

  #lineshapes = c('lorenz'=0, 'gauss'=1, 'pvoigt'=2, 'voigt'=3)
  #lineshape <- as.integer(lineshapes[peak.type])

#  # Initializing baseline basis
#  no1 <- as.integer(nb1 - nk - 1)
#  no2 <- as.integer(nb2 - nk - 1)
#
#  if ( nb1 > 0 ) {
#    .Fortran('init_basis1', x = x, n = n, knots = knots, nk = nk, no = no1)
#  }
#
#  if ( nb2 > 0 ) {
#    .Fortran('init_basis2', x = x, n = n, knots = knots, nk = nk, no = no2)
#  }
#
#  # Initialize convolution vector if required
#  if ( include.convolution && (length(nmrdata@convolution) > 0) ) {
#    
#    # Setting fortran vectors
#    n.convolution <- as.integer(2*n - 1)
#    .Fortran('init_convolution', x = nmrdata@convolution, n = n.convolution)
#
#    # Setting flag
#    convolution <- as.integer(1)
#  } else {
#    convolution <- as.integer(0)
#  }



  f_optim <- function(p) {
    p <- as.numeric(p)

    out <- .Fortran('fit_lineshape_1d', x = x, y = y, n = n, par = p, 
                    grad = gradient, nc = nc, nb1 = nb1, nb2 = nb2, np = np, 
                    eval = evaluation, lineshape = lineshape,
                    convolution = convolution, components = components)
    gradient <- out$grad
    evaluation <- out$eval
    list(objective = evaluation, gradient = gradient)
  }

  # Ignoring bounds if they aren't required
  if ( bounds ) {
    lb <- object@bounds$lower@parameters
    ub <- object@bounds$upper@parameters
  } else {
    lb <- NULL
    ub <- NULL
  }

  # Fitting
  start.time <- Sys.time()
  
  res <- nloptr(x0 = p0,
                eval_f = fit.env$eval,
                eval_g_eq = f_equality,
                eval_g_ineq = f_inequality,
                lb = lb,
                ub = ub,
                opts = d.opts)

  opt.time <- as.numeric(Sys.time() - start.time)

  if ( res$status < 0 ) {
    msg <- sprintf('Optimization issued error code (%i):', res$status)
    warning(paste(msg, res$message, sep = '\n'), call. = FALSE)
  }

  object@parameters <- res$solution

  # Restoring shape
  object <- .spread_parameters(object)

  # Restoring width to Hz
  object <- set_peak_units(object, nmrdata, 'hz')

  # Undoing normalization
  object <- set_normalized(object, nmrdata, FALSE, include.bounds = FALSE)

  # Outputting fit object
  new('NMRFit1D', object, normalized = FALSE, peak_units = 'hz', 
      bounds = list('lower'=lower, 'upper'=upper), nmrdata = nmrdata,
      status = res$status, status_message = res$message, time = opt.time)
}




#========================================================================>
# Basic setter and getter functions
#========================================================================>



#------------------------------------------------------------------------
#' @rdname peaks
#' @export
setMethod("peaks", "NMRFit1D", 
          function(object) object@peaks)

#' @rdname peaks-set
#' @export
setReplaceMethod("peaks", "NMRFit1D",
                 function(object, value) {
                   object@peaks <- value
                   object <- .merge_parameters(object)
                   object <- .drop_bounds(object)
                   validObject(object)
                   object 
                 })

#------------------------------------------------------------------------
#' @rdname baseline
#' @export
setMethod("baseline", "NMRFit1D", 
          function(object) object@baseline)

#' @rdname baseline-set
#' @export
setReplaceMethod("baseline", "NMRFit1D",
                 function(object, value) {
                   if ( is.null(value) ) value <- numeric(0)
                   object@baseline <- value
                   object <- .merge_parameters(object)
                   object <- .drop_bounds(object)
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @rdname baseline_difference
#' @export
setMethod("baseline_difference", "NMRFit1D", 
          function(object) object@baseline_difference)

#' @rdname baseline_difference-set
#' @export
setReplaceMethod("baseline_difference", "NMRFit1D",
                 function(object, value) {
                   if ( is.null(value) ) value <- numeric(0)
                   object@baseline_difference <- value
                   object <- .merge_parameters(object)
                   object <- .drop_bounds(object)
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @rdname baseline_knots
#' @export
setMethod("baseline_knots", "NMRFit1D", 
          function(object) object@knots)

#' @rdname baseline_knots-set
#' @export
setReplaceMethod("baseline_knots", "NMRFit1D",
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
setMethod("phase", "NMRFit1D", 
          function(object) object@phase)

#' @rdname phase-set
#' @export
setReplaceMethod("phase", "NMRFit1D",
                 function(object, value) {
                   if ( is.null(value) ) value <- numeric(0)
                   object@phase <- value
                   object <- .merge_parameters(object)
                   object <- .drop_bounds(object)
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @rdname constraints
#' @export
setMethod("constraints", "NMRFit1D", 
          function(object) object@constraints)

#' @rdname constraints-set
#' @export
setReplaceMethod("constraints", "NMRFit1D",
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
#' @rdname bounds
#' @export
setMethod("bounds", "NMRFit1D", 
          function(object) object@bounds)

#' @rdname bounds-set
#' @export
setReplaceMethod("bounds", "NMRFit1D",
                 function(object, value) {
                   object@bounds<- value
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @rdname lower_bounds
#' @export
setMethod("lower_bounds", "NMRFit1D", 
          function(object) object@bounds$lower)

#' @rdname lower_bounds-set
#' @export
setReplaceMethod("lower_bounds", "NMRFit1D",
                 function(object, value) {
                   object@bounds$lower <- value
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @rdname upper_bounds
#' @export
setMethod("upper_bounds", "NMRFit1D", 
          function(object) object@bounds$upper)

#' @rdname upper_bounds-set
#' @export
setReplaceMethod("upper_bounds", "NMRFit1D",
                 function(object, value) {
                   object@bounds$upper <- value
                   validObject(object)
                   object
                 })

#------------------------------------------------------------------------
#' @rdname nmrdata
#' @export
setMethod("nmrdata", "NMRFit1D", 
          function(object) object@nmrdata)

#' @rdname nmrdata-set
#' @export
setReplaceMethod("nmrdata", "NMRFit1D",
                 function(object, value) {
                   object@nmrdata <- value
                   validObject(object)
                   object 
                 })

#------------------------------------------------------------------------
#' @templateVar slot status
#' @template NMRScaffold_access
#' @name status 
#' @export
setGeneric("status", function(object, ...) 
           standardGeneric("status"))

#' @rdname status 
#' @export
setMethod("status", "NMRFit1D", 
          function(object) object@status)

#------------------------------------------------------------------------
#' @templateVar slot status_message
#' @template NMRScaffold_access
#' @name status_message
#' @export
setGeneric("status_message", function(object, ...) 
           standardGeneric("status_message"))

#' @rdname status_message 
#' @export
setMethod("status_message", "NMRFit1D", 
          function(object) object@status_message)




#========================================================================>
# Compound setter and getter functions (requiring internal calculations)
#========================================================================>




#------------------------------------------------------------------------
#' @rdname peak_type
#' @export
setMethod("peak_type", "NMRFit1D", 
          function(object) object@peak_type)

#' @rdname peak_type-set
#' @export
setReplaceMethod("peak_type", "NMRFit1D",
                 function(object, value) {
                   object <- set_peak_type(object, value)
                   validObject(object)
                   object
                 })

#' @rdname set_peak_type
#' @export
setMethod("set_peak_type", "NMRFit1D", 
  function(object, peak.type, frac.lorenz = 0.9) {
    callNextMethod()
  })

#------------------------------------------------------------------------
#' @rdname normalized 
#' @export
setMethod("normalized", "NMRFit1D", 
          function(object) object@normalized)

#' @rdname normalized-set
#' @export
setReplaceMethod("normalized", "NMRFit1D",
                 function(object, value) {
                   object <- set_normalized(object, normalized = value)
                   validObject(object)
                   object
                 })

#' @rdname set_normalized
#' @export
setMethod("set_normalized", "NMRFit1D", 
  function(object, nmrdata = NULL, normalized = TRUE, include.bounds = FALSE) {
    callNextMethod()
  })

#------------------------------------------------------------------------
#' @rdname peak_units
#' @export
setMethod("peak_units", "NMRFit1D", 
          function(object) object@peak_units)

#' @rdname peak_units
#' @export
setReplaceMethod("peak_units", "NMRFit1D",
                 function(object, value) {
                   object <- set_peak_units(object, peak.units = value)
                   validObject(object)
                   object
                 })

#' @rdname set_peak_units
#' @export
setMethod("set_peak_units", "NMRFit1D", 
  function(object, nmrdata = NULL, peak.units = 'hz', include.bounds = FALSE) {
    callNextMethod()
  })



#========================================================================>
# Plotting  
#========================================================================>



#------------------------------------------------------------------------
#' Plot NMRFit1D object
#'
#' Generates an interactive plot object using the plotly package.
#'
#' Convenience function that generates a graphical representation of the fit.
#' The original data is plotted as a black line, the fit is plotted in red,
#' the baseline is plotted in blue, the residual in red. The fit can be plotted
#' as a composite of all the peaks, or individually.
#'
#' @param x An NMRFit1D object.
#' @param components One of either 'r', 'i', or 'r/i' to include real,
#'                   imaginary or both components. If both components
#'                   are selected, they are displayed in separate subplots.
#' @param sum.lineshapes TRUE to sum together individual peaks.
#' @param sum.baseline TRUE to add the baseline to each fit.
#' @param include.convolution TRUE to apply convolution to the fit, if
#'                            available.
#' @param deconvolve.residual TRUE to remove the effect of convolution
#'                            from the residual, if convolution is not being
#'                            included in the fit.
#' @param apply.phase TRUE to apply the calculated phase to the data.
#' @param nrows Max number of rows to display subplots.
#' @param legend.position One of either 'bottom' or 'side'.
#' @param fit.legend By default, the fit is shown on the legend if 
#'                   sum.lineshapes is TRUE and omitted if sum.lineshapes is
#'                   FALSE. Set fit.legend to TRUE or FALSE to override this
#'                   behaviour.
#'
#' @return A ggplot2 plot.
#'
#' @export
plot.NMRFit1D <- function(x, components = 'r',  apply.phase = TRUE,  
                          sum.lineshapes = TRUE, sum.baseline = TRUE,
                          include.convolution = TRUE,
                          deconvolve.residual = FALSE,
                          nrows = 2, legend.position = 'bottom', 
                          fit.legend = NULL) { 

  # Switching to absolute representation
  nmrfit <- set_peak_units(x, peak.units = 'ppm')
  nmrfit <- set_normalized(nmrfit, normalized = FALSE)

  # Showing legend entries depending on sum.lineshapes
  if ( is.null(fit.legend) ) {
    if ( sum.lineshapes ) fit.legend <- TRUE
    else fit.legend <- FALSE
  }

  # Setting legend options
  if ( legend.position == 'bottom' ) {
    legend.opts <- list(orientation = 'h', xanchor = "center", x = 0.5)
  } else {
    legend.opts <- list()
  }

  # First, performing all calculations
  d <- processed(nmrdata(nmrfit))
  d <- d[order(d$direct.shift), ]

  if ( apply.phase ) {
    phase <- phase(nmrfit)
    d$intensity <- phase_spectrum(d$intensity, phase, degrees = FALSE)
  }

  baseline <- calc_baseline(nmrfit, d$direct.shift)

  # Calculating all line fits
  all.fits <- calc_lineshape(nmrfit, d$direct.shift, include.convolution)

  # Residual is calculated based on sums
  sum.fits <- aggregate(all.fits$intensity, 
                        by = list(all.fits$direct.shift), sum)
  colnames(sum.fits) <- c('direct.shift', 'intensity')
  sum.fits <- sum.fits[order(sum.fits$direct.shift), ]

  # Initializing the plot list
  plots <- list()

  # Checking which components to plot
  plot.r <- grepl('r', components)
  plot.i <- grepl('i', components)

  # Defining function to initialize plot
  f_init <- function(y, color, name) {
    plot_ly(x = d$direct.shift, y = y, color = I(color), 
            name = I(name), type = 'scatter', mode = 'lines',
            legendgroup = 1) %>%
      layout(legend = legend.opts)
  }

  if ( plot.r ) plots$r <- f_init(Re(d$intensity), 'black', 'Real')
  if ( plot.i ) plots$i <- f_init(Im(d$intensity), 'grey', 'Imaginary')

  # Defining function to add traces to initialized plot
  f_add <- function(p, y, color, name, group, showlegend = TRUE) {
    p %>% 
      add_trace(x = d$direct.shift, y = y, color = I(color),
                name = I(name), type = 'scatter', mode = 'lines',
                legendgroup = group, showlegend = showlegend)
  }

  # Adding baseline
  if ( plot.r ) plots$r <- f_add(plots$r, Re(baseline), 'blue', 'Baseline', 2) 
  if ( plot.i ) plots$i <- f_add(plots$i, Im(baseline), 'blue', 'Baseline', 2)

  # Adding residual
  y <- d$intensity - sum.fits$intensity - baseline

  # Removing convolution term from residual if required
  if ( (! include.convolution) && deconvolve.residual ) {
    res.fits <- calc_lineshape(nmrfit, d$direct.shift, TRUE)

    res.fits <- aggregate(res.fits$intensity, 
                          by = list(all.fits$direct.shift), sum)
    colnames(res.fits) <- c('direct.shift', 'intensity')
    res.fits <- res.fits[order(sum.fits$direct.shift), ]
    y <- y - res.fits$intensity + sum.fits$intensity
  }

  if ( plot.r ) plots$r <- f_add(plots$r, Re(y), 'green', 'Residual', 3) 
  if ( plot.i ) plots$i <- f_add(plots$i, Im(y), 'green', 'Residual', 3)

  # Adding fits
  if ( sum.lineshapes ) {
    
    if ( sum.baseline ) y <- sum.fits$intensity + baseline
    else y <- sum.fits$intensity

    if ( plot.r ) plots$r <- f_add(plots$r, Re(y), 'red', 'Fit', 4, fit.legend) 
    if ( plot.i ) plots$i <- f_add(plots$i, Im(y), 'red', 'Fit', 4, fit.legend)

  } else {

    all.fits <- by(all.fits, paste(all.fits$id, all.fits$peak), identity)
    
    # Looping through each peak
    for ( i in 1:length(all.fits) ) {

      d <- all.fits[[i]]

      if ( sum.baseline ) y <- d$intensity + baseline
      else y <- d$intensity

      id <- paste(d$id[1], d$peak[1], sep = '-')

      if ( plot.r ) plots$r <- f_add(plots$r, Re(y), 'red', id, id, fit.legend) 
      if ( plot.i ) plots$i <- f_add(plots$i, Im(y), 'red', id, id, fit.legend)
    }

  }

  if ( length(plots) == 0 ) NULL
  else subplot(plots, shareX = TRUE, shareY = TRUE, 
               nrows = min(length(plots), nrows))
}
