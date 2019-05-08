# Definition of a class structure for lineshape fitting.



#==============================================================================>
#  NMRFit1D
#==============================================================================>



#------------------------------------------------------------------------------
#' Definition of an NMR fit.
#' 
#' Essentially, this class is used to combine an NMRData1D object with multiple
#' NMRSpecies1D objects while also defining baseline and phase corrections
#' terms. There are two primary methods associated with this class: fit() and
#' update(). fit() takes input data and peak definitions, applies a nonlinear
#' least squares fit, and outputs best-fit curves into a slot called results.
#' update() overwrites initial peak definitions with the results of the fit to
#' generate a new set of initial values for further rounds of fitting.
#' 
#' @slot species A list of NMRSpecies1D objects.
#' @slot nmrdata An NMRData1D object used to fit the peaks.
#' @slot knots A vector of chemical shifts corresponds to the internal knots of
#'             the baseline term (two boundary knots are always included).
#' @slot baseline A vectors of complex numeric values representing the baseline.
#'                The actual baseline is generated from a basis spline so these
#'                values roughly represent the intensity values of the baseline
#'                based on the chemical shift locations of the knots. Baseline
#'                degree and default number of internal knots can be set using
#'                nmrsession_1d(baseline = list(degree = 3, n.knots = 0)), with
#'                the baseline vector length equal to the degree + n.knots.
#' @slot phase A number vectors of phase correction terms, starting from 0 order
#'             and up. The default number of phase terms can be set using
#'             nmrsession_1d(n.phase = 1).
#' @slot results A list of results that contains update species list as well as
#'               baseline and phase vectors.
#' 
#' @name NMRFit1D-class
#' @export
NMRFit1D <- setClass("NMRFit1D",
  slots = c(
    species = 'list',
    nmrdata = 'NMRData1D',
    knots = 'numeric',
    baseline = 'complex',
    phase = 'numeric'
  ),
  prototype = prototype(
    species = list(),
    knots = numeric(0), 
    baseline = complex(re = rep(0, 3), im = rep(0, 3)),
    phase = c(0)
  )
)



#==============================================================================>
#  Validation methods
#==============================================================================>



#------------------------------------------------------------------------------
#' NMRFit1D validity test
#'
validNMRFit1D <- function(object) {

  species <- object@species
  nmrdata <- object@nmrdata
  knots <- object@knots
  baseline <- object@baseline
  phase <- object@phase

  # The result list is not validated as it should be relatively safe to assume
  # that this slot would not be touched directly.

  valid <- TRUE
  err <- c()

  #---------------------------------------
  # Checking that all species list items are valid
  for ( specie in species ) {
    logic1 <- class(specie) != 'NMRSpecies1D'
    logic2 <- ! validObject(specie)
    if ( logic1 || logic2 ) {
      valid <- FALSE
      new.err <- paste('All elements of "species" list must be valid',
                       'NMRSpecies1D objects.')
      err <- c(err, new.err)
    }
  }

  #---------------------------------------
  # Checking nmrdata
  if ( (class(nmrdata) != 'NMRData1D') || (! validObject(nmrdata))  ) {

      valid <- FALSE
      new.err <- '"nmrdata" must be a valid NMRData1D object.'
      err <- c(err, new.err)

  }

  #---------------------------------------
  # Checking basline length 
  if ( length(baseline) <= length(knots)  ) {

      valid <- FALSE
      new.err <- paste('"baseline" vector length must be greater than the',
                       '"knots" vector length.')
      err <- c(err, new.err)

  }

  #---------------------------------------
  # Checking phase length 
  if ( length(phase) > 2  ) {

      wrn <- paste('Although "phase" slot lengths of greater than 2 are',
                   'supported, second order phase corrections and above are',
                   'highly unlikely. Results may be misleading.')
      warning(warn)

  }

  #---------------------------------------
  # Output
  if (valid) TRUE
  else err
}

# Add the extended validity testing
setValidity("NMRFit1D", validNMRFit1D)



#==============================================================================>
#  Constructor for initialization
#==============================================================================>



#------------------------------------------------------------------------------
#' Generate an NMRFit1D object
#' 
#' Generates an NMRFit1D object based on a list of NMRSpecies1D objects or
#' other objects that can be converted to NMRSpecies1D objects. See
#' ?nmrresonance_1d and ?nmrspecies_1d for more details about this conversion.
#' Apart from providing peak definitions (via species argument) and data (via
#' nmrdata argument), the main fit decisions are related to baseline and phase
#' correction, as well as how to deal with peaks that are defined outside the
#' range of the data. This latter decision is broken down into two arguments:
#' exclusion.level and exclusion.notification. The exclusion.level parameter
#' determines which part of the overall species to exclude if any of its peaks
#' fall outside the data range: either 'species' for whole species, 'resonance'
#' for just a subset of the species and 'peak' to ignore resonance/species
#' blocks and exclude by specific peak alone. The exclusion.notification
#' parameter determines how to respond when peaks are found to be outside the
#' data range: either 'nothing' to give no notice, 'warning' to issue a
#' warning, and 'error' to issue an error. A function can also be provided with
#' a single input, which corresponds to the message text, which can then be
#' parsed or saved to a log file as desired.
#' 
#' @param species A list of NMRSpecies1D objects or other objects that can be
#'                converted to NMRSpecies1D objects. See ?nmrresonance_1d and
#'                ?nmrspecies_1d for more details about this conversion. If list
#'                elements are named, these names will be use to replace
#'                resonance ids.
#' @param nmrdata An NMRData1D object used to fit the supplied peaks.
#'                automatically generated from the resonance names.
#' @param baseline.order An integer specifying the order of the baseline spline
#'                       function. Note the the internal B-spline implementation
#'                       requires an order of 1 or greater. You can use an order
#'                       of -1 to disable baseline correction.
#' @param n.knots An integer specifying the number of internal knots to use. The
#'                specific position of these knots can be modified later using
#'                knots() function. To modify initial values of the baseline,
#'                use baseline() function.
#' @param phase.order An integer specifying the order of the phase correction
#'                    polynomial. Only 0 and 1st order terms are typically used
#'                    in practice, but a higher order correction is possible.
#' @param delay.fit FALSE to immediately run least squares optimization after
#'                  the NMRFit1D object is initialized, TRUE to skip the
#'                  optimization, enabling more customization. The fit can be
#'                  run manually using fit().
#' @param exclusion.level A string specifying what to do when peaks are found to
#'                        fall outside of the data range: either 'species' to
#'                        exclude the whole species to which the offending peak
#'                        belongs, 'resonance' to exclude the resonance to which
#'                        the offending peak belongs, or 'peak' to exclude just
#'                        the peak itself. resonance/species blocks and exclude
#'                        by specific peak alone.
#' @param exclusion.notification A string specifying how to notify the user when
#'                               peaks are found to be outside the data range:
#'                               either 'none' to give no notice, 'warning' to
#'                               issue a warning, and 'error' to issue an error.
#'                               A function can also be provided with a single
#'                               input, which corresponds to the message text,
#'                               which can then be parsed or saved to a log file
#'                               as desired.
#' @param ... Options passed to nmrspecies_1d if conversion has to be performed.
#'            See ?nmrspecies_1d for more details.
#' 
#' @return An NMRFit1D object.
#' 
#' @export
nmrfit_1d <- function(species, nmrdata,
                      baseline.order = nmrsession_1d$baseline$order,
                      n.knots = nmrsession_1d$baseline$n.knots, 
                      phase.order = nmrsession_1d$phase$order, 
                      delay.fit = FALSE,
                      exclusion.level = 'resonance', 
                      exclusion.notification = 'warning', ...) {

  #---------------------------------------
  # Generating list of species 


  # If the species object is just a single species, add it to list
  if ( class(species) == 'NMRSpecies1D' ) {
    species.list <- list(species)
  }
  # Otherwise, loop through
  else {

    species.list <- list()
  
    for (i in 1:length(species)) {

      specie <- species[[i]]

      # If the object is already an NMRSpecies1D object, add it directly
      if ( class(specie) == 'NMRSpecies1D' ) {
        species.list <- c(species.list, specie)
      }
      # Otherwise, feed it into the nmrspecies_1d constructor
      else {
        species.list <- c(species.list, nmrspecies_1d(specie, ...))
      }
          
      # Modifying id if provided
      specie.id <- names(species)[i]
      if (! is.null(specie.id) ) id(species.list[[i]]) <- specie.id
    }
  }

  #---------------------------------------
  # Checking nmrdata

  if ( (class(nmrdata) != 'NMRData1D') || (! validObject(nmrdata)) ) {
    err <- '"nmrdata" must be a valid NMRData1D object.'
    stop(err)
  }

  #---------------------------------------
  # Baseline and phase

  # The initial value for the baseline is just the median of nmrdata intensity
  if ( baseline.order == -1 ) {
    baseline <- complex(0)
    knots <- numeric(0)
  }
  else {
    n <- baseline.order + n.knots
    baseline <- complex(re = rep(median(Re(nmrdata@processed$intensity)), n),
                        im = rep(median(Im(nmrdata@processed$intensity)), n))

    # Knots are initialized to fall evenly between the chemical shift data
    direct.shift <- range(nmrdata@processed$direct.shift)
    knots <- seq(direct.shift[1], direct.shift[2], length.out = n.knots)
  }

  # The initial value for the phase is always 0
  if ( phase.order == -1 ) phase <- numeric(0)
  else phase <- rep(0, phase.order)

  #---------------------------------------
  # Resulting fit object
  out <- new('NMRFit1D', species = species.list, nmrdata = nmrdata,
                         knots = knots, baseline = baseline, phase = phase)

  # If the fit is delayed, then return current object, otherwise fun fit first
  if ( delay.fit ) out
  else fit(out, exclusion.level = exclusion.level, 
           exclusion.notification = exclusion.notification)

}



#==============================================================================>
#  Display function
#==============================================================================>



#------------------------------------------------------------------------------
#' Display NMRFit1D object
#'
#' Display a quick summary of resonance parameters.
#'
#' @export
setMethod("show", "NMRFit1D", 
  function(object) {

    # Generating compiled data frames
    peaks <- peaks(object)
    bounds <- bounds(object)
    couplings <- couplings(object)

    cat('An object of NMRFit1D class\n\n')

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

#' @rdname peaks
#' @export
setMethod("peaks", "NMRFit1D", 
  function(object) {
    peaks.list <- lapply(object@species, peaks, include.id = TRUE)
    do.call(rbind, peaks.list)
  })



#------------------------------------------------------------------------------
# Couplings

#' @rdname couplings
#' @export
setMethod("couplings", "NMRFit1D", 
  function(object, include.id = FALSE) {
    couplings.list <- lapply(object@species, couplings, include.id = TRUE)
    do.call(rbind, couplings.list)
  })



#------------------------------------------------------------------------------
# Bounds

#' @rdname bounds
#' @export
setMethod("bounds", "NMRFit1D", 
  function(object) {
    f <- function(o, sublist) bounds(o, include.id = TRUE)[[sublist]]
    lower.list <- lapply(object@species, f, sublist = 'lower')
    upper.list <- lapply(object@species, f, sublist = 'upper')

    lower <- do.call(rbind, lower.list)
    upper <- do.call(rbind, upper.list)

    list(lower = lower, upper = upper)
  })



#==============================================================================>
#  Bounds
#==============================================================================>


#------------------------------------------------------------------------------
#' @rdname set_general_bounds
#' @export
setMethod("set_general_bounds", "NMRFit1D",
  function(object, ...) {
    object@species <- lapply(object@speies, set_general_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_offset_bounds
#' @export
setMethod("set_offset_bounds", "NMRFit1D",
  function(object, ...) {
    object@species <- lapply(object@species, set_offset_bounds, ...)
    object
  })



#------------------------------------------------------------------------------
#' @rdname set_conservative_bounds
#' @export
setMethod("set_conservative_bounds", "NMRFit1D",
  function(object, ...) { 
    object@species <- lapply(object@species, set_conservative_bounds, ...)
    object
  })



#========================================================================>
#  Lineshape and area calculations
#========================================================================>



#------------------------------------------------------------------------
#' @rdname f_lineshape
#' @export
setMethod("f_lineshape", "NMRFit1D", 
          getMethod("f_lineshape", "NMRResonance1D"))



#------------------------------------------------------------------------
#' @rdname values
#' @export
setMethod("values", "NMRFit1D",
           getMethod("values", "NMRResonance1D"))



#------------------------------------------------------------------------
#' @rdname areas 
#' @export
setMethod("areas", "NMRFit1D",
           getMethod("areas", "NMRResonance1D"))
