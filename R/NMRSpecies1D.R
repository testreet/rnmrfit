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
    id = 'character',
    resonances = 'list',
    connections = 'data.frame',
    connections.leeway = 'numeric'
  ),
  prototype = prototype(
    id = 'species',
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
  # Checking name
  if ( length(id) != 1 ) {
    valid <- FALSE
    new.err <- '"name" must be a character vector of length 1.'
    err <- c(err, new.err)
  }

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
    if (! identical(colnames(connections), valid.columns) ) {
      valid <- FALSE
      new.err <- sprintf('"connections" must have the following columns: %s',
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
#'                   these names will be use to replace resonance ids.
#' @param areas A vector of areas corresponding to the expected areas of the
#'              resonances. Set to NULL by default, signifying no fixed
#'              constraings.
#' @param id A string specifying species name. If left empty, a name is
#'           automatically generated from the resonance names.
#' @param connections.leeway A value specifying how tightly enforced the
#'                           connection constraints on resonance areas should
#'                           be. E.g. a value of = 0 specifies that the area
#'                           ratios are exact, whereas 0.1 specifies that the
#'                           areas of the resonances may differ by +/- 10
#'                           percent from the specified ratios.
#' @param ... Options passed to nmrresonance_1d if resonances are being
#'            converted from character/numeric vectors. See ?nmrresonance_1d for
#'            more details.
#' 
#' @return An NMRSpecies1D object.
#' 
#' @export
nmrspecies_1d <- function(resonances, areas = NULL, id = NULL, 
                          connections.leeway = 0, ...) {

  #---------------------------------------
  # Generating list of resonances
  resonances.list <- list()

  for (i in 1:length(resonances)) {

    resonance <- resonances[[i]]

    if ( class(resonance) == 'NMRSpecies1D' ) {
      err <- paste("An NMRSpecies1D can't be constructed from other",
                   "NMRSpecies1D objects. Use resonances() to first extract",
                   "the resonance list before creating a new object.")
      stop(err)
    }
    # If the object is already an NMRResonance1D object, add it directly
    else if ( class(resonance) == 'NMRResonance1D' ) {
      resonances.list <- c(resonances.list, resonance)
    }
    # Otherwise, feed it into the nmrresonance_1d constructor
    else {
      resonances.list <- c(resonances.list, nmrresonance_1d(resonance, ...))
    }
        
    # Modifying id if provided
    resonance.id <- names(resonances)[i]
    if (! is.null(resonance.id) ) id(resonances.list[[i]]) <- resonance.id
  }

  #---------------------------------------
  # Defining connections if areas provided

  # Fetching resonance ids from list
  valid.ids <- unlist(lapply(resonances.list, function(o) o@id))

  if (! is.null(areas) ) {
    # Checking that areas correspond to resonances
    err <- paste('Either "areas" vector must have names or the length of',
                 '"areas" vector must match length of resonances.')
    if (! is.null(names(areas)) ) ids <- names(areas)
    else if ( length(areas) == length(valid.ids) ) ids <- valid.ids
    else stop(err) 

    # Checking that area names are valid
    err <- 'Names of "areas" vector must be valid resonance ids.'
    if ( any(! ids %in% valid.ids) ) stop(err)

    # Checking length
    err <- '"areas" vector must be of length 2 or more to add constraints.'
    if ( length(ids) < 2 ) stop(err)

    # Generating connections data frame
    n <- length(areas)
    index.1 <- 1:(n - 1)
    index.2 <- 2:n
    connections <- data.frame(resonance.1 = ids[index.1], 
                              resonance.2 = ids[index.2],
                              area.ratio = areas[index.2]/areas[index.1])
    rownames(connections) <- 1:nrow(connections) 
  } 
  else {
    connections = data.frame()
  }

  # Generating id if it doesn't exist
  if ( is.null(id) ) id <- paste(valid.ids, collapse = '-')

  new('NMRSpecies1D', id = id, resonances = resonances.list, 
                      connections = connections, 
                      connections.leeway = connections.leeway)
}



#==============================================================================>
#  Display function
#==============================================================================>



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
# Id

#' @rdname id
#' @export
setMethod("id", "NMRSpecies1D", 
  function(object) object@id
  )

#' @rdname id-set
#' @export
setReplaceMethod("id", "NMRSpecies1D",
  function(object, value) {
    object@id <- as.character(value)
    validObject(object)
    object 
  })



#------------------------------------------------------------------------------
# Peaks

#' @rdname peaks
#' @export
setMethod("peaks", "NMRSpecies1D", 
  function(object, include.id = FALSE) {
    peaks.list <- lapply(object@resonances, peaks, include.id = TRUE)
    peaks <- do.call(rbind, peaks.list)
    if ( include.id ) cbind(species = object@id, peaks)
    else peaks
  })



#------------------------------------------------------------------------------
# Couplings

#' @rdname couplings
#' @export
setMethod("couplings", "NMRSpecies1D", 
  function(object) object@couplings)

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

#' @rdname bounds
#' @export
setMethod("bounds", "NMRSpecies1D", 
  function(object) object@bounds)

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
