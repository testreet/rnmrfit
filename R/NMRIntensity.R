# Definition of an intensity class to handle various real/imaginary combinations



#------------------------------------------------------------------------------
#' Container for NMR intensity data
#' 
#' Although 1D intensity can be represented by a single real/imaginary complex
#' pair, there is no built-in representation for 2D (or higher dimension) data
#' that may consist of real-real, real-imaginary, imaginary-real, and
#' imaginary-imaginary components. This class serves as a simple container list
#' that can be used to retrieve specific sets of values as required. Re() and
#' Im() can be used to obtain pure real and pure imaginary compenents, whereas
#' other specific real/imaginary can be are accesible via indexing. For 2D
#' data, for example, $rr will output a numeric vector corresponding to the
#' real-real data (equivalent to Re()).
#' 
#' @slot columns Defines the real/imaginary composition of each elemnt as a
#'               series of "r" and "i" values.
#' 
#' @name NMRIntensity-class
#' @export
NMRIntensity <- setClass("NMRIntensity", 
                         slots = c(columns = "character"),
                         contains = "list")



#==============================================================================>
#  Constructor
#==============================================================================>

#------------------------------------------------------------------------------
#' Generate an NMRIntensity object
#' 
#' Generates an NMRIntensity object from a series of vector, matrix, or
#' data.frame object.
#' 
#' @param ... A collection of vector, matrix, or data.frame objects. Vectors
#'            must be specified with names whereas the component names for
#'            matrix and data.frame objects are taken from their column names.
#' 
#' @export
nmrintensity <- function(...) {

  args <- list(...)
  arg.names <- names(args)

  # If there ar eno arguments, return null
  if ( length(args) == 0 ) return(NULL)

  # If there are no names, generate empty spaces to facilitate further parsing
  if ( length(arg.names) == 0 ) arg.names <- rep('', length(args))

  # Initializing a parts list to be combined later
  parts <- list()

  # First, handling the unnamed arguments
  logic.unnamed <- arg.names == ''

  if ( sum(logic.unnamed) > 0 ) {
    unnamed <- args[logic.unnamed]

    # Unnamed arguments should only be matrix or data.frame
    classes <- lapply(unnamed, class)
    f.class <- function(x) any(x %in% c('matrix', 'data.frame'))
    err <- 'Unnamed arguments must have either "matrix" or "data.frame" class.'
    if ( any(! unlist(lapply(classes, f.class))) ) stop(err)

    unnamed.part <- do.call(cbind, unnamed)
    parts <- c(parts, list(unnamed.part))
  }

  # Then the named arguments
  logic.named <- ! logic.unnamed

  if ( sum(logic.named) > 0 ) {
    named <- args[logic.named]

    # Named arguments should all be numeric
    classes <- lapply(named, class)
    f.class <- function(x) any(x %in% c('numeric'))
    err <- 'Named arguments must have a "numeric" class.'
    if ( any(! unlist(lapply(classes, f.class))) ) stop(err)

    named.part <- do.call(cbind, named)
    colnames(named.part) <- names(named)
    parts <- c(parts, list(named.part))
  }

  # Combining and checking column names
  combined <- do.call(cbind, parts)
  columns <- colnames(combined)
  
  err <- 'All component names must be composed of "r" and "i" characters.'
  if (! all(grepl('^[ri]+$', columns)) ) stop(err)

  err <- 'All component names must be unique'
  if ( length(columns) != length(unique(columns)) ) stop(err)

  # Sorting and splitting by row
  columns <- rev(sort(columns))

  combined <- as.matrix(combined[, columns])
  combined <- split(combined, 1:nrow(combined))

  # Outputting NMRIntensity object
  new("NMRIntensity", combined, columns = columns)
}



#==============================================================================>
#  Basic setter and getter functions
#==============================================================================>



#------------------------------------------------------------------------
#' Get the pure real component of the data
#' 
#' Retrieves either the pure real component or the component with the most real
#' terms when dealing with segmented data.
#' 
#' @param z NMRIntensity data.
#' 
#' @return A numeric vector.
#' 
#' @export
setMethod("Re", "NMRIntensity", 
  function(z) {
    out <- unlist(lapply(z, function(x) {x[1]}))
    names(out) <- NULL
    out
  })



#------------------------------------------------------------------------
#' Get the pure imaginary component of the data
#' 
#' Retrieves either the pure imaginary component or the component with the most
#' imaginary terms when dealing with segmented data.
#' 
#' @param z NMRIntensity data.
#' 
#' @return A numeric vector.
#' 
#' @export
setMethod("Im", "NMRIntensity", 
  function(z) {
    out <- unlist(lapply(z, function(x) {x[length(x)]}))
    names(out) <- NULL
    out
  })



#------------------------------------------------------------------------
#' Get a specific real/imaginary component
#' 
#' Specify the components of interest using a series of "r" and "i" characters,
#' such as "ri" or "ir". Internally, the name is passed into grep() so "r.",
#' for example, will be matched as both "rr" and "ri".
#' 
#' @param x NMRIntensity object.
#' @param name A character name for the component. This can be an exact match
#'             such as "ri" or a partial match such as "r.".
#' 
#' @return A numeric vector if only a single component is matched, a complex
#'         vector if two components are matched, and an NMRIntensity object if
#'         three or more components are matched.
#' 
#' @export
`$.NMRIntensity` <- function(x, name) {

  # First, check if name is one of the columns
  if ( name %in% x@columns ) {
    i <- which(name == x@columns)
    out <- lapply(x, function(x) {x[i]})
    names(out) <- NULL
    unlist(out)
  }
  # Otherwise try and match
  else {
    i <- grep(name, x@columns)

    # Output depends on how many matches there are
    if ( length(i) == 0 ) {
      NULL
    }
    else if ( length(i) == 1 ) {
      out <- lapply(x, function(x) {x[i]})
      names(out) <- NULL
      unlist(out)
    }
    else if ( length(i) == 2 ) {
      re = unlist(lapply(x, function(x) {x[i[1]]}))
      im = unlist(lapply(x, function(x) {x[i[2]]}))
      out <- complex(re = re, im = im)
      names(out) <- NULL
      out
    }
    else {
      # First filter the list
      d.list <- lapply(x, function(x) {x[i]})
      d.mat <- do.call(rbind, d.list)
      colnames(d.mat) <- x@columns
      nmrintensity(d.mat)
    }
  }
}
