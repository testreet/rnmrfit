# Definition of an intensity class to handle various real/imaginary combinations



#------------------------------------------------------------------------------
#' Container for NMR intensity data
#' 
#' Although 1D intensity can be represented by a single real/imaginary complex
#' pair, there is no built-in representation for 2D (or higher dimension) data
#' that may consist of real-real, real-imaginary, imaginary-real, and
#' imaginary-imaginary components. This class approximates the extension of a
#' complex values to multiple dimensions be defining a vector, where each
#' element corresponds to a specific real-imaginary component. Re() and Im()
#' functions can still be used to obtain pure real and pure imaginary
#' compenents, whereas other specific real/imaginary can be are accesible with
#' the $ function. For 2D data, for example, $rr will output a numeric vector
#' corresponding to the real-real data (equivalent to Re()).
#' 
#' @name NMRIntensity-class
#' @export
NMRIntensity <- setClass("NMRIntensity", contains = "numeric")



#------------------------------------------------------------------------------
#' Container for NMR intensity data
#' 
#' Multiple NMRintensity vectors can be combined in an NMRIntensities list that
#' supports the same Re(), Im() and $ functions as well as conversions to
#' data.frames.
#' 
#' @name NMRIntensities-class
#' @export
NMRIntensities <- setClass("NMRIntensities", contains = "list")



#==============================================================================>
#  Constructors
#==============================================================================>



#------------------------------------------------------------------------------
#' Shortcut for generating component combinations
.complex_names <- function(n) {
  names <- do.call(crossing, rep(list(c('r','i')), n))
  rev(unlist(apply(names, 1, paste, collapse = '')))
}



#------------------------------------------------------------------------------
#' Shortcut for generating ijk suffixes 
.complex_suffixes <- function(n) {
  names <- .complex_names(n)

  # First convert 2nd i's and up
  while ( n > 1 ) {
    logic <- substr(names, n, n) == 'i'
    substr(names, n, n)[logic] <- letters[8+n]
    n = n - 1
  }

  # And then remove the r
  gsub('r', '', names)
}



#------------------------------------------------------------------------------
#' Generate an NMRIntensity or NMRIntensities object
#' 
#' Generates an NMRIntensity or plural NMRIntensities object from a series of
#' vector, matrix, or data.frame inputs.
#' 
#' @param ... A collection of vector, matrix, or data.frame objects. Vectors
#'            must be specified with names whereas the component names for
#'            matrix and data.frame objects are taken from their column names.
#' 
#' @export
nmrintensity <- function(...) {

  args <- list(...)
  arg.names <- names(args)

  # First, if there ar no arguments, return null
  if ( length(args) == 0 ) return(NULL)

  # All arguments must have the same length or number of columns
  n <- lapply(args, nrow)
  n.alternate <- lapply(args, length)
  logic <- unlist(lapply(n, is.null))
  n[logic] <- n.alternate[logic]
  n <- unique(unlist(n))

  err <- 'All input objects must have the same length or row numbers'
  if ( length(n) > 1 ) stop(err)

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
  new.terms <- do.call(cbind, parts)
  new.columns <- colnames(new.terms)
  n.char <- unique(nchar(new.columns))

  err <- paste('All component names must have the same number of characters,',
               'e.g. rr and ri vs r.')
  if ( length(n.char) > 1 ) stop(err)
  
  err <- 'All component names must be composed of "r" and "i" characters.'
  if (! all(grepl('^[ri]+$', new.columns)) ) stop(err)

  err <- 'All component names must be unique.'
  if ( length(new.columns) != length(unique(new.columns)) ) stop(err)

  # Extending
  all.columns <- .complex_names(n.char)

  all.terms <- matrix(0, nrow = n, ncol = 2^n.char)
  colnames(all.terms) <- all.columns

  # Inserting values
  new.terms <- as.matrix(new.terms[, new.columns])
  all.terms[, new.columns] <- new.terms

  # Splitting by individual rows
  all.terms <- split(all.terms, 1:n)

  # Convert all individual rows to NMRIntensity
  f <- function(x) new("NMRIntensity", x)
  all.terms <- lapply(all.terms, f)

  # If the final list length is greater than 1, make it plural
  if ( length(all.terms) > 1) {
    new("NMRIntensities", all.terms)
  } else {
    all.terms[[1]]
  }
}



#==============================================================================>
# Display methods
#==============================================================================>



#------------------------------------------------------------------------------
#' @importFrom pillar type_sum
#' @export
type_sum.NMRIntensity <- function(x) {'ncmplx'}



#------------------------------------------------------------------------------
#' @importFrom pillar type_sum
#' @export
type_sum.NMRIntensities <- function(x) {'ncmplx'}


#------------------------------------------------------------------------------
#' Show an NMRIntensity object
#' 
#' A concise display for NMRIntensity.
#' 
#' @param object An NMRIntensity object.
#' 
#' @export
setMethod("show", "NMRIntensities", 
  function(object) {
    
    cat('An object of NMRIntensities class\n\n')
    print(summary(as.data.frame(object)))
  })



#==============================================================================>
#  Basic getter functions for NMRIntensity
#==============================================================================>



#------------------------------------------------------------------------
#' Get the pure real component of the data
#' 
#' Retrieves the pure real component of NMRIntensity.
#' 
#' @param z NMRIntensity data.
#' 
#' @return A numeric vector.
#' 
#' @export
setMethod("Re", "NMRIntensity", function(z) z[1]) 



#------------------------------------------------------------------------
#' Get the pure imaginary component of the data
#' 
#' Retrieves the pure imaginary component of NMRIntensity.
#' 
#' @param z NMRIntensity data.
#' 
#' @return A numeric vector.
#' 
#' @export
setMethod("Im", "NMRIntensity", function(z) z[length(z)])



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
#' @return A complex vector if two components are matched, a numeric vector
#'         otherwise.
#' 
#' @export
`$.NMRIntensity` <- function(x, name) {

  err <- 'The total number of characters must correspond to components of x.'
  if ( nchar(name) != sqrt(length(x)) ) stop(err)

  i <- grep(name, .complex_names(sqrt(length(x))))

  # Output depends on how many matches there are
  if ( length(i) == 0 ) {
    NULL
  }
  else if ( length(i) == 2 ) {
    re =x[i[1]]
    im = x[i[2]]
    complex(re = re, im = im)
  }
  else {
    x[i]
  }
}



#==============================================================================>
#  Basic getter functions for NMRIntensities
#==============================================================================>



#------------------------------------------------------------------------
#' Get the pure real component of the data
#' 
#' Retrieves the pure real component of NMRIntensities.
#' 
#' @param z NMRIntensities data.
#' 
#' @return A numeric vector.
#' 
#' @export
setMethod("Re", "NMRIntensities", function(z) unlist(lapply(z, Re))) 



#------------------------------------------------------------------------
#' Get the pure imaginary component of the data
#' 
#' Retrieves the pure imaginary component of NMRIntensities.
#' 
#' @param z NMRIntensities data.
#' 
#' @return A numeric vector.
#' 
#' @export
setMethod("Im", "NMRIntensities", function(z) unlist(lapply(z, Im)))



#------------------------------------------------------------------------
#' Convert NMRIntensities object into a full data.frame
#' 
#' Since an NMRIntensities object is essentially a data.frame broken down by
#' rows so that it can used as a column inside other data.frames, it's
#' sometimes convenient to access the data as a proper data.frame object.
#' 
#' @param x NMRIntensities object.
#' 
#' @return A data.frame with columns labelled as real/imaginary components.
#' 
#' @export
as.data.frame.NMRIntensities <- function(x) {

  d <- do.call(rbind, x)
  colnames(d) <- .complex_names(sqrt(ncol(d)))
  d

}



#------------------------------------------------------------------------
#' Get a specific real/imaginary component
#' 
#' Specify the components of interest using a series of "r" and "i" characters,
#' such as "ri" or "ir". Internally, the name is passed into grep() so "r.",
#' for example, will be matched as both "rr" and "ri".
#' 
#' @param x NMRIntensities object.
#' @param name A character name for the component. This can be an exact match
#'             such as "ri" or a partial match such as "r.".
#' 
#' @return A numeric vectors if a single component is matched, a complex vector
#'         if two components are matched, and a data.frame otherwise.
#' 
#' @export
`$.NMRIntensities` <- function(x, name) {

  err <- 'The total number of characters must correspond to components of x.'
  if ( nchar(name) != sqrt(length(x[[1]])) ) stop(err)

  i <- grep(name, .complex_names(length(x)))

  # Output depends on how many matches there are
  if ( length(i) == 0 ) {
    NULL
  }
  else if ( length(i) == 1 ) {
    unlist(lapply(x, `[`, i))
  }
  else if ( length(i) == 2 ) {
    re = unlist(lapply(x, `[`, i[1]))
    im = unlist(lapply(x, `[`, i[2]))
    complex(re = re, im = im)
  }
  else {
    as.data.frame(x)[,i]
  }
}
