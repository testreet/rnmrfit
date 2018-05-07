# Definition of a class structure for 1D NMR data.

#------------------------------------------------------------------------
#' A class combining NMR data and scan parameters
#'
#' Although it's likely to change in the near future, this class currently
#' serves to combine processed 1D data from the pdata folder of an experiment
#' with the acqus and procs parameter files.
#'
#' @slot processed 1r/1i data stored as a data.frame with direct.shift and
#'                 intensity as columns.
#' @slot acqus A list of acqus parameters.
#' @slot procs A list of procs parameters.
#' @slot product A vector representation of an apodization function that was 
#'               used to generate the intensity data (in this context, 
#'               zero-filling is considered to be the equivalent of 
#'               apodization by a step function). To facilitate convolution, 
#'               the length of this vector is set to 2*n - 1 where n is the 
#'               number of spectrum points.
#' @slot convolution The Fourier transform of the product slot, suitable
#'                   for convolution with the spectrum intensity data.
#'
#' @name NMRData1D-class
#' @export
NMRData1D <- setClass("NMRData1D", 
                      contains = "NMRData",
                      slots = c(product = 'numeric',
                                convolution = 'numeric'))

# Validity testing consists of simply checking the processed data.frame columns
validNMRData1D <- function(object) {

  valid <- TRUE
  msg <- c()

  processed <- object@processed
  product <- object@product
  convolution <- object@convolution

  # Checking column names
  valid.columns <- c('direct.shift', 'intensity')
  new.msg <- sprintf('Processed data must have two columns: %s',
                     paste(valid.columns, collapse = ', '))

  if (! identical(valid.columns, colnames(processed)) ) {
    valid <- FALSE
    msg <- c(msg, new.msg)
  }

  # Checking convolution length
  new.msg <- paste('If not empty, "convolution" and "product" slots must',
                   'have 2*n - 1 elements where n is the number of points in',
                   'the "processed" data.frame')

  if ( (length(product) > 0) || (length(convolution) > 0) ) {

    n <- nrow(processed)
    if ( (length(product) != (2*n-1)) || (length(convolution) != (2*n-1)) ) {
      valid <- FALSE
      msg <- c(msg, new.msg)
    }

  }

  if ( valid ) TRUE
  else msg
}

setValidity("NMRData1D", validNMRData1D)

#------------------------------------------------------------------------
#' Constructor for generating an NMRData1D object
#'
#' Loads data in Bruker folder and uses it to initialize an NMRData1D object.
#'
#' @param path Path to a Bruker scan directory.
#' @param number Processing file number. Defaults to the smallest
#'               number in the pdata directory.
#' 
#' @return An NMRData1D object containing the 1r/1i processed data as well
#'         as the procs and acqus parameters.
#'
#' @export
nmrdata_1d <- function(path, number = NA) {

  # First, loading procs parameters
  procs <- read_procs_1d(path, number)

  # Using the procs file to load the processed data
  processed <- read_processed_1d(path, procs, number)

  # Finally, loading the general acquisition parameters
  acqus <- read_acqus_1d(path)

  # Returning class object
  new("NMRData1D", processed = processed, procs = procs, acqus = acqus)
}



#========================================================================>
# Basic conversions and slot access functions (inherited from NMRData)
#========================================================================>

#------------------------------------------------------------------------

#' @rdname processed 
#' @export
setMethod("processed", "NMRData1D", 
  function(object) callNextMethod())

#' @rdname processed-set
#' @export
setReplaceMethod("processed", "NMRData1D", 
  function(object) callNextMethod())

#------------------------------------------------------------------------

#' @rdname procs 
#' @export
setMethod("procs", "NMRData1D", 
  function(object) callNextMethod())

#' @rdname procs-set
#' @export
setReplaceMethod("procs", "NMRData1D", 
  function(object) callNextMethod())

#------------------------------------------------------------------------

#' @rdname acqus 
#' @export
setMethod("acqus", "NMRData1D", 
  function(object) callNextMethod())

#' @rdname acqus-set
#' @export
setReplaceMethod("acqus", "NMRData1D", 
  function(object) callNextMethod())



#========================================================================>
# Defining list and data.frame like behaviour
#========================================================================>

#' @rdname as.list.NMRData
#' @export
setMethod("as.list", "NMRData1D",  
  function(x) {
    callNextMethod()
  })

#------------------------------------------------------------------------

#' @rdname as.data.frame.NMRData
#' @export
setMethod("as.data.frame", "NMRData1D",  
  function(x) {
    callNextMethod()
  })



#========================================================================>
# Non-inherited methods
#========================================================================>

#------------------------------------------------------------------------
#' Display NMRData1D object
#'
#' Display a quick summary of the processed data.frame.
#'
#' @export
setMethod("show", "NMRData1D", 
  function(object) {
    cat('An object of NMRData1D class\n')
    print(summary(object@processed))
    cat('\n')
  })

#------------------------------------------------------------------------
#' Set convolution applied to processed data
#'
#' Although it's possible to fit the processed data based on its final
#' appearance, it may sometimes be more useful to try and fit the underlying
#' lineshape with an explicit consideration of the convolution that it has
#' undergone. This convolution includes any apodization as well as the process 
#' of zero-filling, which can be considered as the convolution of a 
#' hypothetically longer experiment by a step function. By default, the 
#' applied apodization is read from the procs parameters, but these can
#' be overwritten by either a list of parameters or a custom function.
#'
#' @param object An NMRScaffold1D or NMRScaffold2D object.
#' @param trim Number of point to truncate from the end of the fid's 
#'             acquisition time (not including zero fill). It's common to 
#'             circle shift Bruker's group delay artefact from the start of 
#'             a collected fid to the last 50-100 points. While most 
#'             apodization methods will passively get rid of or diminish the
#'             influence of these points, they may have a stronger impact
#'             with minimal or no apodization.
#' @param param A list of lists combining apodization names and parameters,
#'              (see apodize_signal) for a full list of apodizations and their
#'              parameters. E.g. list('exponential'=list('lb'=1)).
#' @param f A function f(n) capable of generating the convolution vector
#'          of length n (ranging from 0 to the acquisition time). If provided, 
#'          f overrides the param list.
#'
#' @return An NMRData1D object with modified convolution vector. 
#'
#' @name set_convolution
#' @export
setGeneric("set_convolution", 
  function(object, trim = 0, param = list(), f = NULL) {
    standardGeneric("set_convolution")
  })

#' @rdname set_convolution
#' @export
setMethod("set_convolution", "NMRData1D",
  function(object, trim = 0, param = list(), f = NULL) {

    processed <- object@processed

    # Erroring out for unevenly sampled data
    differences <- sort(unique(diff(processed$direct.shift)))
    if ( (length(differences) > 1) && 
         any(abs(diff(differences)) > 1e-10) ) {
      msg <- 'Convolution vectors can only be set for evenly sampled data.'
      stop(msg, call. = TRUE)
    }

    si <- object@procs[['si']]
    td <- object@acqus[['td']]
    sw <- object@acqus[['sw.h']]

    # Correcting group delay
    if ( trim > 0 ) {
      processed <- processed[order(processed$direct.shift), ]
      signal <- ifft(ifftshift(processed$intensity))
      index <- 1:(round(td/2) - trim)
      signal[-index] <- 0
      processed$intensity <- fftshift(fft(signal))

      object@processed <- processed
    }

    # Defining what fraction of overall points are signal (vs. zero-fill)
    frac <- (td/2 - trim)/si
    aq <- 1/(sw)*(td/2 - trim)

    # If param is empty, building up the list from the procs parameters
    if ( length(param) == 0 ) {

      lb <- object@procs[['lb']]
      gb <- object@procs[['gb']]
      ssb <- object@procs[['ssb']]
      tm1 <- object@procs[['tm1']]
      tm2 <- object@procs[['tm2']]

      if ( lb != 0 ) {
        if ( lb > 0 ) param$exponential <- list(lb = lb)
        else param$gaussian <- list(lb = lb, gb = gb)
      }

      if ( ssb != 0 ) {
        msg <- paste('Assuming that "ssb" parameter pertains to sine rather',
                     'than qsin apodization. Use "param" to override.')
        warning(msg, call. = TRUE)

        param$sine <- list(ssb = ssb)
      }

      if ( (tm1 !=0 ) || (tm2 != 0) ) {
        param$trapezoid <- list(tm1 = tm1, tm2 = tm2)
      }
    }

    # If f is not defined, building f from param list
    if ( is.null(f) ) {

      valid.names <- c('exponential', 'gaussian', 'sine', 'sine2', 'trapezoid')
      if ( any(! names(param) %in% valid.names) ) {
        msg <- sprintf('Names of "param" list must only include: %s',
                       paste(valid.names, collapse = ', '))
        stop(msg, call. = TRUE)
      }

      f <- function(n) {
        
        direct.time <- seq(0, aq, length.out = n)
        signal <- rep(1, n)

        # Looping through the apodizations
        for ( method in names(param) ) {
          apod.args <- list(direct.time = direct.time, fid = signal)
          apod.args <- c(apod.args, method = method, as.list(param[[method]]))

          signal <- do.call(apodize_signal, apod.args)
        }

        signal
      }
    }

    # Filling in values
    n.out <- 2*nrow(processed) - 1
    n.signal <- round(frac*n.out)
    n.fill <- n.out - n.signal
    product <- c(f(n.signal), rep(0, n.fill))
    object@product <- product/sum(abs(product))
    object@convolution <- Re(fftshift(fft(object@product)))

    object
  })

#------------------------------------------------------------------------
#' Filter 1D chemical shift data by selecting specific regions
#'
#' Filters processed data to include only that which is contained between
#' a set of lower and upper bounds on chemical shift.
#'
#' @param object An NMRData1D object.
#' @param lower A lower bound for chemical shift (in the direct dimension).
#' @param upper An upper bound for chemical shift (in the direct dimension).
#' @param round.up True to round up the total number of points to a power of 2.
#'              This is useful to ensure rapid convolution using the fft when
#'              convolution is included in the fit.
#'
#' @return An NMRData1D object with filtered processed data.
#'
#' @name filter_1d
#' @export
setGeneric("filter_1d", 
  function(object, lower = NULL, upper = NULL, round.up = FALSE, ...) {
    standardGeneric("filter_1d")
  })

#' @rdname filter_1d 
#' @export
setMethod("filter_1d", "NMRData1D", 
  function(object, lower, upper, round.up) {

    x <- object@processed$direct.shift
  
    if ( is.null(lower) ) lower <- min(x)
    if ( is.null(upper) ) upper <- max(x)

    # Checking bounds
    if ( lower > upper ) {
      msg <- 'The lower bound must be smaller than or equal to the upper bound.'
      stop(msg)
    }

    if ( length(c(lower, upper)) != 2 ) {
      msg <- 'Multiple range definition is not currently supported.'
      stop(msg)
    }

    # If x happens to be empty, don't do anything else
    if ( length(x) == 0 ) return(object)

    # Warning if bounds are outside range of chemical shifts 
    if ( (lower > max(x)) | (upper < min(x)) ) {
      msg <- 'Some bounds are entirely out of chemical shift range'
      warning(msg, call. = FALSE)
    }

    # Generating index
    index <- (x > lower) & (x < upper)

    # Rounding if necessary
    if ( round.up ) {
      n <- sum(index)
      pwr = log(n)/log(2)

      n.new <- 2^ceiling(pwr)
      n.extra <- n.new - n

      # Figuring out if there is space to add new points equally on either side
      rle.index <- rle(index)


      # If there are three entries, they are FALSE, TRUE, FALSE
      if ( length(rle.index$value) == 3 ) {
        left.count <- rle.index$lengths[1]
        i.left <- left.count + 1

        right.count <- rle.index$lengths[3]
        i.right <- sum(rle.index$lengths[1:2])
        
        n.left <- min(left.count, floor(n.extra/2))
        index[(i.left - n.left):(i.left - 1)] <- TRUE
        n.extra <- n.extra - n.left

        n.right <- min(right.count, n.extra)
        index[(i.right + 1):(i.right + n.right)] <- TRUE
        n.extra <- n.extra - n.right
      # If there are two entries, they are TRUE FALSE, or FALSE TRUE
      } else if ( length(rle.index$value) == 2 ) {
        if ( rle.index$value[1] ) {
          right.count <- rle.index$lengths[2]
          i.right <- rle.index$lengths[1:2]
          n.right <- min(right.count, n.extra)
          index[(i.right + 1):(i.right + n.right)] <- TRUE
          n.extra <- n.extra - n.right
        } else {
          left.count <- rle.index$lengths[1]
          i.left <- left.count + 1
          n.left <- min(left.count, n.extra)
          index[(i.left - n.left):(i.left - 1)] <- TRUE
          n.extra <- n.extra - n.left
        }
      } 
    
      # Checking if there are any more indexes unaccounted for
      if ( n.extra > 0 ) {
        msg <- 'It was not possible to round desired range to a power of 2.'
        warning(msg, call. = FALSE)
      }
    
    }

    object@processed <- object@processed[index, ]

    # Adjusting length of convolution vectors
    if ( length(object@product) > 0 ) {
      n <- 2*sum(index) - 1
      product <- object@product
      object@product <- approx(1:length(product), product, n = n)$y
      object@product <- object@product/sum(object@product)
      object@convolution <- Re(fftshift(fft(object@product)))
    }

    object
  })

#========================================================================>
# Defining a plot function
#========================================================================>

#------------------------------------------------------------------------
#' Plot NMRData1D object
#'
#' Generates an interactive plot object using the plotly package.
#'
#' @param x NMRData1D object.
#' @param components One of either 'r', 'i', or 'r/i' to include real,
#'                   imaginary or both components. If both components
#'                   are selected, they are displayed in separate subplots.
#' @param nrows Max number of rows to display subplots.
#'
#' @export
plot.NMRData1D <- function(x, components = 'r', nrows = 2) {
  d <- processed(x)

  plots <- list()
  legend.opts <- list(orientation = 'h', xanchor = "center", x = 0.5)

  if ( grepl('r', components) ) {
    p <- plot_ly(x = d$direct.shift, y = Re(d$intensity), 
                 color = I('black'), name = 'Real') %>%
           add_lines()  %>%
           layout(legend = legend.opts)

    plots$r <- p
  }

  if ( grepl('i', components) ) {
    p <- plot_ly(x = d$direct.shift, y = Im(d$intensity), 
                 color = I('grey'), name = 'Imaginary') %>%
           add_lines()  %>%
           layout(legend = legend.opts)

    plots$i <- p
  }

  if ( length(plots) == 0 ) NULL
  else subplot(plots, shareX = TRUE, shareY = TRUE, 
               nrows = min(length(plots), nrows))
}

setMethod("plot", "NMRData1D", plot.NMRData1D)
