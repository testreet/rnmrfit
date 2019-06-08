# Definition of a class structure for 1D NMR data.



#------------------------------------------------------------------------------
#' A class combining NMR data and scan parameters
#' 
#' An extension of the generic NMRData class to provide 1D-specific methods
#' such as constructors and plotting.
#' 
#' @rdname NMRData1D
#' @export
NMRData1D <- setClass("NMRData1D", contains = "NMRData")

# Validity testing consists of simply checking the processed data.frame columns
validNMRData1D <- function(object) {

  valid <- TRUE
  err <- c()

  processed <- object@processed
  acqus <- object@acqus
  procs <- object@procs

  # Checking processed columns
  valid.columns <- c('direct.shift', 'intensity')
  msg <- sprintf('Processed data may only have the following columns: %s',
                 paste(valid.columns, collapse = ', '))
  if (! identical(valid.columns, colnames(processed)) ) {
    valid <- FALSE
    err <- c(err, msg)
  }

  if ( valid ) TRUE
  else msg
}

setValidity("NMRData1D", validNMRData1D)



#==============================================================================>
#  Constructors and data loading
#==============================================================================>



#------------------------------------------------------------------------------
#' Constructors for generating an NMRData1D object
#' 
#' \code{nmrdata_1d()} can be used as a generic constructor method that handles
#' different types of input (currently limited to Bruker pdata directories and
#' JCAMP=DX files). If the input path points to a file, it is assumed to be a
#' JCAMP-DX file, and if the input path points to a directory, it is assumed to
#' be Bruker scan directory. The specific constructor functions can also be
#' called using \code{nmrdata_1d_from_pdata()} or
#' \code{nmrdata_1d_from_jcamp()}.
#' 
#' @param path Path to a Bruker scan directory or JCAMP-DX file.
#' @param procs.number Specifies pdata directory to load. Defaults to the lowest
#'                     available number. Ignored if loading from a JCAMP-DX
#'                     file.
#' @param blocks.number Specifies block number to use when loading from a JCAMP-
#'                      DX file. Defaults to the first block encountered.
#'                      Ignored if loading from a pdata directory.
#' @param ntuples.number Specifies ntuple entry number to use when loading from
#'                       a JCAMP-DX file. Defaults to the first ntuple entry
#'                       encountered. Ignored if loading from a pdata directory.
#' 
#' @return An NMRData1D object.
#' 
#' @export
nmrdata_1d <- function(path, procs.number = NA, 
                       blocks.number = 1, ntuples.number = 1) {

  # If path is a valid directory, treat as Bruker scan directory
  if ( file.exists(path) && dir.exists(path) ) {
    nmrdata_1d_from_pdata(path, procs.number) 
  }
  # If it's not a directory, it might be a file
  else if ( file.exists(path) )  {
    nmrdata_1d_from_jcamp(path, blocks.number, ntuples.number)
  }
  # Otherwise, error out
  else {
    err <- sprintf('Path "%s" does not point to a file or directory', path)
    stop(err)
  }

}



#------------------------------------------------------------------------------
#' @rdname nmrdata_1d
#' @export
nmrdata_1d_from_pdata <- function(path, procs.number = NA) {

  # First, loading procs parameters
  procs <- read_procs_1d(path, procs.number)

  # Using the procs file to load the processed data
  processed <- read_processed_1d(path, procs, procs.number)

  # Finally, loading the general acquisition parameters
  acqus <- read_acqus_1d(path)

  # Returning class object
  new("NMRData1D", processed = processed, parameters = list(),
                   procs = procs, acqus = acqus)

}



#------------------------------------------------------------------------------
#' @rdname nmrdata_1d
#' @export
nmrdata_1d_from_jcamp <- function(path, blocks.number = 1, ntuples.number = 1) {

  jcamp <- read_jcamp(path, process.tags = TRUE, process.entries = TRUE)

  # Double check that specified block exists
  err <- sprintf("Specified block number not found in %s", path)
  if ( blocks.number > length(jcamp$blocks) ) stop(err)

  jcamp$blocks <- jcamp$blocks[blocks.number]

  # Check to make sure that ntuples exist
  err <- 'Import from JCAMP file currently limited to NTUPLES entries'
  if (! 'ntuples' %in% names(jcamp$blocks[[1]]) ) stop(err)

  # Double check that specified ntuple exists
  err <- sprintf('Specified NTUPLES number not found in block %i of %s', 
                 ntuples.number, path)
  if ( ntuples.number > length(jcamp$blocks[[1]]) ) stop(err)

  jcamp$blocks[[1]]$ntuples <- jcamp$blocks[[1]]$ntuples[ntuples.number]

  # Flattening file
  jcamp.flat <- flatten_jcamp(jcamp)

  # Extracting processed data from ntuples
  descriptors <- jcamp.flat$ntuples$descriptors
  pages <- jcamp.flat$ntuples$pages

  # Checking variable names
  variables <- tolower(descriptors$var.name)
  real.index <- which(str_detect(variables, 'spectrum.*real'))
  imag.index <- which(str_detect(variables, 'spectrum.*imag'))

  # If both real and imaginary data isn't there, abort
  err <- 'Import from JCAMP file currently limited to real/imaginary spectra'
  if ( length(c(real.index, imag.index)) < 2 ) stop(err)

  # Double check that the first entry is frequency
  err <- 'Import from JCAMP file currently limited to frequency abscissa'
  if (! str_detect(variables[1], 'freq') ) stop(err)

  # Proceed to extract data
  real.data <- pages[[real.index - 1]]
  imag.data <- pages[[imag.index - 1]]

  # Checking that frequency is the same
  real.frequency <- real.data[, 1]
  imag.frequency <- imag.data[, 1]
  err <- 'Mismatch in real and imaginary frequency, likely parsing error'
  if ( any(abs(real.frequency - imag.frequency) > 1e-6) ) stop(err)

  # Starting with raw values
  frequency <- real.frequency
  real.data <- real.data[, 2]
  imag.data <- imag.data[, 2]

  # Scaling if factors are non zero
  scale <- descriptors$factor[1]
  if (scale > 1e-6) frequency <- frequency*scale

  scale <- descriptors$factor[real.index]
  if (scale > 1e-6) real.data <- real.data*scale

  scale <- descriptors$factor[imag.index]
  if (scale > 1e-6) imag.data <- imag.data*scale

  # Offsetting if max-min difference is non zero
  d.max <- descriptors$max[1]
  d.min <- descriptors$min[1]
  if ( (d.max - d.min) > 1e-6 ) {
    frequency <- frequency - max(frequency) + d.max
  }

  d.max <- descriptors$max[real.index]
  d.min <- descriptors$min[real.index]
  if ( (d.max - d.min) > 1e-6 ) {
    real.data <- real.data - max(real.data) + d.max
  }

  d.max <- descriptors$max[imag.index]
  d.min <- descriptors$min[imag.index]
  if ( (d.max - d.min) > 1e-6 ) {
    imag.data <- imag.data - max(imag.data) + d.max
  }

  # Doing one final check on the direct shift to check on offset
  direct.shift <- frequency/jcamp.flat$sf
  
  delta <- jcamp.flat$offset - max(direct.shift)
  if ( abs(delta) > 1e-6 ) direct.shift <- direct.shift + delta

  # Finally, combine the data
  intensity <- cmplx1(r = real.data, i = imag.data)
  processed <- tibble(direct.shift = direct.shift, intensity = intensity)

  # Returning class object
  new("NMRData1D", processed = processed, parameters = jcamp.flat,
                   procs = list(), acqus = list()) 
}



#==============================================================================>
#  Formatting and printing
#==============================================================================>



#' @export
format.NMRData1D <- function(x, ...) {
  d <- processed(x)
  components <- paste(colnames(as_tibble(d$intensity)), collapse = ', ') 
  shift.range <- range(d$direct.shift)
  shift.range <- sprintf("%.3f ppm to %.3f ppm", shift.range[1], shift.range[2])
  sprintf('NMRData1D object (%s), %s\n', components, shift.range)
}

#' @export
print.NMRData1D <- function(x, ...) cat(format(x))

#' @export
setMethod("show", "NMRData1D", 
  function(object) cat(format(object))
  )

#' @export
summary.NMRData1D <- function(object, ...) summary(object@processed)



#' @export
is_vector_s3.NMRData1D <- function(x) FALSE

#' @export
type_sum.NMRData1D <- function(x) "NMRData1D"

#' @export
obj_sum.NMRData1D <- function(x) format(x)



#==============================================================================>
#  Non-inherited methods
#==============================================================================>



#------------------------------------------------------------------------------
#' Set convolution applied to processed data
#' 
#' Although it's possible to fit the processed data based on its final
#' appearance, it may sometimes be more useful to try and fit the underlying
#' lineshape with an explicit consideration of the convolution that it has
#' undergone. This convolution includes any apodization as well as the process
#' of zero-filling, which can be considered as the convolution of a
#' hypothetically longer experiment by a step function. By default, the applied
#' apodization is read from the procs parameters, but these can be overwritten
#' by either a list of parameters or a custom function.
#' 
#' @param object An NMRData1D object.
#' @param trim Number of point to truncate from the end of the fid's acquisition
#'             time (not including zero fill). It's common to circle shift
#'             Bruker's group delay artefact from the start of a collected fid
#'             to the last 50-100 points. While most apodization methods will
#'             passively get rid of or diminish the influence of these points,
#'             they may have a stronger impact with minimal or no apodization.
#' @param param A list of lists combining apodization names and parameters, (see
#'              apodize_signal) for a full list of apodizations and their
#'              parameters. E.g. list('exponential'=list('lb'=1)).
#' @param f A function f(n) capable of generating the convolution vector of
#'          length n (ranging from 0 to the acquisition time). If provided, f
#'          overrides the param list.
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

    si <- get_parameter(object, 'si', 'procs')
    td <- get_parameter(object, 'td', 'acqus')
    sw <- get_parameter(object, 'sw.h', 'acqus')

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

      lb <- get_parameter(object, 'lb', 'procs')
      gb <- get_parameter(object, 'gb', 'procs')
      ssb <- get_parameter(object, 'ssb', 'procs')
      tm1 <- get_parameter(object, 'tm1', 'procs')
      tm2 <- get_parameter(object, 'tm2', 'procs')

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



#==============================================================================>
#  Defining a plot function
#==============================================================================>

#------------------------------------------------------------------------------
#' Plot NMRData1D object
#' 
#' Generates an interactive plot object using the plotly package.
#' 
#' @param x NMRData1D object.
#' @param components One of either 'r', 'i', or 'r/i' to include real, imaginary
#'                   or both components. If both components are selected, they
#'                   are displayed in separate subplots.
#' @param nrows Max number of rows to display subplots.
#' @param reverse TRUE to order x-axis with large values on the left and small
#'                values on the right like typical NMR plots.
#' 
#' @export
plot.NMRData1D <- function(x, components = 'r', nrows = 2, reverse = TRUE) {
  d <- processed(x)

  plots <- list()
  legend.opts <- list(orientation = 'h', xanchor = "center", x = 0.5)

  if ( grepl('r', components) ) {
    p <- plot_ly(x = d$direct.shift, y = Re(d$intensity), 
                 color = I('black'), name = 'Real') %>%
           add_lines()  %>%
           layout(legend = legend.opts)

    if ( reverse ) p <- p %>% layout(xaxis = list(autorange = "reversed"))

    plots$r <- p
  }

  if ( grepl('i', components) ) {
    p <- plot_ly(x = d$direct.shift, y = Im(d$intensity), 
                 color = I('grey'), name = 'Imaginary') %>%
           add_lines()  %>%
           layout(legend = legend.opts)

    if ( reverse ) p <- p %>% layout(xaxis = list(autorange = "reversed"))

    plots$i <- p
  }

  if ( length(plots) == 0 ) NULL
  else subplot(plots, shareX = TRUE, shareY = TRUE, 
               nrows = min(length(plots), nrows))
}

setMethod("plot", "NMRData1D", plot.NMRData1D)
