# Functions for import Bruker parameter and data files

#========================================================================>
# fid and ser

#------------------------------------------------------------------------
#' Read 1D Bruker fid file
#'
#' Reads a raw bruker 1D fid file based on specified parameters.
#'
#' @param path Character string that points to either an fid file or
#'             the broader scan directory.
#' @param acqus.list A list of acqus parameters that contains 'td', 'sw',
#'                   'sfo1', and 'bytorda' entries. This list can be
#'                   generated using read_acqus_1d() or through other means.
#'                   These parameters can also be nested within a list
#'                   item called 'acqus' if multiple dimensions are read
#'                   at once.
#'
#' @return A data.frame made of two columns -- "direct.time" and "signal",
#'         corresponding to sampling time and the complex signal,
#'         respectively.
#'
#' @export
read_signal_1d <- function(path, acqus.list) {

  # Checking for required acqus entries
  required.acqus <- c('td', 'bytorda', 'sw', 'sfo1')
  acqus.list <- validate_param(acqus.list, required.acqus, 'acqus')

  # Extracting parameters
  bytorda <- acqus.list$bytorda
  td <- acqus.list$td
  sw <- acqus.list$sw
  sfo1 <- acqus.list$sfo1

  # Setting file path if the path is a directory
  if (dir.exists(path)) {
    path <- file.path(path, 'fid')
  }

  # Reading fid file
  if (bytorda == 1) {
    endian <- 'big'
  } else {
    endian <- 'little'
  }

  signal.data <- safe_read(path, 'bin', size = 4, what = 'integer',
                           n = td, endian = endian)

	if (length(signal.data) < td){
    msg <- sprintf('Error reading fid file, file size does not match data')
		stop(msg)
	}

  signal <- complex(real = signal.data[seq(1, td, by=2)], 
                    imaginary = signal.data[seq(2, td, by=2)]) 

  # Formatting the x-axis
  signal.time <- seq(0, td/(2*sw*sfo1), length.out = td/2)

  # Combining output
  data.frame(direct.time = signal.time, signal = signal)
}

#------------------------------------------------------------------------
#' Read 2D Bruker ser file
#'
#' Reads a raw bruker 2D ser file based on specified parameters.
#'
#' @param path Character string that points to either a ser file or
#'             the broader scan directory.
#' @param acqus.list A list of lists containing acqus parameters with 
#'                   'td', 'sw', 'sfo1', and 'bytorda' entries. 
#'                   This list can be generated using read_acqus() or 
#'                   through other means. Unless the sublists are named
#'                   acqus and acqu2s, the first and second sublits are
#'                   assumed to correspond to direct and indirect dimensions.
#'
#' @return A data.frame made of three columns -- "direct.time",
#'         "indirect.timea," and "signal", corresponding to signal sampling
#'         time, the delay time, and the complex signal, respectively.
#'
#' @export
read_signal_2d <- function(path, acqus.list) {

  # The acqus.list must contain at least two sublists 
  if (length(acqus.list) < 2) {
    msg <- 'acqus.list must contain two sublists'
    stop(msg)
  }

  # Picking off acqus and acqu2s
  if (all(c('acqus', 'acqu2s') %in% names(acqus.list)) ) {
    direct.acqus <- acqus.list$acqus
    indirect.acqus <- acqus.list$acqu2s
  } else {
    direct.acqus <- acqus.list[[1]]
    indirect.acqus <- acqus.list[[2]]
  }

  # Checking for required acqus entries
  direct.required <- c('td', 'bytorda', 'sw', 'sfo1')
  direct.acqus <- validate_param(direct.acqus, direct.required)

  # Checking for required acqu2s entries
  indirect.required <- c('td', 'sw', 'sfo1')
  indirect.acqus <- validate_param(indirect.acqus, indirect.required)

  # Extracting parameters
  bytorda <- direct.acqus$bytorda
  td <- c(direct.acqus$td, indirect.acqus$td)
  sw <- c(direct.acqus$sw, indirect.acqus$sw)
  sfo1 <- c(direct.acqus$sfo1, indirect.acqus$sfo1)

  # Setting file path if the path is a directory
  if (dir.exists(path)) {
    path <- file.path(path, 'ser')
  }

  # Reading ser file
  if (bytorda == 1) {
    endian <- 'big'
  } else {
    endian <- 'little'
  }
  
  signal.data <- safe_read(path, 'bin', size = 4, what = 'integer',
                           n = td[1]*td[2], endian = endian) 

	if (length(signal.data) < td[1]*td[2]){
    msg <- sprintf('Error reading ser file, file size does not match data')
		stop(msg)
	}

  signal <- complex(real = signal.data[seq(1, td[1]*td[2], by=2)], 
                    imaginary = signal.data[seq(2, td[1]*td[2], by=2)]) 

  # Formatting the x- and y- axes
  direct.time <- rep(seq(0, td[1]/(2*sw[1]*sfo1[1]), length.out = td[1]/2),
                     td[2])
  indirect.time <- rep(seq(0, td[2]/(2*sw[2]*sfo1[2]), length.out = td[2]),
                       each = td[1]/2)

  # Combining output
  data.frame(direct.time = direct.time, indirect.time = indirect.time,
             signal = signal)
}

#========================================================================>
# Processed data


#------------------------------------------------------------------------
#' Read 2D Bruker rr/ri/ir/ii files 
#'
#' Reads processed bruker 2D files based on specified parameters.
#' As the processed data can vary considerably based on the quadrature 
#' method used, the resulting data.frame will also vary quite a bit. 
#' If all 2rr/2ri/2ir/2ii files are present, they are stored as 2 complex 
#' valued columns named intensity (2rr/2ri) and intensity2 (2ir/2ii). 
#' If only two files are present, they are stored as a complex valued 
#' column intensity, regardless of which files are actually present, 
#' e.g. 2rr/2ii or 2rr/2ir. 
#'
#' @param path Character string that points to a Bruker scan directory.
#' @param procs.list A list of lists containing procs parameters with 
#'                   'sw.p', 'si', 'sf', 'reverse', and 'offset'  entries. 
#'                   This list can be generated using read_procs() or 
#'                   through other means. Unless the sublists are named
#'                   procs and proc2s, the first and second sublits are
#'                   assumed to correspond to direct and indirect dimensions.
#' @param number The processing file number. Defaults to the smallest number
#'               in the pdata directory.
#' @param hypercomplex TRUE to output full quadrature components
#'                     (2rr, 2ri, 2ir, 2ii) if present, FALSE to omit imaginary 
#'                     components in the indirect dimension (2ir, 2ii). 
#'
#' @return A data.frame made of three or four columns -- "direct.shift"
#'         containing the direct dimension chemical shift, "indirect.shift" 
#'         containing the indirect dimension chemical shift, "intensity",
#'         and possibly "intensity2". See description for more details.
#'
#' @export
read_processed_2d <- function(path, procs.list, number = NA, 
                              hypercomplex = TRUE) {

  # The procs.list must contain at least two sublists 
  if (length(procs.list) < 2) {
    msg <- 'procs.list must contain two sublists'
    stop(msg)
  }

  # Picking off procs and proc2s
  if (all(c('procs', 'proc2s') %in% names(procs.list)) ) {
    direct.procs <- procs.list$procs
    indirect.procs <- procs.list$proc2s
  } else {
    direct.procs <- procs.list[[1]]
    indirect.procs <- procs.list[[2]]
  }

  # Checking for required procs entries
  direct.required <- c('sw.p', 'si', 'sf', 'reverse', 'offset')
  direct.procs <- validate_param(direct.procs, direct.required)

  # Checking for required proc2s entries
  indirect.required <- c('sw.p', 'si', 'sf', 'reverse', 'offset')
  indirect.procs <- validate_param(indirect.procs, indirect.required)

  # Extracting parameters
  sw.p <- c(direct.procs$sw.p, indirect.procs$sw.p)
  si <- c(direct.procs$si, indirect.procs$si)
  sf <- c(direct.procs$sf, indirect.procs$sf)
  rv <- c(direct.procs$reverse, indirect.procs$reverse)
  ofs <- c(direct.procs$offset, indirect.procs$offset)

  n <- si[1]*si[2]

  # Extending and validating path
  path <- validate_pdata_path(path, number)

  # Checking which files actually exist
  all.paths <- file.path(path, c('2rr', '2ri', '2ir', '2ii'))
  existing.paths <- file.exists(all.paths)

  # If at least 2 files present, pick off the first two
  if ( sum(existing.paths) >= 2 ) {

    path.real <- all.paths[existing.paths][1]
    path.imag <- all.paths[existing.paths][2]

    # Reading binary files
    real.data <- safe_read(path.real, 'bin', size = 4, 
                           what = 'integer', n = n)
    imag.data <- safe_read(path.imag, 'bin', size = 4, 
                           what = 'integer', n = n)

  	if ( (length(real.data) < n) | (length(imag.data) < n) ) {
      msg <- 'Error reading processed files, file size does not match data.'
		  stop(msg)
	  }

    real.intensity  <- complex(real = real.data, imaginary = imag.data)

  } else {
    msg <- paste('Error reading processed files, at least two of',
                 '2rr/2ri/2ir/2ii must exist.')
    stop(msg)
  }

  if ( sum(existing.paths) == 4 ) {

    path.real2 <- all.paths[existing.paths][3]
    path.imag2 <- all.paths[existing.paths][4]

    real2.data <- safe_read(path.real2, 'bin', size = 4, 
                            what = 'integer', n = n)
    imag2.data <- safe_read(path.imag2, 'bin', size = 4, 
                            what = 'integer', n = n)

    if ( (length(real2.data) < n) | (length(imag2.data) < n) ) {
      msg <- 'Error reading processed files, file size does not match data'
      stop(msg)
    }

    imag.intensity  <- complex(real = real2.data, imaginary = imag2.data)

  } else {

    imag.intensity <- NULL

  }

  # Formatting direct.shift 
  direct.shift <- seq(ofs[1], ofs[1] - sw.p[1]/sf[1], length.out = si[1])
  
  if (rv[1] == 'yes') {
    direct.shift <- rev(direct.shift)
  }

  # Formatting indirect.shift 
  indirect.shift <- seq(ofs[2], ofs[2] - sw.p[2]/sf[2], length.out = si[2])
  
  if (rv[2] == 'yes') {
    indirect.shift <- rev(indirect.shift)
  }

  # Combining output
  d <- data.frame(direct.shift = rep(direct.shift, si[2]), 
                  indirect.shift = rep(indirect.shift, each = si[1]),
                  intensity = real.intensity)

  if ( hypercomplex & (! is.null(imag.intensity)) ) {
    d$intensity2 <- imag.intensity
  }

  d
}




