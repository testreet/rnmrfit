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






