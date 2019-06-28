# Definition of global options

#------------------------------------------------------------------------------
#' Set global options for 1D peak processing and fitting
#' 
#' A number of common parameters are used across different functions. Such
#' global options are assumed to be valid for all (or at least most) analysis
#' carried out during a single session.
#' 
#' @export
nmrsession_1d <- set_opt(
  "sf" = list(
    .value = NULL,
    .length = 1,
    .class = "numeric"
  ),
  
  "baseline" = list(
    .value = list(order = 3, n.knots = 0),
    .length = 2,
    .class = "list",
    .validate = function (x) { all(names(x) %in% c('degree', 'n.knots')) }
  ),
  
  "phase" = list(
    .value = list(order = 0),
    .length = 1,
    .class = "list",
    .validate = function (x) { names(x) == 'order'} 
  ),
  
  "exclusion" = list(
    .value = list(level = 'peak', notification = 'warning'),
    .length = 2,
    .class = "list",
    .validate = function (x) {
      all(names(x) %in% c('level', 'notification')) &&
      (x[['level']] %in% c('peak', 'resonance', 'species')) &&
      (x[['notification']] %in% c('none', 'message', 'warning', 'stop'))
  }),

  "fit" = list(
    .value = list(
      opts = list(),
      init = function (object, ...) {

        args <- list(...)[c('exclusion.level', 'exclusion.notification')]
        args$nmrdata <- object@nmrdata
        args <- c(list(object = object), args)

        object <- do.call(initialize_heights, args)
        object <- set_conservative_bounds(object, nmrdata = object@nmrdata)
        object
      }
    ),
    .length = 2,
    .class = "list",
    .validate = function (x) {
      all(names(x) %in% c('opts', 'init')) &&
      (class(x[['init']]) == 'function')
  })
)

#------------------------------------------------------------------------------
#' Set global options for 2D peak processing and fitting
#' 
#' A number of common parameters are used across different functions. Such
#' global options are assumed to be valid for all (or at least most) analysis
#' carried out during a single session.
#' 
#' @export
nmrsession_2d <- set_opt(
  "sf" = list(
    .value = list(
      direct = NULL,
      indirect = NULL
    ),
    .length = 2,
    .class = "list",
    .validate = function (x) {
      all(names(x) %in% c('direct', 'indirect'))
  })
)


# Some potential plot options to consider in the future:
# "plot" = list(.value = 
#    list('legend.entries' = c('data', 'fit', 'baseline', 'residual'),
#         'nrows'
#          'legend.position'
#          'reverse')
