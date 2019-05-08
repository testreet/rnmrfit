# Definition of global options

#' Set global options
#' 
#' A number of common parameters are used across different functions. Such
#' global options are assumed to be valid for all (or at least most) analysis
#' carried out during a single session.
#' 
#' @param sf Sweep frequncy (MHz). Used to convert peak width and coupling
#'           information from Hz into ppm.
#' 
#' @export
nmrsession_1d <- set_opt(
  "sf" = list(.value = NULL,
              .length = 1,
              .class = "numeric"),
  "baseline" = list(.value = list(order = 3, n.knots = 0),
                    .length = 2,
                    .class = "list",
                    .validate = function(x) {
                      all(names(x) %in% c('degree', 'n.knots'))
                    }),
  "phase" = list(.value = list(order = 0),
                 .length = 1,
                 .class = "list",
                 .validate = function(x) {names(x) == ''})
  )
