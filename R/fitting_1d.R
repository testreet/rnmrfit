# Proto classes handling the actual optimization

#------------------------------------------------------------------------
# Initializing fit proto class
#' Fit proto class
#'
#' Fit proto class.
#'
#' @export
fit.env <- proto()

#------------------------------------------------------------------------
#' Initialization function for fit proto class
#'
#' Generates a lineshape optimization function based on specified
#' parameters by combining multiple individual peaks with phase,
#' baseline, and convolution terms. 
#'
#' @param . Proto object.
#' @param x Normalized x data for peak fit.
#' @param y Normalized y data for peak fit.
#' @param nc Number of curve parameters.
#' @param nb1 Number of baseline points for real baseline.
#' @param nb2 Number of baseline points for baseline difference.
#' @param knots The interior baseline knots.
#' @param np Number of phase parameters.
#' @param lineshape One of either "lorenz", "gauss", "pvoigt", or "voigt".
#' @param convolution Frequency domain convolution window.
#' @param components One of either "r/i", "r", or "i".
#'
#' @export
fit.env$init <- function(., x, y, nc, nb1, nb2, knots, np, 
                 lineshape, convolution, components) {

  # Setting shared data objects
  .$x <- x
  .$y <- y
  n <- length(x)
  .$p <- rep(0, nc + nb1 + nb2 + np)
  .$np <- np
  .$fit <- as.complex(rep(0, n))
  .$grad <- matrix(complex(0, 0), nrow = n, ncol = nc + nb1 + nb2 + np)

  # Defining the fit function as a collection of separate components,
  # each accesing a different part of the input (parameter) vector
  f.list <- list()
  
  j <- ifelse( grepl('voigt', lineshape), 4, 3) 
  
  gen_f <- switch( lineshape, lorenz = .$fit_lorenz, gauss = .$fit_gauss,
                              pvoigt = .$fit_pvoigt, voigt = .$fit_voigt )

  # Incrementally adding lineshape fit terms
  for ( i in seq(1, nc, j) ) {
    index <- i:(i + j - 1)
    f.list <- c(f.list, gen_f(index))
  }

  # Adding lineshape convolution if necessary
  if (! is.null(convolution) ) {
    f.list <- c(f.list, .$fit_convolve(convolution, nc))
  }

  # Adding baseline and baseline difference if needed
  if ( nb1 > 0 ) {
    index <- (nc + 1):(nc + nb1)
    f.list <- c(f.list, .$fit_baseline(index, knots))
  }

  if ( nb2 > 0 ) {
    index <- (nc + nb1 + 1):(nc + nb1 + nb2)
    f.list <- c(f.list, .$fit_baseline_difference(index, knots))
  }

  # Storing f.list and defining fit function
  .$f.list <- f.list

}

#------------------------------------------------------------------------
#' Global fit function for proto class
#'
#' Calculates overall fit by combining all element that the class was
#' initialized with. Returns the least squares value and gradient.
#'
#' @param . Proto object.
#' @param p Vector of parameters
#'
#' @export
fit.env$eval <- function(., p) {
  .$p <- p
  .$fit[] <- 0
  .$grad[,] <- 0
  
  for ( f in .$f.list ) {
    f(.)
  }

  # Applying phase correction if necessary
  if ( .$np > 0 ) {
    angle = .$p[length(p)]

    yp <- complex(re =  Re(.$y)*cos(angle) + Im(.$y)*sin(angle),
                  im = -Re(.$y)*sin(angle) + Im(.$y)*cos(angle)) 
  } else {
    yp <- .$y
  }

  # Calculating objective function
  eps <- yp - .$fit
  obj <- sum(Re(eps)*Re(eps)) + sum(Im(eps)*Im(eps))

  grad <- -2*(apply(Re(.$grad), 2, function(x) {sum(Re(eps) * x)}) + 
              apply(Im(.$grad), 2, function(x) {sum(Im(eps) * x)}))

  # Adding phase derivative
  if ( .$np > 0 ) {
    i <- length(p)
    grad[i] <- 2*(sum(Re(eps) * (-Re(y)*sin(angle) + Im(y)*cos(angle))) + 
                  sum(Im(eps) * (-Re(y)*cos(angle) - Im(y)*sin(angle))))
  }

  list(objective = obj, gradient = grad)
} 

#------------------------------------------------------------------------
#' Generator function for lorenz fit
#'
#' Assigns a set of fit parameters to a single lorenz peak. The resulting
#' closure calculates peak fit based on parameters stored in the proto
#' class. The results are added to internal .$fit and .$grad arrays.
#'
#' @param . Proto object.
#' @param index Index referencing which values of .$p should be used
#'
#' @export
fit.env$fit_lorenz <- function(., index) {

  force(index)

  p.index <- index[1]
  hl.index <- index[2]
  wl.index <- index[3]

  f <- function(.) {
  
    # Extracting parameters
    p <- .$p[p.index]
    hl <- .$p[hl.index]
    wl <- .$p[wl.index]

    # Common vector terms
    z <- (.$x - p)/wl
    z2 <- z*z
    z2p1 <- z2 + 1

    f <- (1 + complex(i = z))/z2p1
    dfdz <- hl*(-2*z + complex(i = 1 - z2))/(z2p1*z2p1)

    # Derivatives of z for chain rule
    dzdp <- complex(r = -1/wl)
    dzdwl <- complex(r = -z/wl)

    #----------
    # Overall fit
    .$fit <- .$fit + hl*f

    #----------
    # Gradients

    # Position
    .$grad[, index[1]] <- dfdz*dzdp

    # Height
    .$grad[, index[2]] <- f 

    # Lorenz width
    .$grad[, index[3]] <- dfdz*dzdwl

  }
  
  f
}

#------------------------------------------------------------------------
#' Generator function for voigt fit
#'
#' Assigns a set of fit parameters to a single voigt peak. The resulting
#' closure calculates peak fit based on parameters stored in the proto
#' class. The results are added to internal .$fit and .$grad arrays.
#'
#' @param . Proto object.
#' @param index Index referencing which values of .$p should be used
#'
#' @export
fit.env$fit_voigt <- function(., index) {

  force(index)

  p.index <- index[1]
  h.index <- index[2]
  wl.index <- index[3]
  wg.index <- index[4]

  f <- function(.) {

    # Extracting parameters
    p <- .$p[p.index]
    h <- .$p[h.index]
    wl <- .$p[wl.index]
    wg <- .$p[wg.index]

    # Common vector terms
    z <- (.$x - p + complex(i = wl))/(sqrt(2)*wg)
    
    f <- Faddeeva_w(z)
    dfdz <- -2*z*f + complex(i = 2/sqrt(pi))

    # Normalizing factors
    z0 <- complex(i = wl)/(sqrt(2)*wg)
    f0 <- Faddeeva_w(z0) 
    f02 <- f0 * f0
    df0dz0 <- -2*z0*f0 + complex(i = 2/sqrt(pi))

    fnrm <- f/f0

    # Derivatives of z for chain rule
    dzdp <- complex(r = -1/(sqrt(2)*wg))
    dzdwl <- complex(i = 1/(sqrt(2)*wg))
    dzdwg <- -z/wg
    dz0dwg <- -complex(i = wl)/(sqrt(2)*wg*wg)

    #----------
    # Overall fit
    .$fit <- .$fit + h*fnrm

    #----------
    # Gradients

    # Position
    .$grad[, index[1]] <- h*dfdz*dzdp/f0

    # Height
    .$grad[, index[2]] <- fnrm

    # Lorenz width
    .$grad[, index[3]] <- h*(dfdz*f0 - df0dz0*f)/f02*dzdwl

    # Gauss width
    .$grad[, index[4]] <- h*(dfdz*dzdwg*f0 - df0dz0*dz0dwg*f)/f02

  }

  f
}

#------------------------------------------------------------------------
#' Generator function for lineshape convolution
#'
#' Convolutes fit and gradient terms with specified convolution window.
#' The convolution window should be 2*n - 1 in length where n is the
#' length of the data vector.
#'
#' @param . Proto object.
#' @param index Index referencing which values of .$p should be used
#' @param nc The number of curve parameters
#'
#' @export
fit.env$fit_convolve <- function(., convolution, nc) {

  force(convolution)
  force(nc)

  n <- length(.$fit)
  index <- (1 + n):(2*n)

  f <- function(.) {

    .$fit <- convolve(.$fit, convolution, type = 'open')[index]

    f_apply <- function(x) { convolve(x, convolution, type = 'open')[index] }
    .$grad[, 1:nc] <- apply(.$grad[, 1:nc], 2, f_apply)

  }

  f
}

#------------------------------------------------------------------------
#' Generator function for baseline fit
#'
#' Assigns a set of fit parameters to the baseline and precalculates
#' basis. Calculates b-spline baseline based on current parameters stored 
#' in the proto class. The results are added to internal .$fit and 
#' .$grad arrays.
#'
#' @param . Proto object.
#' @param index Index referencing which values of .$p should be used
#' @param knots Internal knots used to generate basis
#'
#' @export
fit.env$fit_baseline <- function(., index, knots) {

  force(index)
  force(knots)

  degree <- length(index) - length(knots) - 1

  # B-spline basis 
  basis <- bs(x, degree = degree, knots = knots, intercept = TRUE)

  # The gradient term is a complex version of the basis
  gradient <- complex(re = basis, im = basis)

  f <- function(.) {

    # Extracting parameters
    p <- matrix(.$p[index], ncol = 1)

    #----------
    # Overall fit
    fit <- basis %*% p
    .$fit <- .$fit + complex(re = fit, im = fit)

    #----------
    # Gradients
    .$grad[, index] <- gradient

  }

  f
}

#------------------------------------------------------------------------
#' Generator function for baseline difference fit
#'
#' Assigns a set of fit parameters to the baseline difference and 
#' precalculates basis. Calculates b-spline baseline based on current 
#' parameters stored in the proto class. The results are added to internal 
#' .$fit and .$grad arrays.
#'
#' @param . Proto object.
#' @param index Index referencing which values of .$p should be used
#' @param knots Internal knots used to generate basis
#'
#' @export
fit.env$fit_baseline_difference <- function(., index, knots) {

  force(index)
  force(knots)

  degree <- length(index) - length(knots) - 1

  # B-spline basis 
  basis <- bs(x, degree = degree, knots = knots, intercept = TRUE)

  # The gradient term is a complex version of the basis
  gradient <- complex(im = basis)

  f <- function(.) {

    # Extracting parameters
    p <- matrix(.$p[index], ncol = 1)

    #----------
    # Overall fit
    .$fit <- .$fit + complex(im = basis %*% p)

    #----------
    # Gradient (just the basis function)
    .$grad[, index] <- gradient

  }

  f
}
