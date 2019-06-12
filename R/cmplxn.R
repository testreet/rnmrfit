# Definition of a new set of complex classes that expand beyond the i notation



#==============================================================================>
# 1D -- basically a copy of complex() to match the 2D version below
#==============================================================================>


#' @export
new_cmplx1 <- function(r = double, i = double()) {
  vec_assert(r, ptype = double())
  vec_assert(i, ptype = double())
  
  new_rcrd(list(r = r, i = i), class = "vctrs_cmplx1")
}

#' @export
cmplx1 <- function(r = 0, i = 0) {
  c(r, i) %<-% vec_cast_common(r, i, .to = double())
  c(r, i) %<-% vec_recycle_common(r, i)
  
  new_cmplx1(r, i)
}

#' @export
format.vctrs_cmplx1 <- function(x, ...) {
  r <- field(x, "r")
  i <- field(x, "i")
  
  signs <- c("1"="+", "0"="+", "-1"="")
  f <- function (x) signs[as.character(sign(x))]
  out <- paste0(r, f(i), i, "i")
  out[is.na(r) | is.na(i)] <- NA
  
  out
}

#' @export
vec_ptype_abbr.vctrs_cmplx1 <- function(x) "cmplx1"

#' @export
vec_ptype_full.vctrs_cmplx1 <- function(x) "complex1d"

#------------------------------------------------------------------------------
# New subsetting functions

#' @export
`$.vctrs_cmplx1` <- function(x, name) field(x, name)
`$<-.vctrs_cmplx1` <- function(x, name, value) {
  field(x, name) <- value
  x
}

#------------------------------------------------------------------------------
# Traditional Re()/Im()

#' @export
Re.vctrs_cmplx1 <- function(z) z$r

#' @export
Im.vctrs_cmplx1 <- function(z) z$i

#------------------------------------------------------------------------------
# Summary

#' @export
as_tibble.vctrs_cmplx1 <- function(x, ...) as_tibble(unclass(x))

#' @export
summary.vctrs_cmplx1 <- function(object, ..., 
                                 digits = max(3, getOption("digits") - 3)) {
  summary(as_tibble(object))

}



#==============================================================================>
# 2D 
#==============================================================================>


#' @export
new_cmplx2 <- function(rr = double, ri = double(), ir = double(), ii = double) {
  vec_assert(rr, ptype = double())
  vec_assert(ri, ptype = double())
  vec_assert(ir, ptype = double())
  vec_assert(ii, ptype = double())
  
  new_rcrd(list(rr = rr, ri = ri, ir = ir, ii = ii), class = "vctrs_cmplx2")
}

#' @export
cmplx2 <- function(rr = 0, ri = 0, ir = 0, ii = 0) {
  c(rr, ri, ir, ii) %<-% vec_cast_common(rr, ri, ir, ii, .to = double())
  c(rr, ri, ir, ii) %<-% vec_recycle_common(rr, ri, ir, ii)
  
  new_cmplx2(rr, ri, ir, ii)
}

#' @export
format.vctrs_cmplx2 <- function(x, ...) {
  rr <- field(x, "rr")
  ri <- field(x, "ri")
  ir <- field(x, "ir")
  ii <- field(x, "ii")
  
  signs <- c("1"="+", "0"="+", "-1"="")
  f <- function (x) signs[as.character(sign(x))]
  out <- paste0(rr, f(ri), ri, "j", f(ir), ir, "i", f(ii), ii, "ji")
  out[is.na(rr) | is.na(ri) | is.na(ir) | is.na(ii)] <- NA
  
  out
}

#' @export
vec_ptype_abbr.vctrs_cmplx2 <- function(x) "cmplx2"

#' @export
vec_ptype_full.vctrs_cmplx2 <- function(x) "complex2d"

#------------------------------------------------------------------------------
# New subsetting functions

#' @export
`$.vctrs_cmplx2` <- function(x, name) field(x, name)
`$<-.vctrs_cmplx2` <- function(x, name, value) {
  field(x, name) <- value
  x
}

#------------------------------------------------------------------------------
# Traditional Re()/Im()

#' @export
Re.vctrs_cmplx2 <- function(z) z$rr

#' @export
Im.vctrs_cmplx2 <- function(z) z$ii

#------------------------------------------------------------------------------
# Summary

#' @export
as_tibble.vctrs_cmplx2 <- as_tibble.vctrs_cmplx1

#' @export
summary.vctrs_cmplx2 <- summary.vctrs_cmplx1
