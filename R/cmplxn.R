# Definition of a new set of complex classes that expand beyond the i notation



#==============================================================================>
# 1D -- basically a copy of complex() to match the 2D version below
#==============================================================================>



new_cmplx1 <- function(r = double, i = double()) {
  vec_assert(r, ptype = double())
  vec_assert(i, ptype = double())
  
  new_rcrd(list(r = r, i = i), class = "vctrs_cmplx1")
}

cmplx1 <- function(r = 0, i = 0) {
  c(r, i) %<-% vec_cast_common(r, i, .to = double())
  c(r, i) %<-% vec_recycle_common(r, i)
  
  new_cmplx1(r, i)
}

format.vctrs_cmplx1 <- function(x, ...) {
  r <- field(x, "r")
  i <- field(x, "i")
  
  out <- paste0(r, "+", i, "i")
  out[is.na(r) | is.na(i)] <- NA
  
  out
}

vec_ptype_abbr.vctrs_cmplx1 <- function(x) "cmplx1"
vec_ptype_full.vctrs_cmple1 <- function(x) "complex1d"

#------------------------------------------------------------------------------
# New subsetting functions

`$.vctrs_cmplx1` <- function(x, name) field(x, name)
`$<-.vctrs_cmplx1` <- function(x, name, value) {
  field(x, name) <- value
  x
}

#------------------------------------------------------------------------------
# Traditional Re()/Im()

Re.vctrs_cmplx1 <- function(z) z$r
Im.vctrs_cmplx1 <- function(z) z$i



#==============================================================================>
# 2D 
#==============================================================================>



new_cmplx2 <- function(rr = double, ri = double(), ir = double(), ii = double) {
  vec_assert(rr, ptype = double())
  vec_assert(ri, ptype = double())
  vec_assert(ir, ptype = double())
  vec_assert(ii, ptype = double())
  
  new_rcrd(list(rr = rr, ri = ri, ir = ir, ii = ii), class = "vctrs_cmplx2")
}

cmplx2 <- function(rr = 0, ri = 0, ir = 0, ii = 0) {
  c(rr, ri, ir, ii) %<-% vec_cast_common(rr, ri, ir, ii, .to = double())
  c(rr, ri, ir, ii) %<-% vec_recycle_common(rr, ri, ir, ii)
  
  new_cmplx2(rr, ri, ir, ii)
}

format.vctrs_cmplx2 <- function(x, ...) {
  rr <- field(x, "rr")
  ri <- field(x, "ri")
  ir <- field(x, "ir")
  ii <- field(x, "ii")
  
  out <- paste0(rr, "+", ri, "j+", ir, "i+", ii, "ji")
  out[is.na(rr) | is.na(ri) | is.na(ir) | is.na(ii)] <- NA
  
  out
}

vec_ptype_abbr.vctrs_cmplx2 <- function(x) "cmplx2"
vec_ptype_full.vctrs_cmple2 <- function(x) "complex2d"

#------------------------------------------------------------------------------
# New subsetting functions

`$.vctrs_cmplx2` <- function(x, name) field(x, name)
`$<-.vctrs_cmplx2` <- function(x, name, value) {
  field(x, name) <- value
  x
}

#------------------------------------------------------------------------------
# Traditional Re()/Im()

Re.vctrs_cmplx2 <- function(z) z$rr
Im.vctrs_cmplx2 <- function(z) z$ii
