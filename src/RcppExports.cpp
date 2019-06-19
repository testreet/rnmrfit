// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// fit_lineshape_1d
double fit_lineshape_1d(const Rcpp::NumericVector x, const Rcpp::ComplexVector y, Rcpp::NumericVector par);
RcppExport SEXP _rnmrfit_fit_lineshape_1d(SEXP xSEXP, SEXP ySEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::ComplexVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_lineshape_1d(x, y, par));
    return rcpp_result_gen;
END_RCPP
}
// lineshape_1d
std::vector< std::complex<double> > lineshape_1d(const NumericVector x, const NumericMatrix par);
RcppExport SEXP _rnmrfit_lineshape_1d(SEXP xSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(lineshape_1d(x, par));
    return rcpp_result_gen;
END_RCPP
}
// test_list
void test_list(const List x);
RcppExport SEXP _rnmrfit_test_list(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type x(xSEXP);
    test_list(x);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rnmrfit_fit_lineshape_1d", (DL_FUNC) &_rnmrfit_fit_lineshape_1d, 3},
    {"_rnmrfit_lineshape_1d", (DL_FUNC) &_rnmrfit_lineshape_1d, 2},
    {"_rnmrfit_test_list", (DL_FUNC) &_rnmrfit_test_list, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rnmrfit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
