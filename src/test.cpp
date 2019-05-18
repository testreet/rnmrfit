#include <Rcpp.h>
#include "Faddeeva.h"

using namespace Rcpp;

// [[Rcpp::export]]
void test_list(const List x) {
  
  NumericVector y;
  int i, j;
  Rcout << "\n" << std::endl;
  for (i = 0; i < x.size(); i++) {
    y = x[i];
    for (j = 0; j < y.size(); j++) {
      Rcout << y[j] << ", " << std::endl;
    }
    Rcout << "\n" << std::endl;
  }

}
