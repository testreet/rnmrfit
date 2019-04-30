#include <Rcpp.h>
#include "Faddeeva.h"

// [[Rcpp::interfaces(r, cpp)]]

//' Generate 1D-NMR peak lineshape
//'
//' Calculates complex lineshapes based on chemical shift and peak parameters.
//'
//' @param direct_shift Vector of chemical shift data in ppm.
//' @param parameters n-by-4 matrix with each row corresponding to a separate
//'                   peak to be added together and the 4 columns corresponding
//'                   to position, width, height, and fraction_gauss.
//' 
//' @return complex vector
//' @export
// [[Rcpp::export]]
std::vector< std::complex<double> > lineshape_1d(
  const std::vector<double>& direct_shift,    
  const std::vector< std:vector<double> >& parameters
  ) { 

  int n = direct_shift.size();
  std::vector< std::complex<double> > intensity(N);

  for(int i = 0; i < N; i++) {
    intensity.at(i) = 1;
  }

  return intensity;
}
