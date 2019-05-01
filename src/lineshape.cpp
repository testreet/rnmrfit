#include <Rcpp.h>
#include "Faddeeva.h"

using namespace Rcpp;

//' Generate 1D-NMR peak lineshape
//'
//' Calculates complex lineshapes based on chemical shift and peak parameters.
//'
//' @param x Vector of chemical shift data in ppm.
//' @param par n-by-4 matrix with each row corresponding to a separate
//'            peak to be added together and the 4 columns corresponding
//'            to position, width, height, and fraction_gauss.
//' 
//' @return complex vector
// [[Rcpp::export]]
std::vector< std::complex<double> > lineshape_1d(
  const NumericVector x, const NumericMatrix par) { 

  int n = x.size(), nrow = par.nrow();
  std::vector< std::complex<double> > out(n, std::complex<double> (0,0));

  double p, h, wl, wg, z;
  std::complex<double> a, b;

  // Outer loop iterates over every peak
  for (int i = 0; i < nrow; i++) {

    // If the value is 0, proceed as Lorentz
    if ( par(i, 3) == 0 ) {
      // Unpacking 
      p = par(i, 0);
      wl = par(i, 1);
      h = par(i, 2);

      // Applying
      for (int j = 0; j < n; j++) {
        z = (x[j] - p)/wl;
        out.at(j) = out.at(j) + 
          (std::complex<double> (h, z)) / (std::complex<double> (1 + z*z, 0));
      }
    }
    // If the value is 1, proceed as Gauss
    else if ( par(i, 3) == 1 ) {
      p = par(i, 0);
      wg = par(i, 1);
      h = par(i, 2);

      // Applying
      for (int j = 0; j < n; j++) {
        z = (x[j] - p)/(sqrt(2)*wg);
        out.at(j) = out.at(j) + 
          h * Faddeeva::w( (std::complex<double> (z, 0)) );
      }
    }
    // Otherwise, Voigt
    else {
      // Unpacking
      p = par(i, 0);
      wl = par(i, 1);
      h = par(i, 2);
      wg = wl*par(i, 3)/(1 - par(i, 3));

      // Pre-calculating
      a = ( std::complex<double> (std::sqrt(2)*wg, 0) );
      b = ( std::complex<double> (h, 0) ) / 
          Faddeeva::w( ( std::complex<double> (0, wl) ) / a);

      // Applying
      for (int j = 0; j < n; j++) {
        out.at(j) = out.at(j) + 
          b * Faddeeva::w( (std::complex<double> (x[j] - p, wl)) / a);
      }
    }
  }

  return out;
}
