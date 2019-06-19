#include <vector>
#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>
#include <nloptrAPI.h>
#include "Faddeeva.h"



//' @importFrom Rcpp sourceCpp
//' @import nloptr
//' @useDynLib rnmrfit



//==============================================================================
// Lineshape and constraint structures
typedef struct {
  std::vector< double > x;	
  std::vector< std::complex<double> > y;	
  std::vector< std::complex<double> > temp_fit;
  std::vector< std::complex<double> > temp_diff;
  std::vector< std::vector < std::complex<double> > > temp_grad;
  int n_points, n_par, n_peaks, count;
} stored_data;


//==============================================================================
// Lineshape functions

//------------------------------------------------------------------------------
void lorentz_1d(int j, const double *par, double *grad, void *lineshape_data) {

	using namespace std;

  // Unpacking lineshape_data structure
  stored_data *sd = (stored_data *) lineshape_data;
  vector< double > x = sd->x;
  vector< complex<double> > *temp_fit = &(sd->temp_fit);
  vector< vector < complex<double> > > *temp_grad = &(sd->temp_grad);
  int n_points = sd->n_points;
  int n_peaks = sd->n_peaks;

  // Extracting parameters
  double p = par[j*4];
  double wl = par[j*4+1];
  double h = par[j*4+2];

  // Main terms
  complex<double> fi, dfdz, dfdp, dfdwl, dzdwl, dzdp;

  // Commonly used intermediate values
  double z, z2, z2p1;

  // Looping through each x value
  for(int i = 0; i < n_points; i++ ) {

    // Pre-calculating common values
    z = (x.at(i)-p)/wl;
    z2 = z*z;
    z2p1 = z2 + 1;

    // Calculating an intermediate f with no h term
    fi = ( complex<double> (1, z) ) / 
         ( complex<double> (z2p1, 0) );
       
    // Tacking on the h
    (*temp_fit).at(i) = ( complex<double> (h, 0) ) * fi;

    // Only compute derivatives if required
    if ( grad ) {
      dfdz = ( complex<double> (-2*z, 1-z2) ) * 
             ( complex<double> (h/(z2p1*z2p1), 0) );

      // Derivatives of z for chain rule
      dzdp = complex<double> (-1/wl, 0);
      dzdwl = complex<double> (-z/wl, 0);

      // Position gradient
      dfdp = dfdz*dzdp;
      (*temp_grad).at(i).at(0) = dfdp;

      // Lorentz width gradient
      dfdwl = dfdz*dzdwl;
      (*temp_grad).at(i).at(1) = dfdwl;

      // Height gradient
      (*temp_grad).at(i).at(2) = fi;

      // Fraction gradient
      (*temp_grad).at(i).at(3) = 0;
    }
  }

  return; 
}



//===============================================================================
// NLOPT objective function

double f_obj_1d(unsigned n, const double *par, double *grad, void *lineshape_data) {
    
	using namespace std;
  
  // Unpacking lineshape_data structure
  stored_data *sd = (stored_data *) lineshape_data;
  vector< double > x = sd->x;
  vector< complex<double> > y = sd->y;
  vector< complex<double> > *temp_fit = &(sd->temp_fit);
  vector< complex<double> > *temp_diff = &(sd->temp_diff);
  vector< vector < complex<double> > > *temp_grad = &(sd->temp_grad);
  int n_points = sd->n_points;
  int n_par = sd->n_par;
  int n_peaks = sd->n_peaks;
  int *count = &(sd->count);

  *count += 1;

  //--------------------
  // Calculate lineshape temp_fit, and gradients grad
	
  // Outer loop iterates over every peak
  for (int j = 0; j < n_peaks; j++) {

    // If the fraction gauss value is 0, proceed as Lorentz (w becomes wl)
    if ( par[j*4+3] == 0 ) {
		  lorentz_1d(j, par, grad, lineshape_data);	
    }
  }
	
	// Sum of squares
	double eval = 0;
  
  for (int i = 0; i < n_points; i++) {
	  (*temp_diff).at(i) = y.at(i) - (*temp_fit).at(i); 
	  eval += (*temp_diff).at(i).real() * (*temp_diff).at(i).real() +
            (*temp_diff).at(i).imag() * (*temp_diff).at(i).imag();
  }

  // Gradients
  if ( grad ) {
    for (int j = 0; j < n_par; j++) {
      grad[j] = 0;
    
      for (int i = 0; i < n_points; i++) {
        grad[j] += -2*( (*temp_diff).at(i).real() * 
                        (*temp_grad).at(i).at(j).real() + 
                        (*temp_diff).at(i).imag() * 
                        (*temp_grad).at(i).at(j).imag() );
      }
    }
  }

  return eval;	
}



//===========================================================================
//Main function called from R

//' @export
// [[Rcpp::export]]
double fit_lineshape_1d(
       const Rcpp::NumericVector x, const Rcpp::ComplexVector y,
       Rcpp::NumericVector par) {

  using namespace std;

  //--------------------
  // Converting and storing data
  
  // Defining some common lengths to avoid confustion
  int n_points = y.size();
  int n_par = par.size();
  int n_peaks = par.size()/4;

  // note that while n_peaks looks a little weird now, it should preempt future
  // issues with baseline/phase terms

  // x and y data  
  vector< complex<double> > y_store(n_points, complex<double> (0,0));
  vector< double > x_store(n_points, 0);

  Rcpp::Function Re("Re");
  Rcpp::Function Im("Im");

  for (int i = 0; i < n_points; i++) {
    y_store.at(i) = complex<double> (Rcpp::as<double>(Re(y.at(i))),
                                     Rcpp::as<double>(Im(y.at(i))));
    x_store.at(i) = x.at(i);
  }

  // Gradient and fit terms
  vector< complex<double> > temp_fit(n_points, complex<double> (0,0));
  vector< complex<double> > temp_diff(n_points, complex<double> (0,0));
  vector< vector < complex<double> > > 
    temp_grad(n_points, vector< std::complex<double> > 
        ( n_peaks*4, std::complex<double> (0,0) ) 
    );

  // Packing lineshape structure
  stored_data lineshape_data[9] = {
    x_store, y_store, temp_fit, temp_diff, temp_grad, 
    n_points, n_par, n_peaks, 0};

  //--------------------
  // Setting up the optimizer

  nlopt_opt opt;

  opt = nlopt_create(NLOPT_LD_SLSQP, 4);    
  //opt = nlopt_create(NLOPT_LN_COBYLA, 4); 
  nlopt_set_min_objective(opt, f_obj_1d, lineshape_data);
  nlopt_set_xtol_rel(opt, 1e-4);

  // bounds
  double lb[4] = { -1, 1e-6, 0, 0 };
  double ub[4] = { 1, 1, 2, 0 };
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);
    
  // Minimum objective value, upon return
  double minf; 						

  //--------------------
  // Running it

  stored_data *sd = (stored_data *) lineshape_data;
  int *count = &(sd->count);

  if ( nlopt_optimize(opt, &(par[0]), &minf ) < 0) {
    Rcpp::Rcout << "nlopt failed!" << std::endl;
  } else {
    Rcpp::Rcout << "nlopt succeeded!" << std::endl;
    	//Checking minf
	Rcpp::Rcout << "Found minimum with" << std::endl;  
    Rcpp::Rcout << "par(" << par[0] << ", " << par[1] << ", " << par[2] << ", " << par[3] << ")" << std::endl;
  }
	Rcpp::Rcout << "Iterations: " << *count << std::endl;  
  nlopt_destroy(opt);

	return minf;
}




       
