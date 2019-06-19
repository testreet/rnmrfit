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
  int n_points, n_par, n_peaks, n_baseline, n_phase, count;
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
    (*temp_fit).at(i) += ( complex<double> (h, 0) ) * fi;

    // Only compute derivatives if required
    if ( grad ) {
      dfdz = ( complex<double> (-2*z, 1-z2) ) * 
             ( complex<double> (h/(z2p1*z2p1), 0) );

      // Derivatives of z for chain rule
      dzdp = complex<double> (-1/wl, 0);
      dzdwl = complex<double> (-z/wl, 0);

      // Position gradient
      dfdp = dfdz*dzdp;
      (*temp_grad).at(i).at(j*4) = dfdp;

      // Lorentz width gradient
      dfdwl = dfdz*dzdwl;
      (*temp_grad).at(i).at(j*4+1) = dfdwl;

      // Height gradient
      (*temp_grad).at(i).at(j*4+2) = fi;

      // Fraction gradient
      (*temp_grad).at(i).at(j*4+3) = 0;
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

  // Reseting temporary fit to zero
  for(int i = 0; i < n_points; i++ ) {
    (*temp_fit).at(i) = complex<double> (0, 0);
  }
  
  // Outer loop iterates over every peak
  for (int j = 0; j < n_peaks; j++) {
    
    // If the fraction gauss value is 0, proceed as Lorentz (w becomes wl)
    // if ( par[j*4+3] == 0 ) {
    if ( true ) {
		  lorentz_1d(j, par, grad, lineshape_data);	
    }
  }
	
	// Sum of squares
	double eval = 0;
  
  for (int i = 0; i < n_points; i++) {
	  (*temp_diff).at(i) = y.at(i) - (*temp_fit).at(i); 
	  eval += (*temp_diff).at(i).real() * (*temp_diff).at(i).real();
    eval += (*temp_diff).at(i).imag() * (*temp_diff).at(i).imag();
  }

  // Gradients
  if ( grad ) {
    for (int j = 0; j < n_par; j++) {
      grad[j] = 0;
    
      for (int i = 0; i < n_points; i++) {
        grad[j] += -2*( (*temp_diff).at(i).real() * 
                        (*temp_grad).at(i).at(j).real() );

        grad[j] += -2*( (*temp_diff).at(i).imag() * 
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
  Rcpp::NumericVector par, Rcpp::NumericVector lb, Rcpp::NumericVector ub,
  int n_peaks, int n_baseline, int n_phase) {

  using namespace std;

  Rcpp::Function R_is_fin("is.finite");

  //--------------------
  // Converting and storing data
  
  int n_points = y.size();
  int n_par = par.size();

  // x and y data  
  vector< complex<double> > y_store(n_points, complex<double> (0,0));
  vector< double > x_store(n_points, 0);

  Rcpp::Function R_re("Re");
  Rcpp::Function R_im("Im");

  for (int i = 0; i < n_points; i++) {
    y_store.at(i) = complex<double> (Rcpp::as<double>(R_re(y.at(i))),
                                     Rcpp::as<double>(R_im(y.at(i))));
    x_store.at(i) = x.at(i);
  }

  // Gradient and fit terms
  vector< complex<double> > temp_fit(n_points, complex<double> (0,0));
  vector< complex<double> > temp_diff(n_points, complex<double> (0,0));
  vector< vector < complex<double> > > 
    temp_grad(n_points, vector< std::complex<double> > 
        ( n_par, std::complex<double> (0,0) ) 
    );

  // Packing lineshape structure
  stored_data lineshape_data[11] = {
    x_store, y_store, temp_fit, temp_diff, temp_grad, 
    n_points, n_par, n_peaks, n_baseline, n_phase, 0};

  //--------------------
  // Setting up the optimizer

  nlopt_opt opt;

  opt = nlopt_create(NLOPT_LD_SLSQP, n_par);    
  //opt = nlopt_create(NLOPT_LN_COBYLA, n_par); 
  nlopt_set_min_objective(opt, f_obj_1d, lineshape_data);
  nlopt_set_xtol_rel(opt, 1e-4);

  // bounds
  nlopt_set_lower_bounds(opt, &(lb[0]));
  nlopt_set_upper_bounds(opt, &(ub[0]));
    
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




       
