#include <vector>
#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>
#include <nloptrAPI.h>
#include "Faddeeva.h"
#include "rnmrfit.h"


//' @importFrom Rcpp sourceCpp
//' @import nloptr
//' @useDynLib rnmrfit



//==============================================================================
// Constraint functions
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
double constrain_position(unsigned n, const double *x, 
                          double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_constraint *d = (data_constraint *) data;
  
  vector< bool > sign = d->sign;	
  vector< int > peak_number = d->peak_number_1;	
  double offset = d->offset;

  // Ensuring that the gradient is zero for all unaffected terms
  if ( grad ) {
    for (int i = 0; i < n; i++) {
      grad[i] = 0;
    }
  }

  // Evaluating
  double eval = 0;

  for (int i = 0; i < peak_number.size(); i++) {

    // Converting peak number to index
    int j = (peak_number.at(i) - 1)*4;

    // Extracting out sign
    bool neg = sign.at(i);

    if ( neg ) { eval -= x[j]; }
    else { eval += x[j]; }

    if (grad) {
      if ( neg ) { grad[j] = -1; }
      else { grad[j] = 1; }
    }
  }
  
  return (eval - offset);
}



//------------------------------------------------------------------------------
double constrain_width(unsigned n, const double *x, 
                      double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_constraint *d = (data_constraint *) data;
  
  vector< bool > sign = d->sign;	
  vector< int > peak_number = d->peak_number_1;	
  double offset = d->offset;

  // Ensuring that the gradient is zero for all unaffected terms
  if ( grad ) {
    for (int i = 0; i < n; i++) {
      grad[i] = 0;
    }
  }

  // Unlike position the width constraints divides the sum
  // of all "positive" indexes by the sum of all "negative" indexes so
  // it's necessary to do 2 passes: one calculate the evaluation and another
  // to calculate the gradients.

  // First pass to calculate evaluation
  double numerator = 0;
  double denominator = 0;

  for (int i = 0; i < peak_number.size(); i++) {

    // Converting peak number to index of width
    int j = (peak_number.at(i) - 1)*4 + 1;

    // Extracting out sign
    bool neg = sign.at(i);

    if ( neg ) { denominator += x[j]; }
    else { numerator += x[j]; }
  }

  double eval = numerator/denominator;

  // Second pass to calculate gradients
  if ( grad ) {

    // All gradient terms are the same
    double grad_pos = 1/denominator;
    double grad_neg = -eval/denominator;

    for (int i = 0; i < peak_number.size(); i++) {

      // Converting peak number to index of width
      int j = (peak_number.at(i) - 1)*4 + 1;

      // Extracting out sign
      bool neg = sign.at(i);

      if ( neg ) { grad[j] = grad_neg; }
      else {  grad[j] = grad_pos; }
    }
  }
  
  return (eval - offset);
}

//------------------------------------------------------------------------------
double constrain_area(unsigned n, const double *x, 
                      double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_constraint *d = (data_constraint *) data;
  
  vector< bool > sign = d->sign;	
  vector< int > peak_number = d->peak_number_1;	
  double offset = d->offset;

  // Ensuring that the gradient is zero for all unaffected terms
  if ( grad ) {
    for (int i = 0; i < n; i++) {
      grad[i] = 0;
    }
  }

  // Unlike position the area constraints divides the sum
  // of all "positive" indexes by the sum of all "negative" indexes so
  // it's necessary to do 2 passes: one calculate the evaluation and another
  // to calculate the gradients.

  // First pass to calculate evaluation
  double numerator = 0;
  double denominator = 0;

  for (int i = 0; i < peak_number.size(); i++) {

    // Converting peak number to index
    int j = (peak_number.at(i) - 1)*4;

    double w = x[j+1];
    double h = x[j+2];
    double f = x[j+3];
    double area = 0;

    // Area calculation depends on fraction
    if ( f < 1e-6 ) { area = M_PI * w * h; }
    else { Rcpp::stop("Gauss and Voigt peaks are not currently supported"); }

    // Extracting out sign
    bool neg = sign.at(i);

    if ( neg ) { denominator += area; }
    else { numerator += area; }
  }

  double eval = numerator/denominator;

  // Second pass to calculate gradients
  if ( grad ) {

    // Calculating some common terms
    double dfda_pos = 1/denominator;
    double dfda_neg = -eval/denominator;
    double dfda = 0;

    for (int i = 0; i < peak_number.size(); i++) {

      // Converting peak number to index
      int j = (peak_number.at(i) - 1)*4;

      double w = x[j+1];
      double h = x[j+2];
      double f = x[j+3];

      // Extracting out sign
      bool neg = sign.at(i);

      if ( neg ) { dfda = dfda_neg; }
      else { dfda = dfda_pos; }

      // Gradient calculations depend on fraction
      if ( f < 1e-6 ) {
        // Width
        grad[j+1] = dfda * M_PI * h;
        // Height
        grad[j+2] = dfda * M_PI * w;
      } else { 
        Rcpp::stop("Gauss and Voigt peaks are not currently supported"); 
      }
    }
  }
  
  return (eval - offset);
}




//==============================================================================
// Lineshape functions
//------------------------------------------------------------------------------
// All lineshape functions take a series of peak parameters and use those
// parameters to increment fit data in a data_lineshape structure. To reuse the
// same functions for 1D and 2D analysis, the fits are divided by resonance
// index, which should be kept at 0 for 1D usage.



//------------------------------------------------------------------------------
void lorentz(double p, double wl, double h, 
             int i_res, int i_par, double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_lineshape *d = (data_lineshape *) data;
  
  vector< double > x = d->x;
  
  vector< complex<double> > *temp_fit = &(d->temp_fit.at(i_res));
  vector< vector < complex<double> > > *temp_grad = &(d->temp_grad);

  // Main terms
  complex<double> fo, dfdz, dfdp, dfdwl, dzdwl, dzdp;

  // Commonly used intermediate values
  double z, z2, z2p1;

  // Looping through each x value
  for(int i = 0; i < x.size(); i++ ) {

    // Pre-calculating common values
    z = (x.at(i)-p)/wl;
    z2 = z*z;
    z2p1 = z2 + 1;

    // Calculating an intermediate f with no h term
    fo = ( complex<double> (1, z) ) / 
         ( complex<double> (z2p1, 0) );
       
    // Tacking on the h
    (*temp_fit).at(i) += ( complex<double> (h, 0) ) * fo;

    // Only compute derivatives if required
    if ( grad ) {
      dfdz = ( complex<double> (-2*z, 1-z2) ) * 
             ( complex<double> (h/(z2p1*z2p1), 0) );

      // Derivatives of z for chain rule
      dzdp = complex<double> (-1/wl, 0);
      dzdwl = complex<double> (-z/wl, 0);

      // Position gradient
      dfdp = dfdz*dzdp;
      (*temp_grad).at(i).at(i_par) = dfdp;

      // Lorentz width gradient
      dfdwl = dfdz*dzdwl;
      (*temp_grad).at(i).at(i_par+1) = dfdwl;

      // Height gradient
      (*temp_grad).at(i).at(i_par+2) = fo;

      // Fraction gradient
      (*temp_grad).at(i).at(i_par+3) = 0;
    }
  }

  return; 
}

//==============================================================================
// 1D fitting

//------------------------------------------------------------------------------
// Objective function

double f_obj_1d(unsigned n, const double *par, double *grad, void *data) {
    
	using namespace std;

  // Unpacking lineshape_data structure
  data_1d *d = (data_1d *) data;

  vector< vector<double> > y = d->y;
  vector< vector<double> > *y_diff = &(d->y_diff);
  int n_par = d->n_par;
  int n_peaks = d->n_peaks;
  int *count = &(d->count);

  data_lineshape *data_direct = &(d->lineshape);
  data_lineshape *d1 = (data_lineshape *) data_direct;

  vector< complex<double> > *temp_fit = &(d1->temp_fit.at(0));
  vector< vector < complex<double> > > *temp_grad = &(d1->temp_grad);

  *count += 1;

  //--------------------
  // Calculate lineshape temp_fit, and gradients grad

  // Reseting temporary fit to zero
  for(int i = 0; i < (*temp_fit).size(); i++ ) {
    (*temp_fit).at(i) = complex<double> (0, 0);
  }
  
  // Loop over peaks
  for (int i = 0; i < n_peaks; i++) {

    double p = par[i*4];
    double wl = par[i*4+1];
    double h = par[i*4+2];
    double f = par[i*4+3];
    
    // If the fraction gauss value is 0, proceed as Lorentz (w becomes wl)
    if ( true ) {
		  lorentz(p, wl, h, 0, i*4, grad, data_direct);	
    }
  }
	
	// Sum of squares
	double eval = 0;
  
  // For 1D fit, assume that indexes of y and lineshape data match up
  for (int i = 0; i < y.size(); i++) {
	  (*y_diff).at(i).at(0) = y.at(i).at(0) - (*temp_fit).at(i).real(); 
	  eval += (*y_diff).at(i).at(0) * (*y_diff).at(i).at(0);

	  (*y_diff).at(i).at(1) = y.at(i).at(1) - (*temp_fit).at(i).imag(); 
	  eval += (*y_diff).at(i).at(1) * (*y_diff).at(i).at(1);
  }

  // Gradients
  if ( grad ) {
    for (int j = 0; j < n_par; j++) {
      grad[j] = 0;
    
      for (int i = 0; i < y.size(); i++) {
        grad[j] += -2*( (*y_diff).at(i).at(0) * 
                        (*temp_grad).at(i).at(j).real() );

        grad[j] += -2*( (*y_diff).at(i).at(1) * 
                        (*temp_grad).at(i).at(j).imag() );
      }
    }
  }

  return eval;	
}



//------------------------------------------------------------------------------
// Main function called from R

//' @export
// [[Rcpp::export]]
double fit_lineshape_1d(
  const Rcpp::NumericVector x_val, const Rcpp::ComplexVector y_val,
  Rcpp::NumericVector par, Rcpp::NumericVector lb, Rcpp::NumericVector ub,
  Rcpp::List eq, Rcpp::List ineq, int n_peaks, int n_baseline, int n_phase) {

  using namespace std;

  //--------------------
  // Converting and storing data
  
  int n_points = y_val.size();
  int n_par = par.size();

  // x and y data  
  vector< vector<double> > y(n_points, vector<double> (2,0) );
  vector< vector<double> > y_fit(n_points, vector<double> (2,0) );
  vector< vector<double> > y_diff(n_points, vector<double> (2,0) );
  vector< double > x(n_points, 0);

  Rcpp::Function R_re("Re");
  Rcpp::Function R_im("Im");

  // For 1D fit, x and y can be treated as paired
  for (int i = 0; i < n_points; i++) {
    y.at(i).at(0) = Rcpp::as<double>(R_re(y_val.at(i)));
    y.at(i).at(1) = Rcpp::as<double>(R_im(y_val.at(i)));
    x.at(i) = x_val.at(i);
  }

  // First, packing lineshape data
  data_lineshape data_direct;
  
  data_direct.x = x;

  vector< vector < complex<double> > > 
    direct_temp_fit( 1, vector< complex<double> > 
                     ( n_points, complex<double> (0,0) ) );
  data_direct.temp_fit = direct_temp_fit;
  
  vector< vector < complex<double> > > 
    direct_temp_grad( n_points, vector< complex<double> > 
                      ( n_par, complex<double> (0,0) ) );
  data_direct.temp_grad = direct_temp_grad;

  // And then packing general data
  data_1d data[1];

  data[0].lineshape = data_direct;
  data[0].y = y;
  data[0].y_fit = y_fit;
  data[0].y_diff = y_diff;
  data[0].n_par = n_par;
  data[0].n_peaks = n_peaks;
  data[0].n_baseline = n_baseline;
  data[0].n_phase = n_phase;
  data[0].count = 0;

  //--------------------
  // Initializing the optimizer

  nlopt_opt opt;

  opt = nlopt_create(NLOPT_LD_SLSQP, n_par);    
  //opt = nlopt_create(NLOPT_LN_COBYLA, n_par); 
  nlopt_set_min_objective(opt, f_obj_1d, &data[0]);
  nlopt_set_xtol_rel(opt, 1e-4);

  // Bounds
  nlopt_set_lower_bounds(opt, &(lb[0]));
  nlopt_set_upper_bounds(opt, &(ub[0]));

  //--------------------
  // Adding constraints

  // First the equality constraints
  int n_eq = eq.size();
  data_constraint data_eq[n_eq];

  for (int i = 0; i < n_eq; i++) {

    vector< double > constraints = eq.at(i); 
    int n_terms = constraints.size() - 2;

    data_eq[i].offset = constraints.at(1);

    data_eq[i].sign = vector< bool > ( n_terms, 0);
    data_eq[i].peak_number_1 = vector< int > ( n_terms, 0);

    for (int j = 0; j < n_terms; j++) {
      bool neg = signbit(constraints.at(j+2));
      data_eq[i].sign[j] = neg;
      if ( neg ) { data_eq[i].peak_number_1[j] = -round(constraints.at(j+2)); }
      else { data_eq[i].peak_number_1[j] = round(constraints.at(j+2)); }
    }

    // Adding constraint based on flag
    int flag = round(constraints.at(0));
    if ( flag == 0 ) {
      nlopt_add_equality_constraint(opt, constrain_position, &data_eq[i], 1e-8);
    } else if ( flag == 1) {
      nlopt_add_equality_constraint(opt, constrain_width, &data_eq[i], 1e-8);
    } else if ( flag == 2) {
      nlopt_add_equality_constraint(opt, constrain_area, &data_eq[i], 1e-8);
    } else {
      Rcpp::stop("Encountered undefined equality constraint. Aborting.");
    }
  }

  // The inequality constraint proceed exactly the same as the equality ones
  int n_ineq = ineq.size();
  data_constraint data_ineq[n_ineq];

  for (int i = 0; i < n_ineq; i++) {

    vector< double > constraints = ineq.at(i); 
    int n_terms = constraints.size() - 2;

    data_ineq[i].offset = constraints.at(1);

    data_ineq[i].sign = vector< bool > ( n_terms, 0);
    data_ineq[i].peak_number_1 = vector< int > ( n_terms, 0);

    for (int j = 0; j < n_terms; j++) {
      bool neg = signbit(constraints.at(j+2));
      data_ineq[i].sign[j] = neg;
      if ( neg ) { data_ineq[i].peak_number_1[j] = -round(constraints.at(j+2));}
      else { data_ineq[i].peak_number_1[j] = round(constraints.at(j+2)); }
    }

    // Adding constraint based on flag
    int flag = round(constraints.at(0));
    if ( flag == 0 ) {
      nlopt_add_inequality_constraint(opt, constrain_position, 
                                      &data_ineq[i], 1e-8);
    } else if ( flag == 1) {
      nlopt_add_inequality_constraint(opt, constrain_width, 
                                      &data_ineq[i], 1e-8);
    } else if ( flag == 2) {
      nlopt_add_inequality_constraint(opt, constrain_area, 
                                      &data_ineq[i], 1e-8);
    } else {
      Rcpp::stop("Encountered undefined inequality constraint. Aborting.");
    }
  }
    
  // Minimum objective value, upon return
  double minf; 						

  //--------------------
  // Running it

  data_1d *d = (data_1d *) data;
  int *count = &(d->count);

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




       
