#include "PriorDistrFunction_test_2nd.hpp"
#include <gsl/gsl_randist.h>
#include <cmath>

using namespace std; 

double PriorDistrFunction_test_2nd::log_pdf(const TDenseVector &x2)
{
	double log_prior_pdf = 0.0; 
	
	unsigned int count = 0; 
	log_prior_pdf += log(gsl_ran_gaussian_pdf((x2(count)-1.0), 1.0) ); 
	// Interval [-1, 1]; giota (level adjustment for ZLB drop in the interest rate)
	count ++; 
	
	if ( x2(count) < 0.0 || x2(count) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x2(count), 1.0, 2.0) ); 
	// Triangle prior on [0 1], grhomu (persistence to markup)
	count ++; 

	if ( x2(count) < 0.0 || x2(count) >= 1.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x2(count), 1.0, 2.0) ); 
	// Triangle prior on [0 1], grhorn (persistence to natural rate of interest)	
	count ++; 

	if ( x2(count) < 0.0 || x2(count) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x2(count), 1.0, 2.0) ); 
	// Triangle prior on [0 1], grhoi (persistence to monetary policy)	
	count ++; 

	if ( x2(count) < 0.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x2(count), 3.878809415314558e-01, 1.628787676246154e-04) ); 
	// Inveal [0.0001 0.5], gsigmamu (s.d. of markup shock)
	count ++; 
	
	if ( x2(count) <= 0.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x2(count), 3.878809415314558e-01, 1.628787676246154e-04) ); 
	// Inveal [0.0001 0.5], gsigmarn (s.d. of demand shock)
	count ++; 
	
	if ( x2(count) <= 0.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x2(count), 4.5897, 21.4621) ); 
	// Inveal [2.5, 12.5] for x(j_), NOT gsigmai; gsigmai = x(j_)/cl4gsigmai in ZLB (s.d. of policy shock)  
	// Note: 1) 1 (or 2.5) basis point annually = 2.5000e-05 (or 6.2500e-05) quarterly.  
	// 2) 2.5000e-04 quarterly = 10 basis points annually in the ZLB. 
	// where this has already been divided by the scale factor: makeAB_inputs.scl4gsigmai = 1.0e+05.   
	// Because it has been divided by makeAB_inputs.scl4gsigmai = 1.0e+05, we set the interval [2.5, 12.5], corresponding to 1 and 5 basis points annually.
	count ++; 
	
	// Gamma_{2,1}: can only be identified through transiitons
	
	log_prior_pdf += log(gsl_ran_gaussian_pdf((x2(count)-1.0), 1.0) ); 
	// Interval [-1, 1], Sunspot component
	count ++;
	log_prior_pdf += log(gsl_ran_gaussian_pdf((x2(count)-1.0), 1.0) ); 
	// Interval [-1, 1], Sunspot component
	count ++; 
	log_prior_pdf += log(gsl_ran_gaussian_pdf((x2(count)-1.0), 1.0) ); 
	// Interval [-1, 1], Sunspot component
	count ++;

	// p1 and p2: can only be identified through transitions
	
	if (count != x2.dim)
	{
		cerr << "*******PriorDistrFunction_test_2nd::log_pdf(): dimension inconsistent*******\n"; 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	}
	else 
		return log_prior_pdf; 
}
