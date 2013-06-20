#include "PriorDistrFunction_test_2nd.hpp"
#include <gsl/gsl_randist.h>
#include <cmath>

using namespace std; 

/*
 * Parameters of gsl_ran_gamma_pdf are not the same as the parameters of fn_logpdf_gamma
 *      In gsl_ran_gamma_pdf, p(x) = 1/(Gamma(a) * b^a) x^{a-1} e^{-x/b}
 *      In fn_logpdf_gamma, p(x) = b^a/Gamma(a) x^{a-1} e^{-bx}
 *      That is, the two b parameters are reciprocal.
 *
 * Inverse gamma vs. Gamma distribution
 *      In fn_logpdf_gamma, p1(x) = b^a/Gamma(a) x^{a-1} e^{-bx}
 *      In fn_logpdf_invgam, p2(x) = b^a/Gamma(a) x^{-a-1} e^{-b/x}
 *      p2(x) = p1(1/x) / x^2
 *      That is, log(p2(x)) = log(p1(1/x))-2*log(x)
 *
 * Using gsl_ran_gamma_pdf to implement fn_logpdf_invgam
 *      fn_logpdf_invgam(x, a, b) <=> log(gsl_ran_gamma_pdf(1.0/x, a, b))-2.0*log(x)
 *
 *
 */

double PriorDistrFunction_test_2nd::log_pdf(const TDenseVector &x2)
{
	double log_prior_pdf = 0.0; 
	
	unsigned int count = 0; 
	// Per Dan's suggestion on 06/18/2013, giota ~ N(0,1.0)
	log_prior_pdf += log(gsl_ran_gaussian_pdf(x2(count), 1.0) );
	// log_prior_pdf += log(gsl_ran_gaussian_pdf((x2(count)-1.0), 1.0) ); 
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
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x2(count), 3.878809415314558e-01, 1.0/1.628787676246154e-04) )-2.0*log(x2(count)); 
	// Inveal [0.0001 0.5], gsigmamu (s.d. of markup shock)
	count ++; 
	
	if ( x2(count) <= 0.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x2(count), 3.878809415314558e-01, 1.0/1.628787676246154e-04) )-2.0*log(x2(count)); 
	// Inveal [0.0001 0.5], gsigmarn (s.d. of demand shock)
	count ++; 
	
	if ( x2(count) <= 0.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x2(count), 4.5897, 1.0/21.4621) )-2.0*log(x2(count)); 
	// Inveal [2.5, 12.5] for x(j_), NOT gsigmai; gsigmai = x(j_)/cl4gsigmai in ZLB (s.d. of policy shock)  
	// Note: 1) 1 (or 2.5) basis point annually = 2.5000e-05 (or 6.2500e-05) quarterly.  
	// 2) 2.5000e-04 quarterly = 10 basis points annually in the ZLB. 
	// where this has already been divided by the scale factor: makeAB_inputs.scl4gsigmai = 1.0e+05.   
	// Because it has been divided by makeAB_inputs.scl4gsigmai = 1.0e+05, we set the interval [2.5, 12.5], corresponding to 1 and 5 basis points annually.
	count ++; 
	
	// Gamma_{2,1}: can only be identified through transiitons

	// Per Dan's suggestions on 06/18/2013, all sunspot component parameters ~ N(0,5.0)	
	// Per Dan's suggestions on 06/19/2013, all sunspot component parameters ~ N(0,1.0)
	log_prior_pdf += log(gsl_ran_gaussian_pdf(x2(count), 1.0) );
	// log_prior_pdf += log(gsl_ran_gaussian_pdf((x2(count)-1.0), 1.0) ); 
	// Interval [-1, 1], Sunspot component
	count ++;
	log_prior_pdf += log(gsl_ran_gaussian_pdf(x2(count), 1.0) );
	// log_prior_pdf += log(gsl_ran_gaussian_pdf((x2(count)-1.0), 1.0) ); 
	// Interval [-1, 1], Sunspot component
	count ++; 
	log_prior_pdf += log(gsl_ran_gaussian_pdf(x2(count), 1.0) );
	// log_prior_pdf += log(gsl_ran_gaussian_pdf((x2(count)-1.0), 1.0) ); 
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
