#include "PriorDistrFunction_test_1st.hpp"
#include <cmath>
#include <gsl/gsl_randist.h>

/*
 * Parameters of gsl_ran_gamma_pdf are not the same as the parameters of fn_logpdf_gamma
 * 	In gsl_ran_gamma_pdf, p(x) = 1/(Gamma(a) * b^a) x^{a-1} e^{-x/b}
 * 	In fn_logpdf_gamma, p(x) = b^a/Gamma(a) x^{a-1} e^{-bx}
 * 	That is, the two b parameters are reciprocal.
 */

using namespace std; 

double PriorDistrFunction_test_1st::log_pdf(const TDenseVector &x1)
{
	double log_prior_pdf = 0; 
	unsigned int i_ = 0; 

	if (x1(i_) <= 0.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(x1(i_), 6.0471, 1.0/1057.6212) ); 
		// Interval (0.0025, 0.01); Steady state quarterly rate of inflation (decimal) (e.g. 0.5/100 (0.005) quarterly corresponding to log(1.02) or 2% annually). 
	i_ ++; 

	if(x1(i_) <= 0.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(x1(i_), 9.3869, 1.0/995.1952) );
		// Interval (0.005, 0.015); Steady state quarterly nominal interest rate, corresponding to i in Zha's notes (demical) (e.g., 0.015 quarterly corresponds to 6% annually).
	i_ ++; 

	if (x1(i_) <= 0.0 || x1(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x1(i_), 2.0, 2.0) ); 
		// Interval (0.14, 0.86); gthetaf (IS coefficient on expected output)
	i_ ++; 

	if (x1(i_) <= 0.0 || x1(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf +=  log(gsl_ran_beta_pdf(x1(i_), 2.0, 2.0) ); 
		// Interval (0.14, 0.86); gtau (IS coefficient on real interest rate)
	i_ ++;

	if (x1(i_) <= 0.0 || x1(i_) >= 1.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x1(i_), 2.0, 2.0) ); 
		// Interval (0.14, 0.86); gomegaf (AS coefficient on expected inflation)
	i_ ++; 

	if (x1(i_) <= 0.0 || x1(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 	
		log_prior_pdf += log(gsl_ran_beta_pdf(x1(i_), 2.0, 2.0) ); 
		// Interval (0.14, 0.86); gkappac (slope of Phillips curve on current output)
	i_ ++; 

	if (x1(i_) <= 0.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(x1(i_), 2.294989312971724e+01,1.0/1.045201956680500e+01) ); 
		// Interval [1.5 3]; gphipi (policy response to inflation)
	i_ ++; 

	if (x1(i_) <= 0.0 || x1(i_) >= 1.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else
		log_prior_pdf += log(gsl_ran_beta_pdf(x1(i_), 2.0, 2.0) ); 
		// Interval (0.14, 0.86); gphiy (plicy response to output gap)
	i_ ++; 			 

	if (x1(i_) <= 0.0 || x1(i_) >= 1.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x1(i_), 1.0, 2.0) ); 	
		// Triangle prior on [0 1]; grhomu (persistence to markup) 
	i_ ++; 

	if (x1(i_) <= 0.0 || x1(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x1(i_), 1.0, 2.0) ); 
		// Triangle prior on [0 1]; grhorn (persistence to natural rate of interest)
	i_ ++; 

	if (x1(i_) <= 0.0 || x1(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x1(i_), 1.0, 2.0) ); 
		// Triangle prior on [0 1]; grhoi (persistence to monetary policy)
	i_ ++; 

	// Gamma distribution vs. Inverse gamma distribution
	// 	X ~ Gamma(alpha, beta)
	// 	Y=1/X ~ InverseGamma(alpha, 1/beta) 
	if (x1(i_) <= 0.0 )
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x1(i_), 3.878809415314558e-01, 1.628787676246154e-04) ); 
		// Inveal [0.0001 0.5]; gsigmamu (s.d. of markup shock)
	i_ ++; 
	
	if (x1(i_) <= 0.0 )
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x1(i_), 3.878809415314558e-01, 1.628787676246154e-04) ); 
		// Inveal [0.0001 0.5]; gsigmarn (s.d. of demand shock)
	i_ ++; 

	if (x1(i_) <= 0.0 )
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x1(i_), 3.878809415314558e-01, 1.628787676246154e-04) );
		// Inveal [0.0001 0.5]; gsigmai (s.d. of policy shock)
	i_ ++; 

	// double check 
	if (i_ != x1.dim)
	{
		cerr << "******* PiorDistrFunction_test_1st::log_pdf(): dimension inconsist *******\n"; 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	}
	else 
		return log_prior_pdf; 
}		

