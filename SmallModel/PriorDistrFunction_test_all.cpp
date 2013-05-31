#include "PriorDistrFunction_test_all.hpp"
#include <gsl/gsl_randist.h>
#include <cmath>

using namespace std; 

double PriorDistrFunction_test_all::log_pdf(const TDenseVector &x)
{
	double log_prior_pdf = 0.0; 
	size_t n_x = 15; 
	unsigned int i_ = 0, j_ = n_x; 

	if (x(i_) <= 0.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf +=  log(gsl_ran_gamma_pdf(x(i_), 6.0471,1.0/1057.6212) ); 
	// Interval (0.0025, 0.01); Steady state quarterly rate of inflation (decimal) (e.g. 0.5/100 (0.005) quarterly corresponding to log(1.02) or 2% annually)
	i_ ++; 

	if (x(i_) <= 0.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else
		log_prior_pdf += log(gsl_ran_gamma_pdf(x(i_), 9.3869,1.0/995.1952) ); 
	// Interval (0.005, 0.015); Steady state quarterly nominal interest rate, corresponding to i in Zha's notes (demical) (e.g., 0.015 quarterly corresponds to 6% annually).
	i_ ++; 

	if (x(i_) <= 0.0 || x(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x(i_), 2.0, 2.0) ); 
	// Interval (0.14, 0.86); gthetaf (IS coefficient on expected output)
	i_ ++; 

	if (x(i_) <= 0.0 || x(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x(i_), 2.0, 2.0) ); 
	// Interval (0.14, 0.86); gtau (IS coefficient on real interest rate)
	i_ ++; 

	if (x(i_) <= 0.0 || x(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x(i_), 2.0, 2.0) ); 
	// Interval (0.14, 0.86); gomegaf (AS coefficient on expected inflation)
	i_ ++; 	
	
	if (x(i_) <= 0.0 || x(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x(i_), 2.0, 2.0) ); 
	// Interval (0.14, 0.86); gkappac (slope of Phillips curve on current output)
	i_ ++; 

	if (x(i_) <= 0.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf +=  log(gsl_ran_gamma_pdf(x(i_), 2.294989312971724e+01, 1.0/1.045201956680500e+01) ); 
	// Interval [1.5 3]; gphipi (policy response to inflation)
	i_ ++; 

	if (x(i_) <= 0.0 || x(i_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
		log_prior_pdf += log(gsl_ran_beta_pdf(x(i_), 2.0, 2.0) ); 
	// Interval (0.14, 0.86); gphiy (plicy response to output gap)
	i_ ++; 

	if (x(i_) < 0.0 || x(i_) >= 1.0 || x(j_) < 0.0 || x(j_) >= 1.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
	{
		log_prior_pdf +=  log(gsl_ran_beta_pdf(x(i_), 1.0, 2.0) ); 
		// Triangle prior on [0 1]; grhomu (persistence to markup)
		log_prior_pdf +=  log(gsl_ran_beta_pdf(x(j_), 1.0, 2.0) ); 
		// Triangle prior on [0 1]; grhomu (persistence to markup)
	}
	i_ ++; 
	j_ ++; 

	if (x(i_) < 0.0 || x(i_) >= 1.0 || x(j_) < 0.0 || x(j_) >= 1.0)
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
	{
		log_prior_pdf += log(gsl_ran_beta_pdf(x(i_), 1.0, 2.0) ); 
		// Triangle prior on [0 1]; grhorn (persistence to natural rate of interest)
		log_prior_pdf += log(gsl_ran_beta_pdf(x(j_), 1.0, 2.0) ); 
		// Triangle prior on [0 1]; grhorn (persistence to natural rate of interest)
	}
	i_ ++; 
	j_ ++; 

	if (x(i_) < 0.0 || x(i_) >= 1.0 || x(j_) < 0.0 || x(j_) >= 1.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
	{
		log_prior_pdf += log(gsl_ran_beta_pdf(x(i_), 1.0, 2.0) ); 
		// Triangle prior on [0 1]; grhoi (persistence to monetary policy)
		log_prior_pdf += log(gsl_ran_beta_pdf(x(j_), 1.0, 2.0) ); 
		// Triangle prior on [0 1]; grhoi (persistence to monetary policy)
	}
	i_ ++; 
	j_ ++; 

	if (x(i_) <= 0.0 || x(j_) <= 0.0 )
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
	{
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x(i_), 3.878809415314558e-01, 1.628787676246154e-04) ); 
		// Inveal [0.0001 0.5]; gsigmamu (s.d. of markup shock)
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x(j_), 3.878809415314558e-01, 1.628787676246154e-04) ); 
		// Inveal [0.0001 0.5]; gsigmamu (s.d. of markup shock)
	}	
	i_ ++; 
	j_ ++; 

	if (x(i_) <= 0.0 || x(j_) <= 0.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
	{
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x(i_), 3.878809415314558e-01, 1.628787676246154e-04) ); 
		// Inveal [0.0001 0.5]; gsigmarn (s.d. of demand shock)
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x(j_), 3.878809415314558e-01, 1.628787676246154e-04) ); 
		// Inveal [0.0001 0.5]; gsigmarn (s.d. of demand shock)
	}
	i_ ++; 
	j_ ++; 

	if (x(i_) <= 0.0 || x(j_) <= 0.0) 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	else 
	{
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x(i_), 3.878809415314558e-01, 1.628787676246154e-04) ); 
		// Inveal [0.0001 0.5]; gsigmai (s.d. of policy shock)
		log_prior_pdf += log(gsl_ran_gamma_pdf(1.0/x(j_), 4.5897, 21.4621) ); 
		// Inveal [2.5, 12.5] for x(j_), NOT gsigmai; gsigmai = x(j_)/scl4gsigmai in ZLB (s.d. of policy shock)
		// Note: 1) 1 (or 2.5) basis point annually = 2.5000e-05 (or 6.2500e-05) quarterly.  2) 2.5000e-04 quarterly = 10 basis points annually in the ZLB.
		// where this has already been divided by the scale factor: makeAB_inputs.scl4gsigmai = 1.0e+05.
		// Because it has been divided by makeAB_inputs.scl4gsigmai = 1.0e+05, we set the interval [2.5, 12.5], corresponding to 1 and 5 basis points annually.
	}	
	i_ ++; 
	j_ ++; 

	log_prior_pdf += log(gsl_ran_gaussian_pdf((x(i_)-1.0), 1.0) ); 
	// Interval [-1, 1]; giota (level adjustment for ZLB drop in the interest rate)
	i_ ++; 	

	// Gamma_{2,1}: Can be only identified through transitions.  The first three is not useful because there is only one endogenous error for this case.
	
	log_prior_pdf += log(gsl_ran_gaussian_pdf((x(j_)-1.0), 1.0) ); 
	// Sunspot component
	j_ ++; 

	// Delta_{2,2}: Can be only identified through staying in ZLB regime.  The first three is not useful because there is only one endogenous error for this case.
	
	log_prior_pdf += log(gsl_ran_gaussian_pdf((x(j_)-1.0), 1.0) ); 
	// Sunspot component
	j_ ++; 

	log_prior_pdf += log(gsl_ran_gaussian_pdf((x(j_)-1.0), 1.0) ); 
	// Sunspot component
	j_ ++; 

	log_prior_pdf += log(gsl_ran_gaussian_pdf((x(j_)-1.0), 1.0) ); 
	// Sunspot component
	j_ ++; 

	// p1 and p2: p1 can be only identified through transitions from normal to ZLB; but p2 can NOT be identified because we don't have data from ZLB to normal.
	log_prior_pdf += log(gsl_ran_beta_pdf(x(j_), 2.0, 2.0) ); 
	// Interval (0.14, 0.86), probability of staying in regime 1
	j_ ++; 

	// double check
	if (j_ != x.dim)
	{
		cerr << "*******PriorDistrFunction_test_all::log_pdf(): dimension inconsistent*******\n"; 
		return CMSSM::MINUS_INFINITY_LOCAL; 
	}
	else 
		return log_prior_pdf; 
}

