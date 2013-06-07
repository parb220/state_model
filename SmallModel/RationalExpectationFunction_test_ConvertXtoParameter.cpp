#include "RationalExpectationFunction_test.hpp"

////////////////////////////////////////////////////
// ftd_convert_x2meaningfulpars_const_hybrid.m
////////////////////////////////////////////////////

using namespace std;
void RationalExpectationFunction_test :: ConvertXtoParameter(const TDenseVector &x)
{
	scl4Sigmai2 = fixed_parameter(0); 
	ilow =  fixed_parameter(1);	// Net rate: 10 basis points (0.1%) interest rate annually at zero bound or 0.025% quarterly at zero bound.
	
	size_t n_x = 15; 
	unsigned int i_=0, j_=n_x; 

	DoubledPistar = x(i_);		i_++;	// Steady state quarterly rate of inflation (decimal) (e.g. 0.5/100 (0.005) quarterly corresponding to log(1.02) or 2% annually)
	iss = x(i_); 			i_++;		// Rstar, steady state quarterly nominal interest rate, corresponding to i in Zha's notes (demical) (e.g., 0.015 quarterly corresponds to 6% annually).
	Thetaf = x(i_); 		i_++;	// gthetaf (IS coefficient of expected output)
	Tau = x(i_); 			i_++;		// gtau (IS coefficient of real interest rate)
	Omegaf = x(i_); 		i_++;       // gomegaf (AS coefficient of expected inflation)
	Kappac = x(i_); 		i_++;       // gkappac (slope of Phillips curve for current output)
	
	PhiDoubledPi = x(i_); 		i_++; 
	Phiy = x(i_); 			i_++; 
	RhoMu1 = x(i_); 		i_++;	 	// Persistence parameter of Markup
	Rhorn1 = x(i_);			i_++; 		// Persistence parameter of Natural rate of interest
	Rhoi1 = x(i_); 			i_++; 		// Persistence parameter of interest rate (monetary policy)
	SigmaMu1 = x(i_); 		i_++;		// S.d. of markup
	Sigmarn1 = x(i_);  		i_++;	 	// S.d. of natural rate of interest
	Sigmai1 = x(i_); 		i_++;		// S.d. of interest rate (monetary policy) 
	Iota = x(i_); 					// ZLB adjustment for level of fall in natural rate of interest. 

	RhoMu2 = x(j_);		j_++; 		// Persistence parmaeter of markup.
	Rhorn2 = x(j_); 	j_++; 		// Persistence parmaeter of the natural rate of interest.
	Rhoi2 = x(j_);		j_++; 		// Persistence parmaeter of the interest rate (monetary policy).
	SigmaMu2 = x(j_);	j_++; 		// S.d. parmaeter of markup.
	Sigmarn2 = x(j_);	j_++; 		// S.d. parmaeter of the natural rate of interest.
	Sigmai2 = x(j_)/scl4Sigmai2;		// S.d. parmaeter of the interest rate (monetary policy).
						// Scale for s.d. of interest rate in ZLB. Scale for s.d. of interest rate in ZLB. 1) 1 (or 2.5) basis point annually = 2.5000e-05 (or 6.2500e-05) quarterly.  2) 2.5000e-04 quarterly = 10 basis points annually in the ZLB. 
						// The scale is needed because, numerically, the prior of inverse gamma cannot cover the range from 0.0002/4 (2 basis points quarterly to 0.01 (1%)
						// Note that this scale will be rescaled at the end when all estimation is completed.
}
