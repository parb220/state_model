#include "InitializeParameter_test.hpp"
#include "CMSSM_Error_Code.hpp"
#include "dw_rand.h"


////////////////////////////////////////////////////////////
//ftd_initialize_meaningfulpars_const_hybrid + ftd_convert_meaningfulpars2x_const_hybrid
////////////////////////////////////////////////////////////

using namespace std; 

TDenseVector InitializeParameter(size_t n, const TDenseVector &fixed_parameter)
{
	TDenseVector DoubledPistar(1, 0.005);	 // Steady state quarterly rate of inflation (decimal) (e.g. 0.5/100 (0.005) quarterly corresponding to log(1.02) or 2% annually)
	TDenseVector iss(1, 0.01);	// Steady state quarterly nominal interest rate, corresponding to i in Zha's notes (demical) (e.g., 0.015 quarterly corresponds to 6% annually).
	TDenseVector Thetaf(1, 0.6);  	// 0.9 gthetaf (IS coefficient of expected output)
	TDenseVector Tau(1, 0.15); 	// gtau (IS coefficient of real interest rate)
	TDenseVector Omegaf(1, 0.4);  	// gomegaf (AS coefficient of expected inflation)
	TDenseVector Kappac(1, 0.05);  // gkappac (slope of Phillips curve for current output)
	TDenseVector PhiDoubledPi(1, 2.5);  
	TDenseVector Phiy(1, 0.5);
	TDenseVector RhoMu1(1, 0.85), RhoMu2(1,0.7);  
	TDenseVector Rhorn1(1, 0.83), Rhorn2(1,0.7); 	// Natural rate of interest.
	TDenseVector Rhoi1(1,0.84), Rhoi2(1,0.7); 	// Monetary policy.
	TDenseVector SigmaMu1(1,0.1), SigmaMu2(1,0.1);  // Markup shock.
	TDenseVector Sigmarn1(1,0.004), Sigmarn2(1, 0.004);  // Natural rate of interest.
	TDenseVector Sigmai1(1,0.004), Sigmai2(1,2.500e-05);   // Monetary policy. %Monetary policy.  1) 1 (or 2.5) basis point annually = 2.5000e-05 (or 6.2500e-05) quarterly.  2) 2.5000e-04 quarterly = 10 basis points annually in the ZLB.
	TDenseVector Iota(1, 0.0);	// ZLB adjustment for level of fall in natural rate of interest.
	TDenseVector p1(1, 0.5);  
	TDenseVector p2(1, 0.5);

	TDenseVector x = RandomNormalVector(n); 
	size_t n_x = 15; 
	unsigned int i_ = 0, j_ = n_x; 	

	x(i_) = DoubledPistar(0);	i_++;	// Steady state quarterly rate of inflation (decimal) (e.g. 0.5/100 (0.005) quarterly corresponding to log(1.02) or 2% annually)                
	x(i_) = iss(0);	i_++;	// Steady state quarterly nominal interest rate, corresponding to i in Zha's notes (demical) (e.g., 0.015 quarterly corresponds to 6% annually).
	x(i_) = Thetaf(0);  	i_++;        
	x(i_) = Tau(0);	i_++;        
	x(i_) = Omegaf(0);	i_++;        
	x(i_) = Kappac(0);	i_++;        
	x(i_) = PhiDoubledPi(0);	i_++;        
	x(i_) = Phiy(0);	i_++;        
	x(i_) = RhoMu1(0);	x(j_) = RhoMu2(0);	i_++;	j_++;
	x(i_) = Rhorn1(0);	x(j_) = Rhorn2(0);	i_++;	j_++;    // Natural rate of interest.
	x(i_) = Rhoi1(0);  	x(j_) = Rhoi2(0);	i_++;	j_++;    // Monetary policy.
	x(i_) = SigmaMu1(0);	x(j_)  = SigmaMu2(0);	i_++;  j_++;
	x(i_) = Sigmarn1(0);	x(j_)  = Sigmarn2(0);	i_++;  j_++;    // Natural rate of interest.
	x(i_) = Sigmai1(0);	x(j_)  = Sigmai2(0)*fixed_parameter[0]; 	i_++;  j_++;    // Monetary policy.
	x(i_) = Iota(0);                                                                 // ZLB adjustment for level of fall in natural rate of interest.
	x(n-2) = p1(0);
	x(n-1) = p2(0);

	return x; 
}

int RandomInit_test_1st(TDenseVector &x, const vector<int> &pos)
{
	if (x.dim < (int)pos.size())
		return ERROR_OCCURRED; 
	
	x(0)  = 0.0025 + dw_uniform_rnd()/200;   	// Steady state quarterly rate of inflation (decimal) (e.g. 0.5/100 (0.005) quarterly corresponding to log(1.02) or 2% annually)
	x(1)  = x(0) + dw_uniform_rnd()/300;   	// Steady state quarterly nominal interest rate, corresponding to i in Zha's notes (demical) (e.g., 0.015 quarterly corresponds to 6% annually).
	x(2)  = dw_uniform_rnd();             	// gthetaf (IS coefficient on expected output)        
	x(3)  = dw_uniform_rnd();             	// gtau (IS coefficient on real interest rate)        
	x(4)  = dw_uniform_rnd();             	// gomegaf (AS coefficient on expected inflation)     
	x(5)  = dw_uniform_rnd();             	// gkappac (slope of Phillips curve on current output)
	x(6)  = 3*dw_uniform_rnd() + 1.3;     	// gphipi (policy response to inflation)   
	x(7)  = dw_uniform_rnd();             	// gphiy (plicy response to output gap)      
	x(8)  = dw_uniform_rnd()/5+0.75;      	// grhomu (persistence to markup)
	x(9)  = dw_uniform_rnd()/5+0.75;     	// grhorn (persistence to natural rate of interest)   
	x(10) = dw_uniform_rnd();             	// grhoi (persistence to monetary policy)      
	x(11) = dw_uniform_rnd()/100;         	// gsigmamu (s.d. of markup shock)      
	x(12) = dw_uniform_rnd()/100;         	// gsigmarn (s.d. of demand shock)      
	x(13) = dw_uniform_rnd()/100;         	// gsigmai (s.d. of policy shock) 
	return SUCCESS; 
}

int RandomInit_test_2nd(TDenseVector &x, const vector<int> &pos)
{
	if (x.dim < (int)pos.size())
		return ERROR_OCCURRED; 
	x(0) = dw_gaussian_rnd();	// giota (level adjustment for ZLB drop in the interest rate)
	x(1) = dw_uniform_rnd(); 	// [0 1]; grhomu (persistence to markup)
	x(2) = dw_uniform_rnd(); 	// [0 1]; grhorn (persistence to natural rate of interest)
	x(3) = dw_uniform_rnd(); 	// [0 1]; grhoi (persistence to monetary policy)
	x(4) = dw_uniform_rnd()/100; 	// [0.0 0.01]; gsigmamu (s.d. of markup shock)
	x(5) = dw_uniform_rnd()/100;	// [0.0 0.01]; gsigmarn (s.d. of demand shock)
	x(6) = 1.0+dw_uniform_rnd()*1.5; 	// [1.0 2.5]; gsigmai (scl4gsigmai * s.d. of policy shock)
	x(7) = dw_gaussian_rnd(); 	// Interval [-1, 1]; Sunspot component

        x(8) = dw_gaussian_rnd();      // Interval [-1, 1]; Sunspot component
	
        x(9) = dw_gaussian_rnd();      // Interval [-1, 1]; Sunspot component

	return SUCCESS; 
}
