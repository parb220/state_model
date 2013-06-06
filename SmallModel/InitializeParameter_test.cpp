#include "InitializeParameter_test.hpp"
#include "CMSSM_Error_Code.hpp"
#include "dw_rand.h"

using namespace std; 

TDenseVector InitializeParameter(size_t n, const TDenseVector &fixed_parameter)
{
	TDenseVector gpistar(1, 0.005);	 // Steady state quarterly rate of inflation (decimal) (e.g. 0.5/100 (0.005) quarterly corresponding to log(1.02) or 2% annually)
	TDenseVector Rstar(1, 0.01);	// Steady state quarterly nominal interest rate, corresponding to i in Zha's notes (demical) (e.g., 0.015 quarterly corresponds to 6% annually).
	TDenseVector gthetaf(1, 0.6);  	// 0.9 gthetaf (IS coefficient of expected output)
	TDenseVector gtau(1, 0.15); 	// gtau (IS coefficient of real interest rate)
	TDenseVector gomegaf(1, 0.4);  	// gomegaf (AS coefficient of expected inflation)
	TDenseVector gkappac(1, 0.05);  // gkappac (slope of Phillips curve for current output)
	TDenseVector gphipi(2); gphipi(0) = 2.5; gphipi(1) = 0.0; 
	TDenseVector gphiy(2); gphiy(0) = 0.5; gphiy(1) = 0.0;
	TDenseVector grhomu(2); grhomu(0) = 0.85; grhomu(1) = 0.7;
	TDenseVector grhorn(2); grhorn(0) = 0.83; grhorn(1) = 0.7;	// Natural rate of interest.
	TDenseVector grhoi(2); grhoi(0) = 0.84; grhoi(1) = 0.7;		// Monetary policy.
	TDenseVector gsigmamu(2); gsigmamu(0) = 0.1; gsigmamu(1) = 0.1; // Markup shock.
	TDenseVector gsigmarn(2); gsigmarn(0) = 0.004; gsigmarn(1) = 0.004;	// Natural rate of interest.
	TDenseVector gsigmai(2); gsigmai(0) = 0.004; gsigmai(1) = 2.500e-05;   // Monetary policy. %Monetary policy.  1) 1 (or 2.5) basis point annually = 2.5000e-05 (or 6.2500e-05) quarterly.  2) 2.5000e-04 quarterly = 10 basis points annually in the ZLB.
	TDenseVector giota(1, 0.0);	// ZLB adjustment for level of fall in natural rate of interest.
	TDenseVector p1(1, 0.5);  
	TDenseVector p2(1, 0.5);

	TDenseVector x = RandomNormalVector(n); 
	size_t n_x = 15; 
	unsigned int i_ = 0, j_ = n_x; 	

	x(i_) = gpistar(0);	i_=i_+1;	// Steady state quarterly rate of inflation (decimal) (e.g. 0.5/100 (0.005) quarterly corresponding to log(1.02) or 2% annually)                
	x(i_) = Rstar(0);	i_=i_+1;	// Steady state quarterly nominal interest rate, corresponding to i in Zha's notes (demical) (e.g., 0.015 quarterly corresponds to 6% annually).
	x(i_) = gthetaf(0);  	i_=i_+1;        
	x(i_) = gtau(0);	i_=i_+1;        
	x(i_) = gomegaf(0);	i_=i_+1;        
	x(i_) = gkappac(0);	i_=i_+1;        
	x(i_) = gphipi(0);	i_=i_+1;        
	x(i_) = gphiy(0);	i_=i_+1;        
	x(i_) = grhomu(0);	x(j_) = grhomu(1);	i_=i_+1;	j_=j_+1;
	x(i_) = grhorn(0);	x(j_) = grhorn(1);	i_=i_+1;	j_=j_+1;    // Natural rate of interest.
	x(i_) = grhoi(0);  	x(j_) = grhoi(1);	i_=i_+1;	j_=j_+1;    // Monetary policy.
	x(i_) = gsigmamu(0);	x(j_)  = gsigmamu(1);	i_=i_+1;  j_=j_+1;
	x(i_) = gsigmarn(0);	x(j_)  = gsigmarn(1);	i_=i_+1;  j_=j_+1;    // Natural rate of interest.
	x(i_) = gsigmai(0);	x(j_)  = gsigmai (1)*fixed_parameter[0]; 	i_=i_+1;  j_=j_+1;    // Monetary policy.
	x(i_) = giota(0);                                                                 // ZLB adjustment for level of fall in natural rate of interest.
	x(n-2) = p1(0);
	x(n-1) = p2(0);

	return x; 
}

int RandomInit_test_1st(TDenseVector &x, const vector<int> &pos)
{
	if (x.dim < pos.size())
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
	unsigned int count = 0; 
	if (pos[count] >= x.dim)
		return ERROR_OCCURRED; 
	x(pos[count]) = dw_gaussian_rnd();	// giota (level adjustment for ZLB drop in the interest rate)
	count ++; 
	if (pos[count] >= x.dim)
		return ERROR_OCCURRED; 
	x(pos[count]) = dw_uniform_rnd(); 	// [0 1]; grhomu (persistence to markup)

	count ++; 
	if (pos[count] >= x.dim)
                return ERROR_OCCURRED;
	x(pos[count]) = dw_uniform_rnd(); 	// [0 1]; grhorn (persistence to natural rate of interest)

	count ++; 
	if (pos[count] >= x.dim)
		return ERROR_OCCURRED; 
	x(pos[count]) = dw_uniform_rnd(); 	// [0 1]; grhoi (persistence to monetary policy)

	count ++; 
	if (pos[count] >= x.dim)
                return ERROR_OCCURRED;
	x(pos[count]) = dw_uniform_rnd()/100; 	// [0.0 0.01]; gsigmamu (s.d. of markup shock)

	count ++; 
	if (pos[count] >= x.dim)
		return 	ERROR_OCCURRED; 
	x(pos[count]) = dw_uniform_rnd()/100;	// [0.0 0.01]; gsigmarn (s.d. of demand shock)

	count ++; 
	if (pos[count] >= x.dim)
		return ERROR_OCCURRED; 
	x(pos[count]) = 1.0+dw_uniform_rnd()*1.5; 	// [1.0 2.5]; gsigmai (scl4gsigmai * s.d. of policy shock)
	count ++; 
	if (pos[count] >= x.dim)
		return ERROR_OCCURRED; 
	x(pos[count]) = dw_gaussian_rnd(); 	// Interval [-1, 1]; Sunspot component

	count ++;
        if (pos[count] >= x.dim)
                return ERROR_OCCURRED;
        x(pos[count]) = dw_gaussian_rnd();      // Interval [-1, 1]; Sunspot component
	
	count ++;
        if (pos[count] >= x.dim)
                return ERROR_OCCURRED;
        x(pos[count]) = dw_gaussian_rnd();      // Interval [-1, 1]; Sunspot component

	count ++; 
	if (count != pos.size())
		return ERROR_OCCURRED; 

	return SUCCESS; 
}
