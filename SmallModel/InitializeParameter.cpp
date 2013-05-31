#include "dw_dense_matrix.hpp"

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
