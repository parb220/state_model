#include "CMSSM.hpp"
#include "CMSSM_Error_Code.hpp"
#include "MakeLbUb_test.hpp"

//
// Set lower and upper bounds for parameter vector x
//
// Returns
// 	SUCCESS upon success (no error)
// 	ERROR_OCCURRED if failed (number of bounds set is not equal to the number of parameters)
//



int MakeLbUb_test(TDenseVector &lb, TDenseVector &ub, size_t bDim, int choice)
{
	const double INFINITY = -CMSSM::MINUS_INFINITY_LOCAL; 
	size_t n_x = 15; 
	unsigned int i_ = 0; 
	unsigned int j_ = n_x; 

	lb.Zeros(bDim); 
	ub.Zeros(bDim); 
	if (choice == 1)
	{
		lb(i_) = 0.0025;	ub(i_) = 0.01;	i_=i_+1;	// Steady state quarterly rate of inflation (decimal) (e.g. 0.5/100 (0.005) quarterly corresponding to log(1.02) or 2% annually)
		lb(i_) = 0.005;	ub(i_) = 0.015;	i_=i_+1;	 // Steady state quarterly nominal interest rate, corresponding to i in Zha's notes (demical) (e.g., 0.015 quarterly corresponds to 6% annually).
		lb(i_) = 0.0;	ub(i_) = 1.0;	i_=i_+1;	// gthetaf (IS coefficient of expected output)
		lb(i_) = 0.0;	ub(i_) = 1.0;	i_=i_+1;	// gtau (IS coefficient of real interest rate)
		lb(i_) = 0.0;	ub(i_) = 1.0;	i_=i_+1;	// gomegaf (AS coefficient of expected inflation)
		lb(i_) = 0.0;	ub(i_) = 1.0;	i_=i_+1;	// gkappac (slope of Phillips curve for current output)
		lb(i_) = 1.0;	ub(i_) = 4.0;	i_=i_+1;	// gphipi (policy response to inflation)
		lb(i_) = 0.0;	ub(i_) = 1.0;	i_=i_+1;	// gphiy (plicy response to output gap)
		lb(i_) = 0.0;	lb(j_) = 0.0;	ub(i_) = 0.9999;	ub(j_) = 0.9999;	i_=i_+1;	j_=j_+1;   // grhomu (persistence to markup)
		lb(i_) = 0.0;   lb(j_) = 0.0;   ub(i_) = 0.9999;	ub(j_) = 0.9999;   	i_=i_+1;	j_=j_+1;   // grhorn (persistence to natural rate of interest)
		lb(i_) = 0.0;	lb(j_) = 0.0;	ub(i_) = 0.9999;	ub(j_) = 0.9999;	i_=i_+1;	j_=j_+1;   // grhoi (persistence to monetary policy)
		lb(i_) = 0.0;	lb(j_) = 0.0;	ub(i_) = 1.0;	ub(j_) = 1.0; 	i_=i_+1;	j_=j_+1;   	// gsigmamu (s.d. of markup shock)
		lb(i_) = 0.0;   lb(j_) = 0.0;   ub(i_) = 1.0;   ub(j_) = 1.0;   i_=i_+1;	j_=j_+1;   	// gsigmarn (s.d. of demand shock)
		lb(i_) = 0.0;   lb(j_) = 0.0;   ub(i_) = 1.0;   ub(j_) = 25.0;  i_=i_+1;	j_=j_+1;   	// x(j_), NOT gsigmai (s.d. of policy shock).  gsigmai = x(j_)/scl4gsigmai for 2nd regime while ub(j_) deals with x(j_), NOT gsigmai. 25.0 corresponds to 10 basis points annually.
		lb(i_) = -5.0;	ub(i_) = 5.0;			// giota (level adjustment for ZLB drop in the interest rate)
	
		lb(21) = -INFINITY;	ub(21) = INFINITY;	// Sunspot component
   		lb(22) = -INFINITY; 	ub(22) = INFINITY;      // Sunspot component
   		lb(23) = -INFINITY; 	ub(23) = INFINITY;     	// Sunspot component
   		lb(24) = -INFINITY; 	ub(24) = INFINITY; 	// Sunspot component
   		lb(25) = -INFINITY; 	ub(25) = INFINITY; 	// Sunspot component
   		lb(26) = -INFINITY; 	ub(26) = INFINITY; 	// Sunspot component
   		lb(27) = -INFINITY; 	ub(27) = INFINITY;	// Sunspot component
   		lb(28) = -INFINITY; 	ub(28) = INFINITY;	// Sunspot component
		lb(29) = 0.0; 	ub(29) = 1.0;			// probability of staying in regime 1   
   		lb(30) =  0.0;	ub(30) = 1.0;			// probability of staying in regime 2
	} 
	return SUCCESS; 
}

int MakeLbUb_test(TDenseVector &lb, TDenseVector &ub, const vector<int> &locs, size_t bDim, int choice)
{
	TDenseVector lb_all, ub_all; 
	int return_code; 
	if ((return_code=MakeLbUb_test(lb_all, ub_all, bDim, choice)) != SUCCESS)
		return 	return_code; 
	lb = lb_all.SubVector(locs); 
	ub = ub_all.SubVector(locs); 
	if (lb.dim == (int)locs.size() && ub.dim == (int)locs.size())
		return SUCCESS; 
	else 
		return ERROR_OCCURRED; 
}
