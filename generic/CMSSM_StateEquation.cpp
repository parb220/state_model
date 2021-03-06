#include <cmath>
#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"

using namespace std; 

void CMSSM:: ClearStateEquation()
{
	b=vector<TDenseVector>(0); 
	F=vector<TDenseMatrix>(0); 
	Phi_e=vector<TDenseMatrix>(0); 
	V=vector<TDenseMatrix>(0); 
}

int CMSSM::UpdateStateEquationParameter(unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x)
// Return:
// 	-1: rational_expectation_function or state_equation_function is not properly st
// 	0: success
// 	>0: error code returned by rational_expectation_function or state_equation_function
{
	if (current_x.dim != x.dim || !(current_x == x)) 
	{
		if (!rational_expectation_function || !state_equation_function)
		{
			// cerr << "RationalExpectationFunction or StateEquationFunction is not properly set up.\n"; 
			ClearStateEquation(); 
			return MODEL_NOT_PROPERLY_SET;
		}
		int error_code; 
		error_code = rational_expectation_function->convert(A,B,Psi,Pi,C,x, nZ,nY, nU, nE, nExpectation); 
		if (error_code != SUCCESS)
		{
			// cerr << "Error occurred during RationalExpectationFunction call: " << error_code << endl; 
			ClearStateEquation(); 
			return error_code; 
		}
		
		error_code = state_equation_function->convert(b,F,Phi_e,V,A,B,Psi,Pi,C,x,nZ,nE,nExpectation,nNu); 
		if (error_code != SUCCESS)
		{
			// cerr << "Error occurred during StateEquationFunction call: " << error_code << endl; 
			ClearStateEquation(); 
			return error_code;  
		}
	}
	return SUCCESS; 
}
