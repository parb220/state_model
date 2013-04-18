#ifndef _MEASUREMENT_EQUATION_FUNCTION_TEST
#define _MEASUREMENT_EQUATION_FUNCTION_TEST

#include "CMSSM.hpp"
#include "findequilibrium.hpp"

class MeasurementEquationFunction_test : public MeasurementEquationFunction
{
public:
	virtual bool convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &free_parameter, const TDenseVector &input_x);
};

bool MeasurementEquationFunction_test::convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &free_parameter, const TDenseVector &input_x)
{
	for (unsigned int i=0; i<a.size(); i++)
	{
		a[i].Zeros(4); 
		a[i].SetElement(free_parameter[GPISTAR_MEASUREMENT],1);	// free_parameter[0]: gpistar

		H[i].Zeros(3,7); 
		H[i].SetElement(1.0,0,0); 
		H[i].SetElement(1.0,1,1); 
		H[i].SetElement(1.0,2,2); 

		Phi_u[i].Zeros(0,0); 
		R[i].Zeros(0,0); 
	}
	return false;	// no error
}

#endif
