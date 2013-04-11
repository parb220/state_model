#ifndef _MEASUREMENT_EQUATION_FUNCTION_TEST
#define _MEASUREMENT_EQUATION_FUNCTION_TEST

#include "CMSSM.hpp"
#include "findequilibrium.hpp"

class MeasurementEquationFunction_test : public MeasurementEquationFunction
{
public:
	virtual bool convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &x, const TDenseVector &free_parameter);
};

bool MeasurementEquationFunction_test::convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &x, const TDenseVector &free_parameter)
{
	for (unsigned int i=0; i<a.size(); i++)
	{
		a[i].Zeros(a[i].dim); 
		a[i].SetElement(free_parameter[GPISTAR_MEASUREMENT],1);	// free_parameter[0]: gpistar

		H[i].Zeros(H[i].rows,H[i].cols); 
		H[i].SetElement(1.0,0,0); 
		H[i].SetElement(1.0,1,1); 
		H[i].SetElement(1.0,2,2); 

		Phi_u[i].Zeros(Phi_u[i].rows,Phi_u[i].cols); 
	}
	return false;	// no error
}

#endif
