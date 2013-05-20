#ifndef _MEASUREMENT_EQUATION_FUNCTION_STSTM1
#define _MEASUREMENT_EQUATION_FUNCTION_STSTM1

#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"

const unsigned int GPISTAR_MEASUREMENT = 0;

class MeasurementEquationFunction_ststm1 : public MeasurementEquationFunction
{
public:
	virtual int convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &input_x);
	MeasurementEquationFunction_ststm1() : MeasurementEquationFunction() {}
	MeasurementEquationFunction_ststm1(const TDenseVector &p) : MeasurementEquationFunction(p) {}
	virtual ~MeasurementEquationFunction_ststm1() {}
};

int MeasurementEquationFunction_ststm1::convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &input_x)
{
	for (unsigned int i=0; i<a.size(); i++)
	{
		a[i].Zeros(3); 
		a[i].SetElement(fixed_parameter[GPISTAR_MEASUREMENT],1);	// fixed_parameter[0]: gpistar

		H[i].Zeros(3,7); 
		H[i].SetElement(1.0,0,0); 
		H[i].SetElement(1.0,1,1); 
		H[i].SetElement(1.0,2,2); 

		Phi_u[i].Zeros(0,0); 
		R[i].Zeros(0,0); 
	}
	return SUCCESS;	// no error
}

#endif
