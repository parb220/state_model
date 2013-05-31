#ifndef _MEASUREMENT_EQUATION_FUNCTION_STSTM1
#define _MEASUREMENT_EQUATION_FUNCTION_STSTM1

#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"

const unsigned int GPISTAR_MEASUREMENT = 0;

class MeasurementEquationFunction_ststm1 : public MeasurementEquationFunction
{
protected: 
	virtual void ConvertXtoParameter(const TDenseVector &x) {}
public:
	virtual int convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &input_x, size_t nZ, size_t nY, size_t nU, size_t nNu);
	MeasurementEquationFunction_ststm1() : MeasurementEquationFunction() {}
	MeasurementEquationFunction_ststm1(const TDenseVector &p) : MeasurementEquationFunction(p) {}
	virtual ~MeasurementEquationFunction_ststm1() {}
};

int MeasurementEquationFunction_ststm1::convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &input_x, size_t nZ, size_t nY, size_t nU, size_t nNu)
{
	a = vector<TDenseVector>(nNu); 
	H = vector<TDenseMatrix>(nNu); 
	Phi_u = vector<TDenseMatrix>(nNu); 
	R = vector<TDenseMatrix>(nNu); 
	for (unsigned int i=0; i<nNu; i++)
	{
		a[i].Zeros(nY); 
		a[i].SetElement(fixed_parameter[GPISTAR_MEASUREMENT],1);	// fixed_parameter[0]: gpistar

		H[i].Zeros(nY,nZ); 
		H[i].SetElement(1.0,0,0); 
		H[i].SetElement(1.0,1,1); 
		H[i].SetElement(1.0,2,2); 

		Phi_u[i].Zeros(nU,nU); 
		R[i].Zeros(nU,nU); 
	}
	return SUCCESS;	// no error
}

#endif
