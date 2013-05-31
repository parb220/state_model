#ifndef _MAKE_A_B_PSI_PI_C_STSTM1
#define _MAKE_A_B_PSI_PI_C_STSTM1

#include "CMSSM.hpp"

const unsigned int GLAMBDAZ_RE=0;
const unsigned int RN_RE = 1;
const unsigned int GPISTAR_RE = 2;
const unsigned int SC_RE = 3;
const unsigned int RLOW_RE = 4;

class RationalExpectationFunction_ststm1 : public RationalExpectationFunction
{
protected: 
	virtual void ConvertXtoParameter(const TDenseVector &x) {}; 
public:
        virtual int convert(vector<vector<TDenseMatrix> > &A, vector<vector<TDenseMatrix> > &B, vector<vector<TDenseMatrix> > &Psi, vector<vector<TDenseMatrix> >&Pi, vector<vector<TDenseVector> > &C, const TDenseVector &input_x, size_t nZ, size_t nY, size_t nU, size_t nE, size_t nExpectation);
	RationalExpectationFunction_ststm1() : RationalExpectationFunction() {}
	RationalExpectationFunction_ststm1(const TDenseVector &p) : RationalExpectationFunction(p) {}
	virtual ~RationalExpectationFunction_ststm1() {}	
};

#endif
