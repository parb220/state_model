#ifndef _C_STATE_EQUATION_FUNCTION_STSTM1_
#define _C_STATE_EQUATION_FUNCTION_STSTM1_

#include "CMSSM.hpp"

class StateEquationFunction_ststm1 : public StateEquationFunction
{
protected: 
	virtual void ConvertXtoParameter(const TDenseVector &x) {}
public:
	virtual int convert(vector<TDenseVector> &b, vector<TDenseMatrix> &F, vector<TDenseMatrix> &Phi_e, vector<TDenseMatrix> &V, const vector<vector<TDenseMatrix> > &A, const vector<vector<TDenseMatrix> > &B, const vector<vector<TDenseMatrix> > &Psi, const vector<vector<TDenseMatrix> >&Pi, const vector<vector<TDenseVector> > &C, const TDenseVector &x, size_t nZ, size_t nE, size_t nExpectation, size_t nNu); 
	StateEquationFunction_ststm1() : StateEquationFunction() {}
	StateEquationFunction_ststm1(const TDenseVector &p): StateEquationFunction(p) {}
	virtual ~StateEquationFunction_ststm1() {}
}; 

#endif
