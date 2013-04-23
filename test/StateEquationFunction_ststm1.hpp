#ifndef _C_STATE_EQUATION_FUNCTION_STSTM1_
#define _C_STATE_EQUATION_FUNCTION_STSTM1_

#include "CMSSM.hpp"
#include "findequilibrium.hpp"

class StateEquationFunction_ststm1 : public StateEquationFunction
{
	virtual int convert(vector<TDenseVector> &b, vector<TDenseMatrix> &F, vector<TDenseMatrix> &Phi_e, vector<TDenseMatrix> &V, const vector<vector<TDenseMatrix> > &A, const vector<vector<TDenseMatrix> > &B, const vector<vector<TDenseMatrix> > &Psi, const vector<vector<TDenseMatrix> >&Pi, const vector<vector<TDenseVector> > &C, const TDenseVector &fixed_parameter, const TDenseVector &x, size_t, size_t, size_t); 
}; 

#endif
