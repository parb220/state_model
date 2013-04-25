#ifndef _TRANSITION_MATRIX_STSTM1
#define _TRANSITION_MATRIX_STSTM1

#include "CMSSM.hpp"
#include "findequilibrium.hpp"
using namespace std; 

class TransitionProbMatrixFunction_ststm1 : public TransitionProbMatrixFunction
{
public:
	virtual int convert(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, size_t nS, size_t nTL, const TDenseVector &fixed_parameter, const TDenseVector &x);
};

#endif
