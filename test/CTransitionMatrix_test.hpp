#ifndef _TRANSITION_MATRIX_TEST
#define _TRANSITION_MATRIX_TEST

#include "CMSSM.hpp"
#include "findequilibrium.hpp"
using namespace std; 

class CTransitionProbMatrixFunction_Test : public TransitionProbMatrixFunction
{
public:
	virtual void convert(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, size_t nS, size_t nTL, const TDenseVector &free_parameter, const TDenseVector &x);
};

#endif
