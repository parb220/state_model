#ifndef _TRANSITION_MATRIX_STSTM1
#define _TRANSITION_MATRIX_STSTM1

#include "CMSSM.hpp"
using namespace std; 

const unsigned int CUTOFF_TRANSITION = 0;
class TransitionProbMatrixFunction_ststm1 : public TransitionProbMatrixFunction
{
public:
	virtual int convert(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, size_t nS, size_t nTL, const TDenseVector &x);
	TransitionProbMatrixFunction_ststm1() : TransitionProbMatrixFunction() {}
	TransitionProbMatrixFunction_ststm1(const TDenseVector &p) : TransitionProbMatrixFunction(p) {}
	virtual ~TransitionProbMatrixFunction_ststm1() {}
};

#endif
