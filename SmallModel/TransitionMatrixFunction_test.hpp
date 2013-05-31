#ifndef _TRANSITION_MATRIX_STSTM1
#define _TRANSITION_MATRIX_STSTM1

#include "CMSSM.hpp"
using namespace std; 

const unsigned int CUTOFF_TRANSITION = 0;
class TransitionProbMatrixFunction_test : public TransitionProbMatrixFunction
{
protected:
	virtual void ConvertXtoParameter(const TDenseVector &x) {}
public:
	virtual int convert(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, size_t nS, size_t nTL, const TDenseVector &x);
	TransitionProbMatrixFunction_test() : TransitionProbMatrixFunction() {}
	TransitionProbMatrixFunction_test(const TDenseVector &p) : TransitionProbMatrixFunction(p) {}	// fixed_parameter[2]: cutoff
	virtual ~TransitionProbMatrixFunction_test() {}
};

#endif
