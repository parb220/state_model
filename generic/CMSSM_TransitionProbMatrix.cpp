#include "CMSSM.hpp"

using namespace std; 

TDenseMatrix CMSSM::GetTranstionProbMatrix(unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x) const
{
	TDenseMatrix Q; 
	if (transition_prob_function) 
		transition_prob_function->convert(Q, t, y, nS, nTL, transition_prob_parameter, x);
	return Q;  
}
