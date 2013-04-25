#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"

using namespace std; 

int CMSSM::GetTranstionProbMatrix(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x) const
{
	if (transition_prob_function == NULL)
		return MODEL_NOT_PROPERLY_SET; 
	else if (transition_prob_function->convert(Q, t, y, nS, nTL, transition_prob_parameter, x) == SUCCESS)
		return SUCCESS;   
	else 
		return ERROR_OCCURRED;  
}
