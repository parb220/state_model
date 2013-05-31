#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"

using namespace std; 

int CMSSM::LogPosterior(double &log_posterior, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob)
// Returns
// 	MODEL_NOT_PROPERLY_SET: if state_equation_function, measurement equation function or transition_prob_function, prior_distr_function not properly set
// 	SUCCESS:	
// 	ERROR_OCCURRED:
{
	if (prior_distr_function == NULL || CheckModelFunctions() != SUCCESS)
	{
		log_posterior = MINUS_INFINITY_LOCAL; 
		return MODEL_NOT_PROPERLY_SET; 
	}
	double log_likelihood; 
	vector<TDenseVector> z_tm1_last; 
	vector<TDenseMatrix> P_tm1_last; 
	TDenseVector p_tm1_last; 
	int likelihood_error = LogLikelihood(log_likelihood, z_tm1_last, P_tm1_last, p_tm1_last, x, y, z_0, P_0, initial_prob); 
	if (likelihood_error == SUCCESS)
	{
		log_posterior = log_likelihood + prior_distr_function->log_pdf(x); 
		return SUCCESS; 
	}	
	else 
	{
		log_posterior = MINUS_INFINITY_LOCAL; 
		return ERROR_OCCURRED; 
	}
}
