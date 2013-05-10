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
		log_posterior = -1.0e30; 
		return MODEL_NOT_PROPERLY_SET; 
	}
	double log_likelihood; 
	int likelihood_error = LogLikelihood(log_likelihood, x, y, z_0, P_0, initial_prob); 
	if (likelihood_error == SUCCESS)
	{
		log_posterior = likelihood_error + prior_distr_function->log_pdf(x,prior_distr_parameter); 
		return SUCCESS; 
	}	
	else 
	{
		log_posterior = -1.0e30; 
		return ERROR_OCCURRED; 
	}
}
