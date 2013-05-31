#include "CMSSM_test_1st.hpp"
#include "CMSSM_Error_Code.hpp"

using namespace std; 

int CMSSM_test_1st::LogPosterior(double &log_posterior, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob)
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
		TDenseVector x1(location.size()); 
		for (unsigned int i=0; i<location.size(); i++)
			x1[i] = x[location[i]]; 
		double log_prior = prior_distr_function->log_pdf(x1); 
		log_posterior = log_likelihood + log_prior; 
		return SUCCESS; 
	}	
	else 
	{
		log_posterior = MINUS_INFINITY_LOCAL; 
		return ERROR_OCCURRED; 
	}
}
