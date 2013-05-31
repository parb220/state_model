#include "CMSSM_test_2nd.hpp"
#include "CMSSM_Error_Code.hpp"

using namespace std; 

int CMSSM_test_2nd::LogLikelihood(double &log_likelihood, vector<TDenseVector> &z_tm1_last, vector<TDenseMatrix> &P_tm1_last, TDenseVector &p_tm1_last, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob)
{
	if (CheckModelFunctions() != SUCCESS)
        {
                log_likelihood = MINUS_INFINITY_LOCAL;
                return MODEL_NOT_PROPERLY_SET;
        }
        if (UpdateStateModelParameters(0, y, x) == SUCCESS)
        {
		vector<vector<TDenseVector> > z_tm1;
                vector<vector<TDenseMatrix> > P_tm1;
                vector<TDenseVector> p_tm1;
	
		int kalman_error = KalmanFilter(log_likelihood, z_tm1, P_tm1, p_tm1, y, z_0, P_0, initial_prob, x);
		if (kalman_error == SUCCESS)
		{
			z_tm1_last = z_tm1.back(); 
			P_tm1_last = P_tm1.back(); 
			p_tm1_last = p_tm1.back(); 
			return SUCCESS; 
		}
	}
	log_likelihood = MINUS_INFINITY_LOCAL; 
	return ERROR_OCCURRED; 
}

int CMSSM_test_2nd::LogPosterior(double &log_posterior, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob)
{
	if (prior_distr_function == NULL || CheckModelFunctions() != SUCCESS)
        {
                log_posterior = CMSSM::MINUS_INFINITY_LOCAL;
                return MODEL_NOT_PROPERLY_SET;
        }
	
	double log_likelihood;
        vector<TDenseVector> z_tm1_last;
        vector<TDenseMatrix> P_tm1_last;
        TDenseVector p_tm1_last;
        int likelihood_error = LogLikelihood(log_likelihood, z_tm1_last, P_tm1_last, p_tm1_last, x, y, z_0, P_0, initial_prob);

	if (likelihood_error == SUCCESS)
        {
                TDenseVector x2(location.size());
                for (unsigned int i=0; i<location.size(); i++)
                        x2[i] = x[location[i]];
                double log_prior = prior_distr_function->log_pdf(x2);
                log_posterior = log_likelihood + log_prior;
                return SUCCESS;
        }
        else
        {
                log_posterior = CMSSM::MINUS_INFINITY_LOCAL;
                return ERROR_OCCURRED;
        }

}
