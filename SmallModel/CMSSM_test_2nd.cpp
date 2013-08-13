#include "CMSSM_test_2nd.hpp"
#include "CMSSM_Error_Code.hpp"

using namespace std; 

int CMSSM_test_2nd::LogLikelihood(double &log_likelihood, vector<TDenseVector> &z_tm1_last, vector<TDenseMatrix> &P_tm1_last, TDenseVector &p_tm1_last, const TDenseVector &x_variable, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob)
{
	if (x_variable.dim != (int)location_variable.size())
	{
		log_likelihood = MINUS_INFINITY_LOCAL; 
		return ERROR_OCCURRED; 
	}
	for (unsigned int i=0; i<location_variable.size(); i++)
		x_variable_constant[location_variable[i]] = x_variable[i]; 

	if (CheckModelFunctions() != SUCCESS)
        {
                log_likelihood = MINUS_INFINITY_LOCAL;
                return MODEL_NOT_PROPERLY_SET;
        }
        if (UpdateStateModelParameters(0, y, x_variable_constant) == SUCCESS)
        {
		vector<vector<TDenseVector> > z_tm1;
                vector<vector<TDenseMatrix> > P_tm1;
                vector<TDenseVector> p_tm1;
	
		int kalman_error = KalmanFilter(log_likelihood, z_tm1, P_tm1, p_tm1, y, z_0, P_0, initial_prob, x_variable_constant);
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

