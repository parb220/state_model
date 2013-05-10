#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"

using namespace std; 

int CMSSM::LogLikelihood(double &log_likelihood, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob)
// Returns
// 	MODEL_NOT_PROPERLY_SET: if state_equation_function, measurement equation function or transition_prob_function not properly set
// 	SUCCESS:	
// 	ERROR_OCCURRED:
{
	if (CheckModelFunctions() != SUCCESS)
	{
		log_likelihood = -1.0e30; 
		return MODEL_NOT_PROPERLY_SET; 
	}
	if (UpdateStateModelParameters(0, y, x) == SUCCESS)
	{
		vector<vector<TDenseVector> > z_tm1;
                vector<vector<TDenseMatrix> > P_tm1;
                vector<TDenseVector> p_tm1;
		
		// Get initial values of z_0 and P_0 from the first few observations
		size_t initial_period = 4; 
		vector<TDenseVector> sub_y(y.begin(), y.begin()+initial_period);
		double initial_log_likelihood; 
		int kalman_error = KalmanFilter(initial_log_likelihood, z_tm1, P_tm1, p_tm1, sub_y, z_0, P_0, initial_prob, x); 
		if (kalman_error == SUCCESS)
		{
			vector<TDenseVector> new_z_0 = z_tm1.back();
			vector<TDenseMatrix> new_P_0 = P_tm1.back(); 

			kalman_error = KalmanFilter(log_likelihood, z_tm1, P_tm1, p_tm1, y, new_z_0, new_P_0, initial_prob, x);
			if (kalman_error == SUCCESS)
				return SUCCESS; 
		}
	}
	log_likelihood = -1.0e30; 
	return 	ERROR_OCCURRED; 
}
