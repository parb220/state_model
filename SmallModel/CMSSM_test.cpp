#include "CMSSM_test.hpp"
#include "CMSSM_Error_Code.hpp"

using namespace std; 

int CMSSM_test::SetConstantPart(const TDenseVector &_x_constant)
{
	if (_x_constant.dim != (int)location_constant.size())
		return ERROR_OCCURRED; 
	for (unsigned int i=0; i<location_constant.size(); i++)
		x_variable_constant[location_constant[i]] =  _x_constant[i]; 
	return SUCCESS; 
}

int CMSSM_test::LogLikelihood(double &log_likelihood, vector<TDenseVector> &z_tm1_last, vector<TDenseMatrix> &P_tm1_last, TDenseVector &p_tm1_last, const TDenseVector &x_variable, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob)
{
	if (x_variable.dim != (int)location_variable.size())
	{
		log_likelihood = MINUS_INFINITY_LOCAL; 
                return ERROR_OCCURRED;
	}
	for (unsigned int i=0; i<location_variable.size(); i++)
		x_variable_constant[location_variable[i]] = x_variable[i]; 
	
	return CMSSM::LogLikelihood(log_likelihood, z_tm1_last, P_tm1_last, p_tm1_last, x_variable_constant, y, z_0, P_0, initial_prob); 
}

int CMSSM_test::LogPosterior(double &log_posterior, const TDenseVector &x_variable, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob)
{
	if (x_variable.dim != (int)location_variable.size()) 
	{
		log_posterior = MINUS_INFINITY_LOCAL;
		return ERROR_OCCURRED; 
	}
	return CMSSM::LogPosterior(log_posterior, x_variable, y, z_0, P_0, initial_prob); 
}
