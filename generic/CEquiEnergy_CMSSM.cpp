#include <vector>
#include "dw_dense_matrix.hpp"
#include "CEquiEnergy_CMSSM.hpp"
#include "InitializeParameter_test.hpp"

using namespace std; 

double CEquiEnergy_CMSSM :: log_posterior_function(CSampleIDWeight &x)
{
	if (!x.calculated)
	{
		target_model->LogPosterior(x.weight, x.data, y, z_0, P_0, initial_prob); 
		x.calculated = true; 
	}
	double bounded_log_posterior; 
	if (if_bounded)
		bounded_log_posterior = -(-x.weight>h_bound ? -x.weight:h_bound)/t_bound; 
	else 
		bounded_log_posterior = x.weight; 
	return bounded_log_posterior; 
}

double CEquiEnergy_CMSSM :: log_likelihood_function(const CSampleIDWeight &x)
{
	double log_likelihood; 
	vector<TDenseVector> z_tm1_last; 
	vector<TDenseMatrix> P_tm1_last; 
	TDenseVector p_tm1_last; 
	target_model->LogLikelihood(log_likelihood, z_tm1_last, P_tm1_last, p_tm1_last, x.data, y, z_0, P_0, initial_prob); 
	return log_likelihood; 	
}

CEquiEnergy_CMSSM :: CEquiEnergy_CMSSM() : CEquiEnergyModel(), target_model(NULL), y(vector<TDenseVector>(0)), z_0(vector<TDenseVector>(0)), P_0(vector<TDenseMatrix>(0)), initial_prob(TDenseVector(0)) 
{ }

CEquiEnergy_CMSSM :: CEquiEnergy_CMSSM(bool _if_bounded, unsigned int eL, double _h, double _t, const CSampleIDWeight &_x, CMetropolis *_metropolis, time_t _time, CMSSM *_model) : CEquiEnergyModel(_if_bounded, eL, _h, _t, _x, _metropolis, _time), target_model(_model), y(vector<TDenseVector>(0)), z_0(vector<TDenseVector>(0)), P_0(vector<TDenseMatrix>(0)), initial_prob(TDenseVector(0))
{}

CEquiEnergy_CMSSM :: ~CEquiEnergy_CMSSM()
{}
