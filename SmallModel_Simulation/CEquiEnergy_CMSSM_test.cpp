#include <vector>
#include <ctime>
#include "dw_dense_matrix.hpp"
extern "C"
{
	#include "dw_rand.h"
}
#include "CEquiEnergy_CMSSM_test.hpp"
#include "InitializeParameter_test.hpp"

using namespace std; 

double CEquiEnergy_CMSSM_test :: log_posterior_function(CSampleIDWeight &x)
{
	if (!x.calculated)
	{
		target_model->LogPosterior(x.weight, x.data, y, z_0, P_0, initial_prob); 
		x.calculated = true; 
	}
	double bounded_log_posterior; 
	if (if_bounded)
		bounded_log_posterior = x.weight/t_bound; 
	else 
		bounded_log_posterior = x.weight; 
	return bounded_log_posterior; 
}

double CEquiEnergy_CMSSM_test :: log_likelihood_function(const CSampleIDWeight &x)
{
	double log_likelihood; 
	vector<TDenseVector> z_tm1_last; 
	vector<TDenseMatrix> P_tm1_last; 
	TDenseVector p_tm1_last; 
	target_model->LogLikelihood(log_likelihood, z_tm1_last, P_tm1_last, p_tm1_last, x.data, y, z_0, P_0, initial_prob); 
	return log_likelihood; 	
}

CEquiEnergy_CMSSM_test :: CEquiEnergy_CMSSM_test() : CEquiEnergyModel(), target_model(NULL), y(vector<TDenseVector>(0)), z_0(vector<TDenseVector>(0)), P_0(vector<TDenseMatrix>(0)), initial_prob(TDenseVector(0)) 
{ }

CEquiEnergy_CMSSM_test :: CEquiEnergy_CMSSM_test(bool _if_bounded, unsigned int eL, double _t, const CSampleIDWeight &_x, CMetropolis *_metropolis, time_t _time, CMSSM_test *_model) : CEquiEnergyModel(_if_bounded, eL, _t, _x, _metropolis, _time), target_model(_model), y(vector<TDenseVector>(0)), z_0(vector<TDenseVector>(0)), P_0(vector<TDenseMatrix>(0)), initial_prob(TDenseVector(0))
{}

CEquiEnergy_CMSSM_test :: ~CEquiEnergy_CMSSM_test()
{}

void CEquiEnergy_CMSSM_test::SaveSampleToStorage(CStorageHead &storage, const CSampleIDWeight &sample)
{
	vector<int> locs_variable = target_model->locs_variable(); 
	vector<int> locs_constant = target_model->locs_constant(); 
	CSampleIDWeight save_sample;
	save_sample.data.Resize(locs_variable.size()+locs_constant.size()); 
	save_sample.data.SetSubVector(locs_variable, sample.data); 
	save_sample.data.SetSubVector(locs_constant, target_model->GetConstantPart()); 
	save_sample.weight = sample.weight; 
	save_sample.id = sample.id; 
	storage.DepositSample(energy_level, storage.BinIndex(energy_level, -save_sample.weight),save_sample); 
}

void CEquiEnergy_CMSSM_test::Take_Sample_Just_Drawn_From_Storage(const CSampleIDWeight &x_complete)
{
	target_model->SetConstantPart(x_complete.data.SubVector(target_model->locs_constant())); 
	current_sample.data = x_complete.data.SubVector(target_model->locs_variable()); 
	current_sample.weight = x_complete.weight; 
	current_sample.id = (int)(time(NULL)-timer_when_started);
	current_sample.calculated = x_complete.calculated; 	
}
