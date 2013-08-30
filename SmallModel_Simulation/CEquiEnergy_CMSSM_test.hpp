#ifndef _EQUI_ENERGY_CMSSM_
#define _EQUI_ENERGY_CMSSM_

#include "CEquiEnergyModel.h"
#include "CMSSM_test.hpp"

class TDenseVector; 
class TDenseMatrix; 
class CStorageHead; 
class CEESParameter; 

class CEquiEnergy_CMSSM_test : public CEquiEnergyModel
{
public:
	CMSSM_test *target_model; 
	vector<TDenseVector> y; 
	vector<TDenseVector> z_0;
	vector<TDenseMatrix> P_0; 
	TDenseVector initial_prob;
  
	virtual double log_posterior_function(CSampleIDWeight &x); 
	// x cannot be constant, because x.weight will be set as the real log_posterior calculated
	// from target_model, where the returning value is the bounded log_posterior.
	// that is, return value = -(-x.weight>h_bound ? -x.weight:h_bound)/t_bound
	virtual double log_likelihood_function(const CSampleIDWeight &x); 
	// returning value is the real log_likelihood calculated from target_model
	
	CEquiEnergy_CMSSM_test(); 
	CEquiEnergy_CMSSM_test(bool _if_bounded, int eL, double _t, const CSampleIDWeight &_x, CMetropolis *_metropolis, time_t _time, CMSSM_test *_model); 
	virtual ~CEquiEnergy_CMSSM_test(); 

	virtual void SaveSampleToStorage(CStorageHead &storage, const CSampleIDWeight &sample);
	virtual void Take_Sample_Just_Drawn_From_Storage(const CSampleIDWeight &x_complete); 
};

#endif 
