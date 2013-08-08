#ifndef _EQUI_ENERGY_CMSSM_
#define _EQUI_ENERGY_CMSSM_

#include "CEquiEnergyModel.h"
#include "CMSSM.hpp"

class TDenseVector; 
class TDenseMatrix; 
class CStorageHead; 
class CEESParameter; 

class CEquiEnergy_CMSSM : public CEquiEnergyModel
{
public:
	CMSSM *target_model; 
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
	
	CEquiEnergy_CMSSM(); 
	CEquiEnergy_CMSSM(bool _if_bounded, unsigned int eL, double _h, double _t, const CSampleIDWeight &_x, CMetropolis *_metropolis, time_t _time, CMSSM *_model); 
	virtual ~CEquiEnergy_CMSSM(); 
};

#endif 
