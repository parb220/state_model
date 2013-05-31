#ifndef _CMSSM_TEST_2ND_
#define _CMSSM_TEST_2ND_
#include "CMSSM.hpp"

class CMSSM_test_2nd : public CMSSM
{
protected:
	vector<int> location; 
public:
	CMSSM_test_2nd() : CMSSM(), location() {}
	CMSSM_test_2nd(const vector<int> &_loc) : CMSSM(), location(_loc) {}
	CMSSM_test_2nd(const vector<int> &_loc, size_t _nL, size_t _nTL, size_t _nS, size_t _nZ, size_t _nY, size_t _nU, size_t _nE, size_t _nExpectation, RationalExpectationFunction * ref=NULL, StateEquationFunction * sef=NULL, MeasurementEquationFunction * mef=NULL, TransitionProbMatrixFunction * tpf=NULL, PriorDistributionFunction * pdf= NULL, const TDenseVector & _x=TDenseVector()) : CMSSM(_nL, _nTL, _nS, _nZ, _nY, _nU, _nE, _nExpectation, ref, sef, mef, tpf, pdf, _x), location(_loc) {}

	virtual int LogLikelihood(double &log_likelihood, vector<TDenseVector> &z_tm1_last, vector<TDenseMatrix> &P_tm1_last, TDenseVector &p_tm1_last, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob); 
	virtual int LogPosterior(double &log_posterior, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob);
};

#endif
