#ifndef _CMSSM_TEST_2ND_
#define _CMSSM_TEST_2ND_
#include "CMSSM_test.hpp"

class CMSSM_test_2nd : public CMSSM_test
{
public:
	CMSSM_test_2nd() : CMSSM_test(){}
	CMSSM_test_2nd(const vector<int> &_loc_v, const vector<int> &_loc_c) : CMSSM_test(_loc_v, _loc_c) {}
	CMSSM_test_2nd(const vector<int> &_loc_v, const vector<int> &_loc_c, size_t _nL, size_t _nTL, size_t _nS, size_t _nZ, size_t _nY, size_t _nU, size_t _nE, size_t _nExpectation, RationalExpectationFunction * ref=NULL, StateEquationFunction * sef=NULL, MeasurementEquationFunction * mef=NULL, TransitionProbMatrixFunction * tpf=NULL, PriorDistributionFunction * pdf= NULL, const TDenseVector & _x=TDenseVector()) : CMSSM_test(_loc_v, _loc_c, _nL, _nTL, _nS, _nZ, _nY, _nU, _nE, _nExpectation, ref, sef, mef, tpf, pdf, _x) {}
	virtual ~CMSSM_test_2nd(){}
	virtual int LogLikelihood(double &log_likelihood, vector<TDenseVector> &z_tm1_last, vector<TDenseMatrix> &P_tm1_last, TDenseVector &p_tm1_last, const TDenseVector &x_variable, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob); 
};

#endif
