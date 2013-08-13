#ifndef _CMSSM_TEST_1ST
#define _CMSSM_TEST_1ST
#include "CMSSM.hpp"

class CMSSM_test : public CMSSM
{
protected:
	vector<int> location_variable;
        vector<int> location_constant;
        TDenseVector x_variable_constant;
public:
	CMSSM_test() : CMSSM(), location_variable(), location_constant(), x_variable_constant(TDenseVector(0)) {}
	CMSSM_test(const vector<int> &_loc_v, const vector<int> &_loc_c) : CMSSM(), location_variable(_loc_v), location_constant(_loc_c), x_variable_constant(TDenseVector(_loc_v.size()+_loc_c.size(), 0.0)) {}
	CMSSM_test(const vector<int> &_loc_v, const vector<int> &_loc_c, size_t _nL, size_t _nTL, size_t _nS, size_t _nZ, size_t _nY, size_t _nU, size_t _nE, size_t _nExpectation, RationalExpectationFunction * ref=NULL, StateEquationFunction * sef=NULL, MeasurementEquationFunction * mef=NULL, TransitionProbMatrixFunction * tpf=NULL, PriorDistributionFunction * pdf= NULL, const TDenseVector & _x=TDenseVector()) : CMSSM(_nL, _nTL, _nS, _nZ, _nY, _nU, _nE, _nExpectation, ref, sef, mef, tpf, pdf, _x), location_variable(_loc_v), location_constant(_loc_c), x_variable_constant(TDenseVector(_loc_v.size()+_loc_c.size(),0.0)) {}
	virtual ~CMSSM_test(){}

	const vector<int> &locs_variable() const {return location_variable; }
        const vector<int> &locs_constant() const {return location_constant; }
	const TDenseVector & GetConstantPart() const {return x_variable_constant.SubVector(location_constant);}
        int SetConstantPart(const TDenseVector &_constant_x);

	virtual int LogLikelihood(double &log_likelihood, vector<TDenseVector> &z_tm1_last, vector<TDenseMatrix> &P_tm1_last, TDenseVector &p_tm1_last, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob);

	virtual int LogPosterior(double &log_posterior, const TDenseVector &x_variable, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob);
};

#endif
