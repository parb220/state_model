#ifndef _RATIONALEXPECTATIONFUNCTION_TEST_
#define _RATIONALEXPECTATIONFUNCTION_TEST_

#include "CMSSM.hpp"

class RationalExpectationFunction_test : public RationalExpectationFunction
{
protected:
	virtual void ConvertXtoParameter(const TDenseVector &x); 
	double scl4Sigmai2;
	double ilow;
	double cutoff;
	double DoubledPistar;
	double iss;
	double Thetaf;
	double Tau;
	double Omegaf;
	double Kappac;
	double PhiDoubledPi;
	double Phiy;
	double RhoMu1;
	double Rhorn1;
	double Rhoi1;
	double SigmaMu1;
	double Sigmarn1;
	double Sigmai1;
	double RhoMu2;
	double Rhorn2;
	double Rhoi2;
	double SigmaMu2;
	double Sigmarn2;
	double Sigmai2;
	double Iota;
	double Thetab;
	double Omegab;
	double istar;

public:
	virtual int convert(vector<vector<TDenseMatrix> > &A, vector<vector<TDenseMatrix> > &B, vector<vector<TDenseMatrix> > &Psi, vector<vector<TDenseMatrix> >&Pi, vector<vector<TDenseVector> > &C, const TDenseVector &input_x, size_t nZ, size_t nY, size_t nU, size_t nE, size_t nExpectation);
	RationalExpectationFunction_test() : RationalExpectationFunction() {}
	RationalExpectationFunction_test(const TDenseVector &p) : RationalExpectationFunction(p) {}	// fixed_parameter[0]: scl4gsigmai; fixed_parameter[1]: ilow
	~RationalExpectationFunction_test() {}
};

#endif
