#ifndef _MEASUREMENTEQUATIONFUNCTION_TEST_
#define _MEASUREMENTEQUATIONFUNCTION_TEST_

#include "CMSSM.hpp"

class MeasurementEquationFunction_test : public MeasurementEquationFunction
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
	virtual int convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &input_x, size_t nZ, size_t nY, size_t nU, size_t nNu);
	MeasurementEquationFunction_test() : MeasurementEquationFunction() {}
	MeasurementEquationFunction_test(const TDenseVector &p) : MeasurementEquationFunction(p) {}	// fixed_parameter[0]: scl4gsigmai fixed_parameter[1]: ilow
	virtual ~MeasurementEquationFunction_test() {}
};

#endif
