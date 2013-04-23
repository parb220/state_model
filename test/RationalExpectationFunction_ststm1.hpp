#ifndef _MAKE_A_B_PSI_PI_C_STSTM1
#define _MAKE_A_B_PSI_PI_C_STSTM1

#include "CMSSM.hpp"
#include "findequilibrium.hpp"

class RationalExpectationFunction_ststm1 : public RationalExpectationFunction
{
        virtual int convert(vector<vector<TDenseMatrix> > &A, vector<vector<TDenseMatrix> > &B, vector<vector<TDenseMatrix> > &Psi, vector<vector<TDenseMatrix> >&Pi, vector<vector<TDenseVector> > &C, const TDenseVector &free_parameter, const TDenseVector &input_x);
};

#endif
