#include <vector>
#include "dw_dense_matrix.hpp"
#include "CMSSM_Error_Code.hpp"
#include "CMSSM_test.hpp"
#include "CMSSM_test_2nd.hpp"
#include "ReadWriteFile.hpp"
#include "RationalExpectationFunction_test.hpp"
#include "StateEquationFunction_test.hpp"
#include "TransitionMatrixFunction_test.hpp"
#include "MeasurementEquationFunction_test.hpp"
#include "PriorDistrFunction_test_1st.hpp"
#include "PriorDistrFunction_test_2nd.hpp"
#include "PriorDistrFunction_test_all.hpp"

using namespace std; 

bool ReadOutData(const string &data_file_name, vector<TDenseVector> &y1st, vector<TDenseVector> &y2nd, vector<TDenseVector> &y)
{
	/* Read in data: begin */
	TDenseVector qm_date;		// dates
	vector<TDenseVector> qdata; 	// variables
	if (LoadData(qm_date, qdata, data_file_name) != SUCCESS)
		return false; 

	/* Read in data: end */
	y = qdata; 
	
	// First regime (normal, above 50 basis points at annual rate)
	int yrStart = 1985, qmStart = 1, yrEnd = 2008, qmEnd = 3; 
	int y1_length = (yrEnd-yrStart)*4 + (qmEnd-qmStart) +1; 
	y1st=vector<TDenseVector>(y.begin(), y.begin()+y1_length); 	
	// Second regime (ZLB, about 2.5 basis points quarterly and 10 basis points annually)
	yrStart = 2009, qmStart = 1, yrEnd = 2012, qmEnd = 2; 
	int y2_length = (yrEnd-yrStart)*4 + (qmEnd-qmStart) + 1; 
	y2nd=vector<TDenseVector>(y.end()-y2_length, y.end()); 	
	return true; 
}

void SetUpModel(size_t nFree, size_t nY, TDenseVector &fixed_parameter, CMSSM_test &model_1st, CMSSM_test_2nd &model_2nd, CMSSM_test &model_all) 
{
	// fixed_parameters for (RationalExpectationFunction, StateEquationFunction, MeasurementEquationFunction) and TransitionProbMatrix
	fixed_parameter.Resize(3);
	fixed_parameter(0) = 1.0e+05; 	// scl4gsigmai2: Scale for s.d. of interest rate in ZLB. Scale for s.d. of interest rate in ZLB. 1) 1 (or 2.5) basis point annually = 2.5000e-05 (or 6.2500e-05) quarterly.  2) 2.5000e-04 quarterly = 10 basis points annually in the ZLB.
	fixed_parameter(1) = 0.000375;   	// ilow: Net rate: 10 basis points (0.1%) interest rate annually at zero bound or 0.025% quarterly at zero bound. In 2012Q2, it is 3.83e-04*4 = 0.0015 annually.
	fixed_parameter(2) = 0.005/4;  	// cut-off: Annual rate of 50 baisis points for the interest rate.
	
	// Hubrid NK model (standard reduced form)
	size_t n1 = 14; 	// Number of parameters staying in Regime 1
	vector<int> locs_variable_x1, locs_variable_x2, locs_variable_xall;     // locations of variable parameters for 1st and 2nd regimes, and for all regimes
	vector<int> locs_constant_x1, locs_constant_x2, locs_constant_xall;     // locations of invariable parameters for 1st and 2nd regimes, and for all regimes
	// locs_variable_x1 : [0:n1); 
	for (unsigned int i=0; i<n1; i++)
                locs_variable_x1.push_back(i);
        for (unsigned int i=n1; i<nFree; i++)
                locs_constant_x1.push_back(i);

	// locs_variable_x2 : [n1:n1+7), [n1+7+5: n1+7+5+3)
	for (unsigned int i=0; i<n1; i++)
                locs_constant_x2.push_back(i);
        for (unsigned int i=n1; i<n1+7; i++)
                locs_variable_x2.push_back(i);
        for (unsigned int i=n1+7; i<n1+7+5; i++)
                locs_constant_x2.push_back(i);
        for (unsigned int i=n1+7+5; i<n1+7+5+3; i++)
                locs_variable_x2.push_back(i);
        for (unsigned int i=n1+7+5+3; i<nFree; i++)
                locs_constant_x2.push_back(i);

	// locs_variable_xall : [1:n1+7), n1+7+1, [n1+7+5,n1+7+5+4)
	for (unsigned int i=0; i<n1+7; i++)
                locs_variable_xall.push_back(i);
        locs_constant_xall.push_back(n1+7);
        locs_variable_xall.push_back(n1+7+1);
        for (unsigned int i=n1+7+2; i<n1+7+5; i++)
                locs_constant_xall.push_back(i);
        for (unsigned int i=n1+7+5; i<n1+7+9; i++)
                locs_variable_xall.push_back(i);
        for (unsigned int i=n1+7+9; i<nFree; i++)
                locs_constant_xall.push_back(i);
	
	// sizes
	size_t nS = 4; 
	size_t nL = 1; 
	size_t nTL = 1; 
	size_t nZ = 7; 
	size_t nE = 3; 
	size_t nU = 0; 
	size_t nExpectation = 2; 

	// Setting up the CMSSM model
	RationalExpectationFunction_test *rational_expectation_func = new RationalExpectationFunction_test(fixed_parameter); 
	StateEquationFunction_test *state_equation_func = new StateEquationFunction_test(fixed_parameter); 
	MeasurementEquationFunction_test *measurement_equation_func = new MeasurementEquationFunction_test(fixed_parameter); 
	TransitionProbMatrixFunction_test *transition_prob_func = new TransitionProbMatrixFunction_test(fixed_parameter); 

	model_1st=CMSSM_test(locs_variable_x1, locs_constant_x1, nL, nTL, nS, nZ, nY, nU, nE, nExpectation, rational_expectation_func, state_equation_func, measurement_equation_func, transition_prob_func, new PriorDistrFunction_test_1st(fixed_parameter)); 

	model_2nd=CMSSM_test_2nd(locs_variable_x2, locs_constant_x2, nL, nTL, nS, nZ, nY, nU, nE, nExpectation, rational_expectation_func, state_equation_func, measurement_equation_func, transition_prob_func, new PriorDistrFunction_test_2nd(fixed_parameter) ); 
	
	model_all=CMSSM_test(locs_variable_xall, locs_constant_xall, nL, nTL, nS, nZ, nY, nU, nE, nExpectation, rational_expectation_func, state_equation_func, measurement_equation_func, transition_prob_func, new PriorDistrFunction_test_all(fixed_parameter) );

	return;   
}
	
