#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include "CMSSM_Error_Code.hpp"
#include "findequilibrium.hpp"
#include "RationalExpectationFunction_ststm1.hpp"
#include "StateEquationFunction_ststm1.hpp"
#include "MeasurementEquationFunction_ststm1.hpp"
#include "TransitionMatrixFunction_ststm1.hpp" 
#include "MakeLbUb_ststm1.hpp"
#include "ReadWriteFile.hpp"

using namespace std; 

int main(int argc, char **argv)
{
	static struct option long_options[] = 
	{
		{"data_file", required_argument, 0, 'd'}, 
		{"initial_value_file", required_argument, 0, 'v'}, 
		{0, 0, 0, 0}
	}; 
	int option_index = 0; 
	string data_file_name, initial_value_file_name; 
	size_t n_tries = 10; 
	while (1)
	{
		int c = getopt_long(argc, argv, "d:v:", long_options, &option_index); 
		if (c == -1)
			break; 
		switch(c)
		{
			case 'd':
				data_file_name = string(optarg); break; 
			case 'v':
				initial_value_file_name = string(optarg); break; 
			default:
				break; 
		}
	}
	if (data_file_name.length() == 0 || initial_value_file_name.length() == 0) 
	{
		cerr << "Usage: " << argv[0] << " -d data file -v initial value file.\n"; 
		abort(); 
	}

	// prefixed parameters for state equations
	TDenseVector state_equation_parameter(5);
	state_equation_parameter.SetElement(1.005, GLAMBDAZ_STATE);	// glambdaz, Gross rate: (1+2%) annually per capita or (1+0.5%) quarterly per capita 
	state_equation_parameter.SetElement(0.01, RN_STATE);	// rn, Net rate: rn=4% annually or 1% quarterly (natural rate of interest or steady state real interest rate)
	state_equation_parameter.SetElement(0.005, GPISTAR_STATE);	// gpistar, Net rate: log(1.02) -- 2% annually or 0.5% quarterly
	state_equation_parameter.SetElement(0.7, SC_STATE);	// sc, Steady state share of private consumption in C+G
	state_equation_parameter.SetElement(0.00025, RLOW_STATE);	// Rlow, Net rate 10 basis points (0.1%) interest rate annually at zero bound or 0.025% quarterly at zero bound. In 2012Q2, it is 3.83e04*4-0.0015 annually
	
	// prefixed parameters for measurement equations
	TDenseVector measurement_equation_parameter(1); 
	measurement_equation_parameter.SetElement(0.005,GPISTAR_MEASUREMENT); 	// measurment_equation_parameter[0] = gpistar	

	// prefixed parameters for transition prob matrix
	TDenseVector transition_prob_parameter(1); 
	transition_prob_parameter.SetElement(0.005/4, CUTOFF_TRANSITION);	// transition_prob_parameter[0] = cutoff; 

	/* Read in data: begin */
	TDenseVector qm_date;		// dates
	vector<TDenseVector> qdata; 	// variables
	if (LoadData(qm_date, qdata, data_file_name) != SUCCESS)
	{
		cerr << "------LoadData(): error occurred ------\n"; 
		abort(); 
	}
	size_t nSample = qdata.size(); 	// total sample size
	size_t nY = qdata[0].dim;	// number of observables
	/* Read in data: end */
	
	// sizes
	size_t nS = 4; 
	size_t nL = 1; 
	size_t nTL = 1; 

	// Setting up the CMSSM model
	CMSSM model(nL, nTL, nS, state_equation_parameter, measurement_equation_parameter, transition_prob_parameter, TDenseVector(), new RationalExpectationFunction_ststm1, new StateEquationFunction_ststm1, new MeasurementEquationFunction_ststm1, new TransitionProbMatrixFunction_ststm1, NULL); 

	size_t nFree = 23+8+2; // 23: model parameters; 6=2X3+2: sunspot parameters with 2 endogenous errors and 3 fundamental shocks; 2: probability of staying in ZLB
	size_t max_count; 
	vector<TDenseVector> initialX;

	if (LoadInitialValue(initialX, initial_value_file_name) != SUCCESS)
	{
		cerr << "------LoadInitialValue(): error occurred ------\n"; 
		abort(); 
	}
	model.UpdateStateModelParameters(0,qdata,initialX[0]);

	// Starts in state one with probability one
	TDenseVector initial_prob(model.nXi, 0.0);
	initial_prob.SetElement(1.0, 0); 
	 
	// Initial values for initial Kalman filter 
	vector<TDenseVector> z0(model.nXi);
	vector<TDenseMatrix> P0(model.nXi); 
	for (unsigned int i=0; i<model.nXi; i++)
	{
		z0[i].Zeros(model.nZ); 
		P0[i].Identity(model.nZ);
		P0[i].Multiply(P0[i], 100); 
	}

	/* Calculate log-likelihood */
	double log_likelihood; 
	for (unsigned int i=0; i<initialX.size(); i++)
	{
		if (model.LogLikelihood(log_likelihood, initialX[i], qdata, z0, P0, initial_prob) == SUCCESS)
			cout << "Log Likelihood " << log_likelihood << endl; 
		else 
			cout << "Error occurred when calculating log_likelihood.\n"; 
	
	}
}
