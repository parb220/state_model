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
		{"max_tries", required_argument, 0, 't'},
		{0, 0, 0, 0}
	}; 
	int option_index = 0; 
	string data_file_name, initial_value_file_name; 
	size_t n_tries = 10; 
	while (1)
	{
		int c = getopt_long(argc, argv, "d:v:t:", long_options, &option_index); 
		if (c == -1)
			break; 
		switch(c)
		{
			case 'd':
				data_file_name = string(optarg); break; 
			case 'v':
				initial_value_file_name = string(optarg); break; 
			case 't': 
				n_tries = atoi(optarg); break; 
			default:
				break; 
		}
	}
	if (data_file_name.length() ==0) 
	{
		cerr << "Usage: " << argv[0] << " -d data file -v initial value file -t maximum tries.\n"; 
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
	TDenseVector x0(nFree), x0Valid(nFree), lb, ub, function_value;  

	if (LoadInitialValue(initialX, initial_value_file_name) == SUCCESS)
		x0Valid = initialX[0]; 
	else // Find valid starting value 
	{ 
		// Lower and upper bounds for x 
		if ( MakeLbUb_ststm1(lb, ub, nFree) != SUCCESS ) 	
		{
			cerr << "------MakeLbUb(): error occurred ------\n";
			abort(); 
		}

		// Find valid starting value for x 
		max_count = 500; 
		x0.RandomNormal(nFree); 			// Initial guess of x
		if (model.ValidInitialPoint(function_value, x0Valid, x0, max_count, lb, ub) != SUCCESS)
		{
			cerr << "------ CMSSM::ValidInitialPoint(): no equilibrium exists ------.\n"; 
			abort(); 
		}
	}
	model.UpdateStateModelParameters(0,qdata,x0Valid);

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
	}

	/* Maximize log-likelihood */
	const double MINUS_INFINITY = -1.0e30; 

	unsigned int i=0, number_bad=0; 
	bool bad; 
	double ll;	// log likelihood
	TDenseVector xOptimal(nFree); 

	TDenseVector best_solution(nFree+1, 0.0);
	best_solution.SetElement(MINUS_INFINITY, 0); 
	vector<TDenseVector> solutions; 
	bool if_search_initial;  
	if (initialX.empty())
	{
		initialX.resize(n_tries); 
		for (unsigned int i=0; i<n_tries; i++)
			initialX[i].RandomNormal(nFree);
		solutions.resize(n_tries);
		if_search_initial = true;  
	}
	else 
	{
		solutions.resize(initialX.size()); 
		if_search_initial = false; 
	}
	
	for (unsigned int i=0; i<solutions.size(); i++)
		solutions[i] = TDenseVector(nFree+1,0.0); 

	while (i < initialX.size() && number_bad < 200)
	{
		solutions[i].SetElement(MINUS_INFINITY, 0); 
		bad = true; 
		
		// Find valid starting value 
		x0 = initialX[i]; 
		if ( (x0[x0.dim-1] <=0.0) || (x0[x0.dim-1] > 1.0) )
			x0.SetElement(0.5,x0.dim-1);
		if ( (x0[x0.dim-2] <=0.0) || (x0[x0.dim-2] > 1.0) )
			x0.SetElement(0.5,x0.dim-2); 

		max_count = 10; 
		if (if_search_initial) 
		{
			if ( MakeLbUb_ststm1(lb, ub, nFree) == SUCCESS && model.ValidInitialPoint(function_value, x0Valid, x0, max_count, lb, ub) == SUCCESS && model.Maximize_LogLikelihood_NPSOL(ll,xOptimal,qdata,z0,P0,initial_prob,x0Valid) != MODEL_NOT_PROPERLY_SET ) 
				bad = false; 
			else 
				bad = true; 
		}
		else 
		{
			x0Valid = x0; 
			if (model.Maximize_LogLikelihood_NPSOL(ll,xOptimal,qdata,z0,P0,initial_prob,x0Valid) != MODEL_NOT_PROPERLY_SET )
				bad = false; 
			else 
				bad = true; 
				
		}
		if (bad )
			number_bad ++; 	
		else 
		{
			solutions[i].SetElement(ll, 0); 
			for (unsigned int j=1; j<solutions[i].dim; j++)
				solutions[i].SetElement(xOptimal[j-1], j); 
			if (ll > best_solution[0])
			{
				for (unsigned int j=0; j<solutions[i].dim; j++)
					best_solution.SetElement(solutions[i][j],j); 
			}
			i ++; 
			bad = false; 
		}
	}
	for (unsigned int i=0; i<best_solution.dim; i++)
		printf("%.20g ", best_solution[i]); 
	printf("\n"); 
}
