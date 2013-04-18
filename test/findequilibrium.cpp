#include <cstdlib>
#include "findequilibrium.hpp"
#include "CMakeABPsiPiC_ststm1.hpp"
#include "CMeasurementEquationFunction_test.hpp"
#include "CTransitionMatrix_test.hpp"
#include "MakeLbUb_ststm1.hpp"
#include "ReadWriteFile.hpp"

using namespace std; 

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cerr << "Usage: " << argc << "file-containing-the-data" << endl; 
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
	string filename = string(argv[1]);	
	vector<double> qm_date;		// dates
	vector<TDenseVector> qdata; 	// variables
	if (LoadData(qm_date, qdata, filename))
	{
		cerr << "------ Error occurred when reading data ------\n"; 
		abort(); 
	}
	size_t nSample = qdata.size(); 	// total sample size
	size_t nY = qdata[0].dim;	// number of observables
	/* Read in data: end */
	
	// sizes
	size_t nS = 4; 
	size_t nL = 1; 
	size_t nTL = 1; 
	size_t nFree = 23+8+2; // 23: model parameters; 6=2X3+2: sunspot parameters with 2 endogenous errors and 3 fundamental shocks; 2: probability of staying in ZLB

	// Setting up the CMSSM model
	CMSSM model(nL, nTL, nS, state_equation_parameter, measurement_equation_parameter, transition_prob_parameter); 
	model.state_equation_function = new MakeABPsiPiC_ststm1; 
	model.measurement_equation_function = new MeasurementEquationFunction_test;  
	model.transition_prob_function = new CTransitionProbMatrixFunction_Test; 	

	/* Find valid starting value: begin */
	ObjectiveFunction_Validation::ABPsiPiC_function = &function_ABPsiPiC; 
	ObjectiveFunction_Validation::state_equation_parameter = state_equation_parameter; 
	ObjectiveFunction_Validation::measurement_equation_parameter = measurement_equation_parameter; 
	ObjectiveFunction_Validation::transition_prob_parameter = transition_prob_parameter; 

	// Lower and upper bounds for x 
	TDenseVector lb, ub; 
	if ( MakeLbUb_ststm1(lb, ub, nFree) ) 	
	{
		cerr << "------ Error occurred when setting the lower and upper bounds of the parameter ------\n";
		abort(); 
	}

	// Find valid starting value for x
	TDenseVector x0(nFree);  x0.RandomNormal(nFree); // Initial guess of x
	TDenseVector function_value; 			// To hold function values of all max_count searches
	size_t max_count = 500; 
	if ( model.ValidInitialPoint(function_value, x0, max_count, lb, ub) )
	{
		cerr << "------ Attempting to find valid initial point for estimation: no equilibrium exists ------" << endl; 
		abort(); 
	}
	// Display valid x
	// Display(model.x); 
	/* Find valid start value: end */

	/* Set up state equation after a legitimate x is found*/
	// model.nZ and mode.nE should be set here
	model.StateEquationHelper(&function_ABPsiPiC); 
	
	/* Set up measurement equation */
	if( model.MeasurementEquation(&measurement_function) )
	{
		cerr << "Error occurred while setting up measurment equation.\n"; 
		abort(); 
	} 

	/* Starts in state one with probability one */
	TDenseVector initial_prob(model.nXi, 0.0);
	initial_prob.SetElement(1.0, 0); 
	 
	/* Initial values for initial Kalman filter */
	vector<TDenseVector> z0(model.nXi, TDenseVector(model.nZ,0.0) );
	vector<TDenseMatrix> P0(model.nXi, TDenseMatrix(model.nZ, model.nZ) ); 
	for (unsigned int i=0; i<P0.size(); i++)
		P0[i].Identity(model.nZ); 	
	
	/* Transition probability matrix: */

	/* Maximize log-likelihood */
	const double MINUS_INFINITY = -1.0e30; 
	size_t n_tries = 10; 
	TDenseVector best_solution(model.x.dim+1, 0.0);
	best_solution.SetElement(MINUS_INFINITY, 0); 
	vector<TDenseVector> solutions(n_tries, TDenseVector(model.x.dim+1,0.0) );  

	unsigned int i=0, number_bad=0; 
	while (i < n_tries && number_bad < 200)
	{
		solutions[i].SetElement(MINUS_INFINITY, 0); 
		bool bad = true; 
		
		// Find valid starting value
		x0.RandomNormal(nFree); 
		if ( (x0[x0.dim-1] <=0.0) || (x0[x0.dim-1] > 1.0) )
			x0.SetElement(0.5,x0.dim-1);
		if ( (x0[x0.dim-2] <=0.0) || (x0[x0.dim-2] > 1.0) )
			x0.SetElement(0.5,x0.dim-2); 

		max_count = 10; 
		if ( MakeLbUb_ststm1(lb, ub, nFree) )
        	{
                	cerr << "------ Error occurred when setting the lower and upper bounds of the parameter ------\n";
                	bad = true; 
        	}
		else if ( model.ValidInitialPoint(function_value, x0, max_count, lb, ub) )
        	{
                	cerr << "------ Attempting to find valid initial point for estimation: no equilibrium exists ------" << endl;
                	bad = true; 
        	}
		else 
		{
			// Is likelihood well defined at the initial value
			double ll = -minus_log_likelihood(qdata, z0, P0, &trans_prob_function, initial_prob, model);
			
			// unconstrained optimization
			

			// save solution
			double mll = minus_log_likelihood(qdata, z0, P0, &trans_prob_function, initial_prob, model);  
			ll = -mll; 
			solutions[i].SetElement(ll, 0); 
			for (unsigned int j=1; j<solutions[i].dim; j++)
				solutions[i].SetElement(model.x[j-1], j); 

			if (ll > best_solution[0])
			{
				for (unsigned int j=0; j<solutions[i].dim; j++)
					best_solution.SetElement(solutions[i][j],j); 
			}
			i ++; 
			bad = false; 
		}
		if (bad)
			number_bad ++;
	}
}
