#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <ctime>
#include <getopt.h>
#include <cstdlib>
#include "dw_rand.h"
#include "CMSSM.hpp"
#include "CMSSM_Error_Code.hpp"
#include "RationalExpectationFunction_test.hpp"
#include "StateEquationFunction_test.hpp"
#include "MeasurementEquationFunction_test.hpp"
#include "TransitionMatrixFunction_test.hpp" 
#include "MakeLbUb_test.hpp"
#include "ReadWriteFile.hpp"
#include "CMSSM_test_1st.hpp"
#include "CMSSM_test_2nd.hpp"
#include "PriorDistrFunction_test_1st.hpp"
#include "PriorDistrFunction_test_2nd.hpp"
#include "PriorDistrFunction_test_all.hpp"
#include "InitializeParameter_test.hpp" 

using namespace std; 

int main(int argc, char **argv)
{
	static struct option long_options[] = 
	{
		{"data_file", required_argument, 0, 'd'}, 
		{"initial_value_file", required_argument, 0, 'v'}, 
		{"max_tries", required_argument, 0, 't'},
		{"minus_infinity", required_argument, 0, 'i'},
		{"output_file", required_argument, 0, 'f'},
		{0, 0, 0, 0}
	}; 
	int option_index = 0; 
	string data_file_name, initial_value_file_name, output_file_name; 
	size_t n_tries = 10; 
	double minus_infinity = -1.0e30; 
	while (1)
	{
		int c = getopt_long(argc, argv, "d:v:t:i:f:", long_options, &option_index); 
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
			case 'i':
				minus_infinity = atof(optarg); break; 
			case 'f': 
				output_file_name = string(optarg); break; 
			default:
				break; 
		}
	}
	if (data_file_name.length() ==0) 
	{
		cerr << "Usage: " << argv[0] << " -d data file -v initial value file -t maximum tries.\n"; 
		abort(); 
	}

	// random number generator
	dw_initialize_generator(time(NULL));

	// Hubrid NK model (standard reduced form)
	size_t nFree = 21+8+2;	// 21 (15+6): model parameters; 6=2X3 (Delta in Dan's sunspont notation) + 2=2x1 (gamma in Dan's sunspot notation): sunspot parameters (with 2 endogeneous errors and 3 fundamental shocks); 2 probabilites of staying in ZLB
	size_t n1 = 14; 	// Number of parameters staying in Regime 1
	vector<int> locs_x1, locs_x2, locs_xall;	// locations of parameters for 1st and 2nd regimes, and for all parameters 
	for (unsigned int i=0; i<n1; i++)
	{
		locs_x1.push_back(i); 
		locs_xall.push_back(i); 
	}
	for (unsigned int i=n1; i<n1+7; i++)
	{
		locs_xall.push_back(i); 
		locs_x2.push_back(i); 
	}
	for (unsigned int i=n1+7+5; i<n1+7+5+3; i++)
		locs_x2.push_back(i); 
	locs_xall.push_back(n1+6+2); 
	for (unsigned int i=n1+6+2+4; i<n1+7+2+7; i++)
		locs_xall.push_back(i); 
	
	// fixed_parameters for (RationalExpectationFunction, StateEquationFunction, MeasurementEquationFunction) and TransitionProbMatrix
	TDenseVector fixed_parameter(3);
	fixed_parameter(0) = 1.0e+05; 	// scl4gsigmai2: Scale for s.d. of interest rate in ZLB. Scale for s.d. of interest rate in ZLB. 1) 1 (or 2.5) basis point annually = 2.5000e-05 (or 6.2500e-05) quarterly.  2) 2.5000e-04 quarterly = 10 basis points annually in the ZLB.
	fixed_parameter(1) = 0.000375;   	// ilow: Net rate: 10 basis points (0.1%) interest rate annually at zero bound or 0.025% quarterly at zero bound. In 2012Q2, it is 3.83e-04*4 = 0.0015 annually.
	fixed_parameter(2) = 0.005/4;  	// cut-off: Annual rate of 50 baisis points for the interest rate.

	// Initial parameters
	vector<TDenseVector> initialX; 
	TDenseVector x0; 
	if (LoadInitialValue(initialX, initial_value_file_name) == SUCCESS)
                x0 = initialX[0];
	else 
		x0 = InitializeParameter(nFree, fixed_parameter); 

	/* Not implementing
 * 		Calibrated equilibrium
 * 		Impulse responsese based on x0
 	*/	

	/* Read in data: begin */
	TDenseVector qm_date;		// dates
	vector<TDenseVector> qdata; 	// variables
	if (LoadData(qm_date, qdata, data_file_name) != SUCCESS)
	{
		cerr << "------LoadData(): error occurred ------\n"; 
		abort(); 
	}
	/* Read in data: end */
	size_t nSample = qdata.size(); 	// total sample size
	size_t nY = qdata[0].dim;	// number of observables
	vector<TDenseVector> y = qdata; 
	/* Switch order of qdata, so that the new order is [1 0 2]
  	coressponding to inflation, output and interest rate, respectively.
	vector<TDenseVector> y = vector<TDenseVector>(nSample); 
	vector<int>new_order(nY); 
	new_order[0] = 1; new_order[1] = 0; new_order[2] = 2; 
	for (unsigned int i=0; i<nSample; i++)
		y[i] = qdata[i].SwitchOrder(new_order); 
	Switch order of qdata: end */
	
	// First regime (normal, above 50 basis points at annual rate)
	int yrStart = 1985, qmStart = 1, yrEnd = 2008, qmEnd = 3; 
	int y1_length = (yrEnd-yrStart)*4 + (qmEnd-qmStart) +1; 
	vector<TDenseVector> y1st(y.begin(), y.begin()+y1_length); 	
	// Second regime (ZLB, about 2.5 basis points quarterly and 10 basis points annually)
	yrStart = 2009, qmStart = 1, yrEnd = 2012, qmEnd = 2; 
	int y2_length = (yrEnd-yrStart)*4 + (qmEnd-qmStart) + 1; 
	vector<TDenseVector> y2nd(y.end()-y2_length, y.end()); 	
	
	// sizes
	size_t nS = 4; 
	size_t nL = 1; 
	size_t nTL = 1; 
	size_t nZ = 7; 
	size_t nE = 3; 
	size_t nU = 0; 
	size_t nExpectation = 2; 

	// Setting up solutions
	CMSSM::MINUS_INFINITY_LOCAL = minus_infinity; 

	// Setting up the CMSSM model
	RationalExpectationFunction_test *rational_expectation_func = new RationalExpectationFunction_test(fixed_parameter); 
	StateEquationFunction_test *state_equation_func = new StateEquationFunction_test(fixed_parameter); 
	MeasurementEquationFunction_test *measurement_equation_func = new MeasurementEquationFunction_test(fixed_parameter); 
	TransitionProbMatrixFunction_test *transition_prob_func = new TransitionProbMatrixFunction_test(fixed_parameter); 
	CMSSM_test_1st model_1st(locs_x1, nL, nTL, nS, nZ, nY, nU, nE, nExpectation, rational_expectation_func, state_equation_func, measurement_equation_func, transition_prob_func, new PriorDistrFunction_test_1st(fixed_parameter) ); 
	
	CMSSM_test_2nd model_2nd(locs_x2, nL, nTL, nS, nZ, nY, nU, nE, nExpectation, rational_expectation_func, state_equation_func, measurement_equation_func, transition_prob_func, new PriorDistrFunction_test_2nd(fixed_parameter) ); 
	
	CMSSM_test_1st model_all(locs_xall, nL, nTL, nS, nZ, nY, nU, nE, nExpectation, rational_expectation_func, state_equation_func, measurement_equation_func, transition_prob_func, new PriorDistrFunction_test_all(fixed_parameter) ); 

	// 1st regime (s_t*, s_{t-1}*) = (1*, 1*)
	// Starts in state one with probability one
	TDenseVector initial_prob_1st(model_1st.nXi, 0.0);
	initial_prob_1st(0) = 1.0;  
	// Initial values for initial Kalman filter 
	vector<TDenseVector> z0_1st(model_1st.nXi);
	vector<TDenseMatrix> P0_1st(model_1st.nXi); 
	for (unsigned int i=0; i<model_1st.nXi; i++)
	{
		z0_1st[i].Zeros(model_1st.nZ); 
		P0_1st[i]=Identity(model_1st.nZ)*100.0;
	}
	double log_likelihood_1st, log_likelihood_2nd, log_likelihood_all; 
	double log_posterior_1st, log_posterior_2nd, log_posterior_all; 

	vector<TDenseVector> z_tm1_last_1st, z_tm1_last_2nd, z_tm1_last_all; 
       	vector<TDenseMatrix> P_tm1_last_1st, P_tm1_last_2nd, P_tm1_last_all; 
       	TDenseVector p_tm1_last_1st, p_tm1_last_2nd, p_tm1_last_all; 
	
	TDenseVector initial_prob_2nd(model_2nd.nXi); 
	initial_prob_2nd.Zeros(); 
       	initial_prob_2nd(10) = 1.0;     // 10: the position for (3*, 3*) (staying in the ZLB regime)
	vector<TDenseVector> z0_2nd(model_2nd.nXi); 
	vector<TDenseMatrix> P0_2nd(model_2nd.nXi); 

	/* Just for debugging
	// Calling LogLikelihood to get z_tm1_last, P_tm1_last and p_tm1_last to be used by 2nd regime
	if (model_1st.LogLikelihood(log_likelihood_1st, z_tm1_last_1st, P_tm1_last_1st, p_tm1_last_1st, x0, y1st, z0_1st, P0_1st, initial_prob_1st) == SUCCESS)
		cout << "LogLikelihood 1st regime: " << log_likelihood_1st << endl; 
	else 
		cerr << "LogLikelihood 1st regime: Error occurred \n"; 
	double log_posterior_1st; 
	if (model_1st.LogPosterior(log_posterior_1st, x0, y1st, z0_1st, P0_1st, initial_prob_1st) == SUCCESS)
		cout << "LogPosterior 1st regime: " << log_posterior_1st << endl; 
	else 
		cerr << "LogPosterior 1st regime: Error occurred \n"; 

	// 2nd regime
	for (unsigned int i=0; i<model_2nd.nXi; i++)
	{
		z0_2nd[i].CopyContent(z_tm1_last_1st[0]); 
		P0_2nd[i].CopyContent(P_tm1_last_1st[0]); 
	} 
	// Just for debugging
	if (model_2nd.LogLikelihood(log_likelihood_2nd, z_tm1_last_2nd, P_tm1_last_2nd, p_tm1_last_2nd, x0, y2nd, z0_2nd, P0_2nd, initial_prob_2nd) == SUCCESS) 
		cout << "LogLikelihood 2nd regime: " << log_likelihood_2nd << endl; 
	else 
		cerr << "LogLikelihood 2nd regime: Error occurred\n"; 
	if (model_2nd.LogPosterior(log_posterior_2nd, x0, y2nd, z0_2nd, P0_2nd, initial_prob_2nd) == SUCCESS)
		cout << "LogPosterior 2nd regime: " << log_posterior_2nd << endl; 
	else 
		cerr << "LogPosterior 2nd regime: Error occurred\n";  */

	// all together
	if (model_all.LogLikelihood(log_likelihood_all, z_tm1_last_all, P_tm1_last_all, p_tm1_last_all, x0, y, z0_1st, P0_1st, initial_prob_1st) == SUCCESS) 
		cout << "LogLikelihood all : " << log_likelihood_all << endl; 
	else 
		cerr << "LogLikelihood all : Error occurred\n"; 
	if (model_all.LogPosterior(log_posterior_all, x0, y, z0_1st, P0_1st, initial_prob_1st) == SUCCESS)
		cout << "LogPosterior all: " << log_posterior_all << endl; 
	else 
		cerr << "LogPosterior all: Error occurred\n"; 

	// Searching for optimal parameters
	TDenseVector best_solution(nFree+2, 0.0);
	best_solution(0) = best_solution(1) = CMSSM::MINUS_INFINITY_LOCAL; 
	vector<TDenseVector> solutions(n_tries); 
	for (unsigned int i=0; i<solutions.size(); i++)
	{
		// solutions[i]: (MINUS_INFINITY, MINUS_INFINITY, x0)
		solutions[i] = TDenseVector(nFree+2); 
		solutions[i][0] = solutions[i][1] = CMSSM::MINUS_INFINITY_LOCAL; 
		for (int unsigned j=0; j<x0.dim; j++)
			solutions[i][j+2] = x0[j]; 
	} 

	size_t number_bad = 0; 
	TDenseVector x_input, initial_x1, initial_x2, initial_xall, lb, ub;  
	TDenseVector xOptimal_1st, xOptimal_2nd, xOptimal_all; 
	double lpOptimal_1st, lpOptimal_2nd, lpOptimal_all; 
	unsigned int i=0; 
	bool bad; 
	
	x_input = x0; 
	while (i < n_tries  && number_bad < 200)
	{
		solutions[i].SetElement(CMSSM::MINUS_INFINITY_LOCAL, 0); 
		bad = false; 
	
		//////////////////////////////////////////////////////
		// First block: 1st regime parameters only 
		/////////////////////////////////////////////////////
	
		// initialize
		if (!bad)
		{
			initial_x1 = x_input;
			// only elements at locs_x1 are changed by RandomInit_test_1st, 
			// while all the other elements are as x_input.
			if (!bad && i>0 && RandomInit_test_1st(initial_x1, locs_x1) != SUCCESS)
			{
				cerr << "Randomly initialize parameters of the 1st regime: Error occurred.\n";
				bad = true; 
			}	
		} 
		// bound 
		if(!bad && MakeLbUb_test(lb,ub,initial_x1.dim) != SUCCESS )
		{
			cerr << "MakeLbUb_test: Error occurred.\n"; 
			bad = true; 
		}
        	// log-posterior for initial_x1
		if (!bad && model_1st.LogPosterior(log_posterior_1st, initial_x1, y1st, z0_1st, P0_1st, initial_prob_1st) != SUCCESS)
		{
                	cerr << "LogPosterior 1st regime: Error occurred \n";
			bad = true; 
		}
		// optimizing 
		if (!bad && log_posterior_1st > CMSSM::MINUS_INFINITY_LOCAL)
			if  (model_1st.Maximize_LogPosterior_NPSOL(lpOptimal_1st,xOptimal_1st,y1st,z0_1st,P0_1st,initial_prob_1st,initial_x1) == MODEL_NOT_PROPERLY_SET)
			{
				cerr << " No eqilibrium for the 1st regime.\n";
				bad = true; 
			}
		// Update x_input(locs_x1) 
		if (!bad)
		{
			for (unsigned int j=0; j<locs_x1.size(); j++)
				x_input(locs_x1[j]) = xOptimal_1st(locs_x1[j]); 
		}

		//////////////////////////////////////////////////////
		// Second block: 2nd regime parameters only
		//////////////////////////////////////////////////////

		// obtain z_tm1_last_1st and P_tm1_last_1st to intialize z_02nd and P0_2nd
        	if (!bad && model_1st.LogLikelihood(log_likelihood_1st, z_tm1_last_1st, P_tm1_last_1st, p_tm1_last_1st, x0, y1st, z0_1st, P0_1st, initial_prob_1st) != SUCCESS)
		{
                	cerr << "LogLikelihood 1st regime: Error occurred"  << endl; 
			bad = true; 
		}

		// initiaize
		if (!bad)
		{	
        		for (unsigned int j=0; j<model_2nd.nXi; j++)
        		{
               			z0_2nd[j].CopyContent(z_tm1_last_1st[0]); 
               			P0_2nd[j].CopyContent(P_tm1_last_1st[0]); 
        		}
			initial_x2 = x_input; 
			// Only initial_x2(locs_x2) are changed by RandomInit_test_2nd
			if (i>0 && RandomInit_test_2nd(initial_x2, locs_x2) != SUCCESS)
			{
				cerr << "Randomly initialize parameters for the 2nd regime: error occurred.\n"; 
				bad = true; 
			}
		}
		// bounds
		if( !bad && MakeLbUb_test(lb,ub,initial_x2.dim) != SUCCESS )
		{
                        cerr << "MakeLbUb_test: Error occurred.\n";
			bad = true; 
		}
		// log-posterior for initial_x2
		if (!bad && model_2nd.LogPosterior(log_posterior_2nd, initial_x2, y2nd, z0_2nd, P0_2nd, initial_prob_2nd) != SUCCESS)
		{
                        cerr << "LogPosterior 2nd regime: Error occurred \n";
			bad = true;  
		}
		// optimizing
		if (!bad && log_posterior_2nd > CMSSM::MINUS_INFINITY_LOCAL)
			if (model_2nd.Maximize_LogPosterior_NPSOL(lpOptimal_2nd,xOptimal_2nd,y2nd,z0_2nd,P0_2nd,initial_prob_2nd,initial_x2) == MODEL_NOT_PROPERLY_SET)
			{
				cout << "No equilibrium for the 2nd regime.\n"; 
				bad = true; 
			}
		
		// Update x_input(locs_x2)
		if (!bad)
		{
			for (unsigned int j=0; j<locs_x2.size(); j++)
				x_input(locs_x2[j]) = xOptimal_2nd(locs_x2[j]); 			
		}
	
		////////////////////////////////////////////////////////////
		// Overall: all parameters
		///////////////////////////////////////////////////////////
		
		// initialize and bound
		if (!bad)
		{
			initial_xall = x_input; 
			if( !bad && MakeLbUb_test(lb,ub,initial_xall.dim) != SUCCESS )
			{
                        	cerr << "MakeLbUb_test: Error occurred.\n";	
				bad = true; 
			}
		}
		// log-posterior for initial_xall
		if (!bad && model_all.LogPosterior(log_posterior_all, initial_xall, y, z0_1st, P0_1st, initial_prob_1st) != SUCCESS)
		{
                        cerr << "LogPosterior all: Error occurred \n";
			bad = true; 
		}
		// optimize
		if (!bad && log_posterior_all > CMSSM::MINUS_INFINITY_LOCAL)
			if (model_all.Maximize_LogPosterior_NPSOL(lpOptimal_all,xOptimal_all,y,z0_1st,P0_1st,initial_prob_1st,initial_xall) == MODEL_NOT_PROPERLY_SET)
			{
				cerr << "No equilibrium for the overall regime.\n"; 
				bad = true; 
				number_bad ++; 
			}
	
		if (!bad)
		{
			// update x_input(locs_xall)
			for (unsigned int j=0; j<locs_xall.size(); j++)
				x_input(locs_xall[j]) = xOptimal_all(locs_xall[j]); 
			
			// log-likelihood for final solution x_input
			if (!bad && model_all.LogLikelihood(log_likelihood_all, z_tm1_last_all, P_tm1_last_all, p_tm1_last_all, x_input, y, z0_1st, P0_1st, initial_prob_1st) != SUCCESS)
			{
				cerr << "LogLikelihood all: Error occurred\n"; 
				bad = true; 
			}
			// log_posterior for final solution x_input
			if (!bad && model_all.LogPosterior(log_posterior_all, x_input, y, z0_1st, P0_1st, initial_prob_1st) != SUCCESS)
			{
				cerr << "LogPosterior all: Error occurred\n"; 
				bad = true; 
			}
			solutions[i](0)=log_posterior_all; 
			solutions[i](1)=log_likelihood_all; 
			for (unsigned int j=0; j<x_input.dim; j++)
				solutions[i][j+2] = x_input[j]; 
			if (log_posterior_all > best_solution[0])
				best_solution = solutions[i]; 
		}
		i++; 
	}		

	// Output results
	ofstream output_file; 
	if (output_file_name.length() > 0)
		output_file.open(output_file_name.c_str()); 
	if (output_file_name.length() > 0 && output_file)
	{
		/*output_file << "================== best solution =====================\n";
                for (unsigned int i=0; i<best_solution.dim; i++)
                	output_file << setprecision(20) << best_solution[i] << " ";
                output_file << endl;*/
        
               	 for (unsigned int i=0; i<n_tries; i++)
                {
                        // output_file << "\n ================ solution: " << i << "=====================\n";
                        for (unsigned int j=0; j<solutions[i].dim; j++)
                                output_file << setprecision(20) << solutions[i][j] << " ";
                        output_file << endl;
                }
		output_file.close(); 
	}
	else 
	{
		/*cout << "================== best solution =====================\n"; 
		for (unsigned int i=0; i<best_solution.dim; i++)
			cout << setprecision(20) << best_solution[i] << " "; 
		cout << endl; */

		for (unsigned int i=0; i<n_tries; i++)
		{
			// cout << "\n ================ solution: " << i << "=====================\n"; 
			for (unsigned int j=0; j<solutions[i].dim; j++)
				cout << setprecision(20) << solutions[i][j] << " "; 
			cout << endl; 
		}
	}
}
