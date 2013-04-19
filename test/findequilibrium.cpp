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
	CMSSM model(nL, nTL, nS, state_equation_parameter, measurement_equation_parameter, transition_prob_parameter, new MakeABPsiPiC_ststm1, new MeasurementEquationFunction_test, new CTransitionProbMatrixFunction_Test); 

	size_t nFree = 23+8+2; // 23: model parameters; 6=2X3+2: sunspot parameters with 2 endogenous errors and 3 fundamental shocks; 2: probability of staying in ZLB
	/* Find valid starting value: begin */
	// Lower and upper bounds for x 
	TDenseVector lb, ub; 
	if ( MakeLbUb_ststm1(lb, ub, nFree) ) 	
	{
		cerr << "------MakeLbUb(): error occurred ------\n";
		abort(); 
	}

	TDenseVector x0(nFree), x0Valid(nFree), function_value;  
	size_t max_count = 500; 
	int return_code; 

	/* Find valid starting value for x */
	x0.RandomNormal(nFree); 			// Initial guess of x
	return_code = model.ValidInitialPoint(function_value, x0Valid, x0, max_count, lb, ub);
	if ( return_code < 0 )
	{
		cerr << "------CMSSM::ValidInitialPoint(): state equation function, measurement equation function or transition prob function not properly set ------.\n"; 
		abort();
	}
	else if ( return_code > 0)
	{
		cerr << "------ CMSSM::ValidInitialPoint(): no equilibrium exists ------.\n"; 
		abort(); 
	}
	// Display x0Valid 
	// Display(x0Valid); 
	/* Find valid start value: end */

	/* Starts in state one with probability one */
	TDenseVector initial_prob(model.nXi, 0.0);
	initial_prob.SetElement(1.0, 0); 
	 
	/* Initial values for initial Kalman filter */
	vector<TDenseVector> z0(model.nXi, TDenseVector(model.nZ,0.0) );
	vector<TDenseMatrix> P0(model.nXi, TDenseMatrix(model.nZ, model.nZ) ); 
	for (unsigned int i=0; i<P0.size(); i++)
		P0[i].Identity(model.nZ); 	

	/* Maximize log-likelihood */
	const double MINUS_INFINITY = -1.0e30; 
	size_t n_tries = 10; 
	TDenseVector best_solution(nFree+1, 0.0);
	best_solution.SetElement(MINUS_INFINITY, 0); 
	vector<TDenseVector> solutions(n_tries, TDenseVector(nFree+1,0.0) );  

	unsigned int i=0, number_bad=0; 
	bool bad; 
	double ll;	// log likelihood
	double mll;	// minus log likelihood
	TDenseVector xOptimal(nFree); 

	while (i < n_tries && number_bad < 200)
	{
		solutions[i].SetElement(MINUS_INFINITY, 0); 
		bad = true; 
		
		// Find valid starting value
		x0.RandomNormal(nFree); 
		if ( (x0[x0.dim-1] <=0.0) || (x0[x0.dim-1] > 1.0) )
			x0.SetElement(0.5,x0.dim-1);
		if ( (x0[x0.dim-2] <=0.0) || (x0[x0.dim-2] > 1.0) )
			x0.SetElement(0.5,x0.dim-2); 

		max_count = 10; 
		if ( MakeLbUb_ststm1(lb, ub, nFree) )
        	{
                	cerr << "------ MakeLbUB(): error occurred ------\n";
                	bad = true; 
        	}
		else 
		{
			return_code =  model.ValidInitialPoint(function_value, x0Valid, x0, max_count, lb, ub);
			if (return_code < 0)
			{
				cerr << "------ CMSSM::ValidInitialPoint(): state equation function, measurement equation function or transition prob function not properly set ------.\n";
                		bad = true; 
			}
			else if (return_code > 0)
        		{
                		cerr << "------ CMSSM::ValidInitialPoint(): no equilibrium exists ------" << endl;
                		bad = true; 
        		}
			else 
			{
				// unconstrained minimize minus loglikelihood
				return_code = model.Minimize_MinusLogLikelihood(mll,xOptimal,qdata,z0,P0,initial_prob,x0Valid);  
				if (return_code < 0)
				{
					cerr << "------ CMSSM::Minimize_MinusLogLikelihood(): state equation function, measurement equation function or transition prob function not properly set ------.\n";
					bad = true; 
				}	
				else if (return_code > 0)
				{
					cerr << "------ CMSSM::Minimize_MinusLogLikelihood(): npsol error " << return_code << endl; 
					bad = true; 
				}
				else 
				{
					ll = -mll; 
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
		}
		if (bad)
			number_bad ++;
	}
}
