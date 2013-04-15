#include <cstdlib>
#include "findequilibrium.hpp"
#include "CMakeABPsiPiC_ststm1.hpp"
#include "CMeasurementEquationFunction_test.hpp"
#include "CTransitionMatrix_test.hpp"
#include "MakeLbUb_ststm1.hpp"

using namespace std; 

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cerr << "Usage: " << argc << "file-containing-the-data" << enedl; 
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

	// sizes
	size_t nS = 4; 
	size_t nL = 1; 
	size_t nTL = 1; 
	size_t nFree = 23+8+2; // 23: model parameters; 6=2X3+2: sunspot parameters with 2 endogenous errors and 3 fundamental shocks; 2: probability of staying in ZLB

	// Setting up the CMSSM model
	CMSSM model(nL, nTL, nS, state_equation_parameter, measurement_equation_parameter, transition_prob_parameter); 
	model.x = TDenseVector(nFree); 

	// Find valid starting value
	MakeABPsiPiC_ststm1 function_ABPsiPiC; 
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
	
	// Read in data
	string filename = string(argv[1]);	
	vector<double> qm_date;		// dates
	vector<TDenseVector> qdata; 	// variables
	if (LoadData(qm_date, qdata, filename))
	{
		cerr << "------ Error occurred when reading data ------\n"; 
		abort(); 
	}
	size_t nSample = qdata.size(); 	// total sample size
	size_t nY = qdata[0].dim;		// number of observables
}
