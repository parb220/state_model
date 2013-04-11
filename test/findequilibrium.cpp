#include "findequilibrium.hpp"
#include "CMakeABPsiPiC_ststm1.hpp"
#include "CMeasurementEquationFunction_test.hpp"
#include "CTransitionMatrix_test.hpp"
#include "MakeLbUb_ststm1.hpp"

using namespace std; 

int main()
{
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

	// Loading data; 

	// sizes
	size_t nS = 4; 
	size_t nL = 1; 
	size_t nTL = 1; 
	size_t nZ = 7; 
	size_t nU = 0; 
	

	// CMSSM::CMSSM(size_t _nL, size_t _nTL, size_t _nS, size_t _nZ, size_t _nY, size_t _nE, size_t _nU, const TDenseVector &_x, const TDenseVector &_sfp, const TDenseVector &_mfp, const TDenseVector &_tpp)

	MakeABPsiPiC_ststm1 function_ABPsiPiC; 
	ObjectiveFunction_Validation::ABPsiPiC_function = &function_ABPsiPiC; 
	ObjectiveFunction_Validation::state_equation_parameter = state_equation_parameter; 

	// Omit checking if the model makes sense under 1st regime (Line 36 -- 42)
	
	// Find valid starting value
	// Lower and upper bounds 
	TDenseVector lb, ub; MakeLbUb_ststm1(lb, ub); 	
	size_t nFree = 23+6+1; // 23: model parameters; 6=2X3: sunspot parameters with 2 endogenous errors and 3 fundamental shocks; 1: probability of staying in ZLB
	TDenseVector x0(nFree); 	
	CMSSM model( 
	if (	

}
