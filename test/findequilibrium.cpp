#include "CMSSM.hpp"

using namespace std; 

int main()
{
	// % prefixed parameters
	TDenseVector state_equation_parameter(5);
	state_equation_parameter.SetElement(1.005,0);	// glambdaz, Gross rate: (1+2%) annually per capita or (1+0.5%) quarterly per capita 
	state_equation_parameter.SetElement(0.01,1);	// rn, Net rate: rn=4% annually or 1% quarterly (natural rate of interest or steady state real interest rate)
	state_equation_parameter.SetElement(0.005,2);	// gpistar, Net rate: log(1.02) -- 2% annually or 0.5% quarterly
	state_equation_parameter.SetElement(0.7,3);	// sc, Steady state share of private consumption in C+G
	state_equation_parameter.SetElement(0.00025,4);	// Rlow, Net rate 10 basis points (0.1%) interest rate annually at zero bound or 0.025% quarterly at zero bound. In 2012Q2, it is 3.83e04*4-0.0015 annually

	
}
