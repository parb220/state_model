#include "CMSSM.hpp"

using namespace std; 

void CMSSM::MeasurementEquation(MeasurementEquationFunction *function)
{
	if (function->convert(a, H, Phi_u, R, x, measurement_equation_parameter))
	{
		nY = a[0].dim; 
		nU = Phi_u[0].rows;
	}
}
