#include "CMSSM.hpp"

using namespace std; 

int CMSSM::SetMeasurementEquationParameter(const TDenseVector &input_x)
// Return value
// 	-1: measurement_equation_function not specified
// 	0: success
// 	1: fail 
{
	if (!measurement_equation_function)
		return -1; 

	if (measurement_equation_function->convert(a, H, Phi_u, R, measurement_equation_parameter, input_x))
		return 1;				 

	if (nY != a[0].dim)
		nY = a[0].dim;
	if (nU != Phi_u[0].rows) 
		nU = Phi_u[0].rows;
	return 0; 
}

/* Nothing to do yet */
void CMSSM::UpdateMeasurementEquationParameter(unsigned int t, const vector<TDenseVector> &y)
{
}
/* Nothing to do yet */
