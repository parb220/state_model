#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"

using namespace std; 

void CMSSM::ClearMeasurementEquation()
{
	a = vector<TDenseVector>(0); 
	H = vector<TDenseMatrix>(0); 
	Phi_u = vector<TDenseMatrix>(0); 
	R = vector<TDenseMatrix>(0); 
}

int CMSSM::UpdateMeasurementEquationParameter(unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x)
// Returns: 
// 	-1: measurement_equation_function not properly set
// 	0: success
// 	>0: error code returned by measurement_equation_function 
{
	if (current_x.dim != x.dim || !(current_x == x) )
	{
		if (measurement_equation_function == NULL)
		{
			// cerr << "MeasurementEquationFunction is not properly set up.\n"; 
			ClearMeasurementEquation(); 
			return MODEL_NOT_PROPERLY_SET; 
		}
		int error_code = measurement_equation_function->convert(a, H, Phi_u, R, x);
		if (error_code != SUCCESS)
		{
			// cerr << "Error occurred during MeasurementEquationFunction call: " << error_code << endl; 
			ClearMeasurementEquation(); 
			return error_code; 
		}
		if (nY != a[0].dim)
			nY = a[0].dim;
		if (nU != Phi_u[0].rows) 
			nU = Phi_u[0].rows;
	}
	return SUCCESS; 
}
