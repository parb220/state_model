#include "CMSSM.hpp"

using namespace std; 

void CMSSM::MeasurementEquation(MeasurementEquationFunction *function)
{
	a = vector<TDenseVector>(nNu, TDenseVector(nY,0.0) ); 
	H = vector<TDenseMatrix>(nNu, TDenseMatrix(nY,nZ,0.0) );
	Phi_u = vector<TDenseMatrix>(nNu, TDenseMatrix(nY,nU,0.0) ); 
	R = vector<TDenseMatrix>(nNu,TDenseMatrix(nU,nU,0.0) ); 
	function->convert(a, H, Phi_u, R, x, measurement_equation_parameter); 
}
