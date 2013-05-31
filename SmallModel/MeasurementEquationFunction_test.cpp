#include "MeasurementEquationFunction_test.hpp"
#include "CMSSM_Error_Code.hpp"

using namespace std; 

int MeasurementEquationFunction_test::convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &x, size_t nZ, size_t nY, size_t nU, size_t nNu)
{
	// Proper sizing
	a = vector<TDenseVector>(nNu); 
	H = vector<TDenseMatrix>(nNu); 
	Phi_u = vector<TDenseMatrix>(nNu); 
	R = vector<TDenseMatrix>(nNu); 
	
	ConvertXtoParameter(x); 
	
	// ameasure and Hmeasure
	TDenseVector ameasure(nY, 0.0); 
	TDenseMatrix Hmeasure(nY, nZ, 0.0); 

	// #include "Regime1/MathematicaInputMeasureMatrices.cpp"
	#include "MathematicaInputMeasureMatrices.Regime1"
	// Because measurement equations are the same for all the regimes
	// we only use Regime1's to assign a, H, Phi_u and R for all the regimes
	for (unsigned int i=0; i<nNu; i++)
	{
		a[i].CopyContent(ameasure); 
		H[i].CopyContent(Hmeasure); 
	}

	// #include "Regime2/MathematicaInputMeasureMatrices.cpp"
	return SUCCESS; 
}
