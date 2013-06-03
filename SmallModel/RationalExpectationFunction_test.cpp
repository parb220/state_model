#include "RationalExpectationFunction_test.hpp"
#include "CMSSM_Error_Code.hpp"

using namespace std; 

int RationalExpectationFunction_test::convert(vector<vector<TDenseMatrix> > &A_output, vector<vector<TDenseMatrix> > &B_output, vector<vector<TDenseMatrix> > &Psi_output, vector<vector<TDenseMatrix> >&Pi_output, vector<vector<TDenseVector> > &C_output, const TDenseVector &x, size_t nZ, size_t nY, size_t nU, size_t nE, size_t nExpectation)
{
	int error_code = SUCCESS; 
	// Proper sizing
	size_t nRegime = 2; 
	A_output = vector<vector<TDenseMatrix> >(nRegime,vector<TDenseMatrix>(nRegime)); 
	B_output = vector<vector<TDenseMatrix> >(nRegime,vector<TDenseMatrix>(nRegime)); 
	Psi_output = vector<vector<TDenseMatrix> >(nRegime,vector<TDenseMatrix>(nRegime)); 
	Pi_output = vector<vector<TDenseMatrix> >(nRegime,vector<TDenseMatrix>(nRegime)); 
	C_output = vector<vector<TDenseVector> >(nRegime,vector<TDenseVector>(nRegime));  

	// Convert x to meaningful parameters
	ConvertXtoParameter(x); 
	#include "Regime1/MathematicaInputDerivedSSPars.cexps"

	// Regime 1; 
	TDenseMatrix A(nZ,nZ,0.0), B(nZ,nZ,0.0), CapitalGammau(nZ,nE,0.0), CapitalGammaeta(nZ,nExpectation,0.0); 
	TDenseVector cgensys(nZ,0.0); 
	// #include "Regime1/MathematicaInputGensysForm.cpp"
	#include "MathematicaInputGensysForm.Regime1"
	A_output[0][0].CopyContent(A); 
	A_output[0][1].CopyContent(A); 
	if (Rank(A_output[0][0]) < nZ || Rank(A_output[0][1]) < nZ)
		error_code = ERROR_OCCURRED; 
	B_output[0][0].CopyContent(B); 
	B_output[0][1].CopyContent(B); 
	Psi_output[0][0].CopyContent(CapitalGammau); 
	Psi_output[0][1].CopyContent(CapitalGammau); 
	Pi_output[0][0].CopyContent(CapitalGammaeta); 
	Pi_output[0][1].CopyContent(CapitalGammaeta); 
	C_output[0][0].CopyContent(cgensys); 
	C_output[0][1].CopyContent(cgensys); 

	// Regime 2; 
	//#include "Regime2/MathematicaInputGensysForm.cpp"
	A.Zeros(nZ,nZ); 
	B.Zeros(nZ,nZ); 
	CapitalGammau.Zeros(nZ,nE); 
	CapitalGammaeta.Zeros(nZ,nExpectation); 
	#include "MathematicaInputGensysForm.Regime2"
	A_output[1][0].CopyContent(A); 
	A_output[1][1].CopyContent(A); 
	if (Rank(A_output[1][0]) < nZ || Rank(A_output[1][1]) < nZ)
		error_code = ERROR_OCCURRED; 
	B_output[1][0].CopyContent(B); 
	B_output[1][1].CopyContent(B); 
	Psi_output[1][0].CopyContent(CapitalGammau); 
	Psi_output[1][1].CopyContent(CapitalGammau); 
	Pi_output[1][0].CopyContent(CapitalGammaeta); 
	Pi_output[1][1].CopyContent(CapitalGammaeta); 
	C_output[1][0].CopyContent(cgensys); 
	C_output[1][1].CopyContent(cgensys); 
	
	return error_code; 
} 
