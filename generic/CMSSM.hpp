// Because the copy constructor and overloaded = operator for both TDenseVector and
// TDensematrix are implemented as equating pointers plus incrementing the memory 
// count (rather than newing/copying memory), the constructor function for CMSSM 
// generates vectors of TDenseVector of TDenseMatrix for (a, H, V, b, F, Phi) where
// each element, which is either a TDenseVector or a TDenseMatrix object, is copy-
// constructed from the corresponding element of $right. Therefore the resulting 
// CMSSM object does not completely own its own memory: (1) s, nState, nMeasurement,
// nSNoise and nMNoise are of its own, while (2) a, H, V, b, F, and Phi are pointing
// to  the old memory pointed by $right   

#ifndef _CLASS_CMSSM
#define _CLASS_CMSSM 

#include <vector>
#include "dw_dense_matrix.hpp"

using namespace std; 

class MakeABPsiPiC; 
class TransitionProbMatrixFunction; 
class MeasurementEquationFunction; 
class CMSSM; 
class ObjectiveFunction_Validation; 

class MakeABPsiPiC 
// Derive gensys form of the state equation from x and makeAB_inputs
// return true if some A is not of full rank and false otherwise
//
// A(s(t)*x(t) = C(s(t)) + B(s(t))*x(t-1) + Psi(s(t))*epsilon(t) + Pi(s(t))*eta(t)
//
// A, B, Psi and Pi are matrices (vector<vector<> >) of TDenseMatrix
// C is an matrix (vector<vector<> >) of TDenseVector
{
public:
	virtual bool convert(vector<vector<TDenseMatrix> > &A, vector<vector<TDenseMatrix> > &B, vector<vector<TDenseMatrix> > &Psi, vector<vector<TDenseMatrix> >&Pi, vector<vector<TDenseVector> > &C, const TDenseVector &x, const TDenseVector &free_parameter)=0; 
};

class MeasurementEquationFunction
// Measurement equation y(t) = a(s(t)) + H(s(t))*z(t) +phi_u(s(t))*u(t)
// a, H, Phi_u, x, nY, nZ, nU, nNU correspond to parameters in CSSM model
// free_parameter contains any parameters (e.g., gpistar) that can be passed 
{
public:
	virtual bool convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &x, const TDenseVector &free_parameter)=0; 
}; 


class TransitionProbMatrixFunction
// Transition probability derived from a vector and measurment variables
// s: number of regimes (transition probability matrix is s-by-s)
// Y: measurement variables corresponding to one/multiple time points
// x: parameters of the original model from which CMSSM is derived 
{
public:
	virtual void convert(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x, size_t nS, size_t nTL, const TDenseVector &free_parameter) = 0; 
};

//== Reigme (Markov) swithcing state model ==
class CMSSM
{
public:
	// Parameters
	size_t nL;	// number of lagged regimes that enter the measurement or state equations
	size_t nTL;	// total number of lagged regimes to track
	size_t nS;	// number of regimes
	size_t nNu;	// == nS^lags, number of regimes affecting coefficients
	size_t nXi;	// == nS^(total_lags+1), total number of regimes being tracked
	size_t nZeta;	// == nS^total_lags; 
	size_t nZ;	// number of state variables
	size_t nY;	// number of measurement variables
	size_t nU;	// number of noise variables in measurement
	size_t nE; 	// number of noise variables in state 

	TDenseVector x;	// 	parameters of the original model from which CMSSM is derived

	vector<vector<TDenseMatrix> >A; //
	vector<vector<TDenseMatrix> >B; // parameters of the gensys form of the state equation
	vector<vector<TDenseVector> >C; // A(s(t)) * x(t) = C(s(t)) + B(s(t))*x(t-1) + ... 
	vector<vector<TDenseMatrix> >Psi; //                Psi(s(t))*epsilon(t) + Pi(s(t)*eta(t) 
	vector<vector<TDenseMatrix> >Pi; //
	
	vector <TDenseVector> a;	// 
	vector <TDenseMatrix> H; 	//	y_t=a_{st} + H{st}z_t + Phi_u_{st}u_t
	vector <TDenseMatrix> Phi_u; 	//
	vector <TDenseMatrix> R;	// 	R = Phi_u*Phi_u'
	vector <TDenseVector> b; 	//
	vector <TDenseMatrix> F; 	//	z_t=b_{st} + F_{st}z_{t-1} + Phi_e_{st}e_t
	vector <TDenseMatrix> Phi_e; 	//
	vector <TDenseMatrix> V;	//	V = Phi_e*Phi_e'

	TDenseMatrix Q;	// Transition matrix;

	TDenseVector state_equation_parameter;	// free parameters to be used in state equations
	TDenseVector measurement_equation_parameter;	// free parameters to be used in measurement equations
	TDenseVector transition_prob_parameter;	// free parameters to be used in transition prob

	// Parameter specification for the state equation 
	void ClearStateEquation(); 
	unsigned int StateEquationHelper(MakeABPsiPiC *); 	
	void StateEquation(unsigned int t, const vector<TDenseVector> &y);

	// Parameter specification for the measurement equation
	void MeasurementEquation(MeasurementEquationFunction *function);

	// Valid initial point
	bool ValidInitialPoint(TDenseVector &, const TDenseVector &x0, size_t max_count, TDenseVector& lb, TDenseVector &ub); 

	// Check solution
	bool CheckSolution(MakeABPsiPiC *); 

	// Constructor and destrunctor 
	CMSSM();
	CMSSM(size_t _nL, size_t _nTL, size_t _nS, const TDenseVector &, const TDenseVector &, const TDenseVector &); 
	CMSSM(const CMSSM &right); 
	CMSSM &operator=(const CMSSM &right); 
	~CMSSM() {}

	// Kalman filter
	bool KalmanFilter(double &log_likelihood, vector<vector<TDenseVector> > &z_tm1, vector<vector<TDenseMatrix> > &P_tm1, vector<TDenseVector > &p_tm1, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, TransitionProbMatrixFunction *transition_matrix_function, const TDenseVector &initial_prob);
	double minus_log_likelihood(const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, TransitionProbMatrixFunction *transition_matrix_function, const TDenseVector &initial_prob);

}; 

class ObjectiveFunction_Validation
{
public:
	static MakeABPsiPiC *ABPsiPiC_function; 
	static TDenseVector state_equation_parameter; 
	static TDenseVector measurement_equation_parameter; 
	static TDenseVector transition_prob_parameter;
	static void *function(int *mode, int *n, double *x, double *f, double *g, int *);  
};

#endif
