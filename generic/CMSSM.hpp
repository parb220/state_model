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
class PriorDistributionFunction; 
class RationalExpectationFunction;
class StateEquationFunction;  
class TransitionProbMatrixFunction; 
class MeasurementEquationFunction; 
class CMSSM; 
class MinusLogLikelihood;
class MinusLogPosterior;
class ObjectiveFunction_Validation;

class PriorDistributionFunction
{
public:
	virtual double log_pdf(const TDenseVector &x, const TDenseVector &fixed_parameter) = 0; 
};

class RationalExpectationFunction
// Gensys form of the state equation
//
// A(s(t)*x(t) = C(s(t)) + B(s(t))*x(t-1) + Psi(s(t))*epsilon(t) + Pi(s(t))*eta(t)
//
// Assume:
// 	fixed_parameter
// 	x:	free parameter
//
// Result:
// 	A, B, Psi, Pi and C are specified
//
// Return:
// 	1: some A is not of full rank
// 	0: success
{
public:
	virtual int convert(vector<vector<TDenseMatrix> > &A, vector<vector<TDenseMatrix> > &B, vector<vector<TDenseMatrix> > &Psi, vector<vector<TDenseMatrix> >&Pi, vector<vector<TDenseVector> > &C, const TDenseVector &fixed_parameter, const TDenseVector &x)=0; 
};

class StateEquationFunction
// State equation z(t) = b(s(t)) + F(s(t))*z(t-1) + Phi_e(s(t))*e(t)
//
// Assume:
// 	A, B, Psi, Pi, and C: CMSSM parameters resulted from RationalExpectationFunction
//	fixed_parameter:
//	x: free_parameter
// 
// Result:
// 	b, F, Phi_e and V are specified
//
// Return:
// 	0: success 
// 	1: otherwise
{
public:
	virtual int convert(vector<TDenseVector> &b, vector<TDenseMatrix> &F, vector<TDenseMatrix> &Phi_e, vector<TDenseMatrix> &V, const vector<vector<TDenseMatrix> > &A, const vector<vector<TDenseMatrix> > &B, const vector<vector<TDenseMatrix> > &Psi, const vector<vector<TDenseMatrix> >&Pi, const vector<vector<TDenseVector> > &C, const TDenseVector &fixed_parameter, const TDenseVector &x, size_t, size_t, size_t)=0;
};

class MeasurementEquationFunction
// Measurement equation y(t) = a(s(t)) + H(s(t))*z(t) + Phi_u(s(t))*u(t)
//
// Assume:
// 	fixed_parameter:
// 	x: free_parameter
//
// Result:	
// 	a, H, Phi_u and R are specfied
// 
// Return:
// 	0: success 
// 	1: otherwise 
{
public:
	virtual int convert(vector<TDenseVector> &a, vector<TDenseMatrix> &H, vector<TDenseMatrix> &Phi_u, vector<TDenseMatrix> &R, const TDenseVector &fixed_parameter, const TDenseVector &x)=0; 
}; 

class TransitionProbMatrixFunction
// Transition probability 
// 
// Assume:
// 	t: time
// 	y: observations
// 	nS: number of regimes
// 	nTL: total number of lagged regimes 
// 	fixed_parameter 
// 	x: free parameter
// 
// Result
// 	Q: transition probability matrix
//
// Return: 
// 	SUCCESS: success
// 	ERROR_OCCURRED: error_occurred
{
public:
	virtual int convert(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, size_t nS, size_t nTL, const TDenseVector &fixed_parameter, const TDenseVector &x) = 0; 
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

protected:
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

protected:
	TDenseVector current_x; 	// x that have just been used to setup 
					// state equation, measurement equation
					// and transition probability matrix

public:
	// fixed parameters
	TDenseVector state_equation_parameter;	// parameters to be used in state equations
	TDenseVector measurement_equation_parameter;	// parameters to be used in measurement equations
	TDenseVector transition_prob_parameter;	// parameters to be used in transition prob
	TDenseVector prior_distr_parameter; 	// parameters to be used in specifying prior distribution functions 

	// functions used to specify rationa expectation model, state model, and transition prob matrix
	RationalExpectationFunction *rational_expectation_function;
	StateEquationFunction *state_equation_function; 
	MeasurementEquationFunction *measurement_equation_function; 
	TransitionProbMatrixFunction *transition_prob_function; 
	PriorDistributionFunction *prior_distr_function; 
	
	int UpdateStateModelParameters(unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x); 
	int GetTranstionProbMatrix(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x) const; 	

	// Valid initial point
	int ValidInitialPoint(TDenseVector &, TDenseVector &, const TDenseVector &x, size_t max_count, const TDenseVector& lb, const TDenseVector &ub); 

	// Kalman filter
	// 	Because it calls UpdateStateModelParameters, it cannot be constant
	int KalmanFilter(double &log_likelihood, vector<vector<TDenseVector> > &z_tm1, vector<vector<TDenseMatrix> > &P_tm1, vector<TDenseVector > &p_tm1, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob,  const TDenseVector &x);

	// Calculate log likelihood
	// 	Because it calls KalmanFilter, it cannot be constatn
	int LogLikelihood(double &log_likelihood, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob); 
	
	// Minimize minus log likelihood 
	// 	Because it calls KalmanFilter, it cannot be constant
	int Maximize_LogLikelihood(double &minus_log_likelihood_optimal, TDenseVector &x_optimal, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob, const TDenseVector &x0);	
 
	// Calculate log posterior 
	// 	Because it calls LogLikelihood (which calls KalmanFilter), it cannot be constant
	int LogPosterior(double &log_posterior, const TDenseVector &x, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob); 

	// Maximize log posterior
	// 	Because it calls KalmanFilter, it cannot be constant
	int Maximize_LogPosterior(double &log_posterior_optimal, TDenseVector &x_optimal, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob, const TDenseVector &x0);

	// Constructor and destrunctor 
	CMSSM();
	CMSSM(size_t _nL, size_t _nTL, size_t _nS, const TDenseVector &, const TDenseVector &, const TDenseVector &, const TDenseVector &, RationalExpectationFunction * =NULL, StateEquationFunction * =NULL, MeasurementEquationFunction * =NULL, TransitionProbMatrixFunction * =NULL, PriorDistributionFunction * = NULL, const TDenseVector & =TDenseVector()); 
	CMSSM(const CMSSM &right); 
	CMSSM &operator=(const CMSSM &right); 
	~CMSSM() {}

protected:
	int UpdateStateEquationParameter(unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x); 
	int UpdateMeasurementEquationParameter(unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x);
	void ClearStateEquation(); 
	void ClearMeasurementEquation(); 
	int CheckModelFunctions() const; 
}; 

class MinusLogLikelihood
{
public:
        static CMSSM *model;
        static vector<TDenseVector> y;
        static vector<TDenseVector> z_0;
        static vector<TDenseMatrix> P_0;
        static TDenseVector initial_prob;

        static void *function(int *mode, int*n, double *x, double *f, double *g, int *nstate);
};

class MinusLogPosterior
{
public:
	static CMSSM *model; 
	static vector<TDenseVector> y; 
	static vector<TDenseVector> z_0; 
	static vector<TDenseMatrix> P_0; 
	static TDenseVector initial_prob; 

	static void *function(int *mode, int *n, double *x, double *f, double *g, int *nstate); 
};

class ObjectiveFunction_Validation
{
public:
        static CMSSM *model;
        static void *function(int *mode, int *n, double *x, double *f, double *g, int *nstate);
};

#endif
