#include <cmath>
#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"

using namespace std; 

CMSSM::CMSSM() :
nL(0), nTL(0), 
nS(0), nNu(0), nXi(0), nZeta(0), 
nZ(0), nY(0), nE(0), nU(0), nExpectation(0),  
rational_expectation_function(NULL), 
state_equation_function(NULL),
measurement_equation_function(NULL),
transition_prob_function(NULL),
prior_distr_function(NULL),
current_x() 
{
	A = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0) );   
        B = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0) );   
        C = vector<vector<TDenseVector> >(0,vector<TDenseVector>(0) ); 
        Psi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0) );   
        Pi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0) );
	
	a = vector<TDenseVector>(0);
        H = vector<TDenseMatrix>(0);
        Phi_u = vector<TDenseMatrix>(0);
        R = vector<TDenseMatrix>(0);
        b = vector<TDenseVector>(0);
        F = vector<TDenseMatrix>(0);
        Phi_e = vector<TDenseMatrix>(0);
        V = vector<TDenseMatrix>(0);
}


CMSSM::CMSSM(size_t _nL, size_t _nTL, size_t _nS, size_t _nZ, size_t _nY, size_t _nU, size_t _nE, size_t _nExpectation,  RationalExpectationFunction *_ref, StateEquationFunction *_sef, MeasurementEquationFunction *_mef, TransitionProbMatrixFunction *_tpmf, PriorDistributionFunction *_priordf,  const TDenseVector &_cx) :
nL(_nL), nTL(_nTL),
nS(_nS), 
nZ(_nZ), nY(_nY), nE(_nE), nU(_nU), nExpectation(_nExpectation), 
rational_expectation_function(_ref), 
state_equation_function(_sef),
measurement_equation_function(_mef),
transition_prob_function(_tpmf), 
prior_distr_function(_priordf),
current_x(_cx)
{
	nNu = (size_t)pow(nS, nL+1); 
	nXi = (size_t)pow(nS, nTL+1); 
	nZeta = (size_t)pow(nS, nTL); 
	
	A = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0) ); 
	B = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0) ); 
	C = vector<vector<TDenseVector> >(0,vector<TDenseVector>(0) ); 
	Psi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0) ); 
	Pi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0) );

        a = vector<TDenseVector>(nNu);
        H = vector<TDenseMatrix>(nNu);
        Phi_u = vector<TDenseMatrix>(nNu);
        R = vector<TDenseMatrix>(nNu);
        b = vector<TDenseVector>(nNu);
        F = vector<TDenseMatrix>(nNu);
        Phi_e = vector<TDenseMatrix>(nNu);
        V = vector<TDenseMatrix>(nNu);
}

CMSSM::CMSSM(const CMSSM &right) :
nL(right.nL), nTL(right.nTL),
nS(right.nS), nNu(right.nNu), nXi(right.nXi), nZeta(right.nZeta),
nZ(right.nZ), nY(right.nY), nE(right.nE), nU(right.nU), nExpectation(right.nExpectation), 
A(right.A), B(right.B), C(right.C), Psi(right.Psi), Pi(right.Pi), 
a(right.a), H(right.H), Phi_u(right.Phi_u), R(right.R),
b(right.b), F(right.F), Phi_e(right.Phi_e), V(right.V),
rational_expectation_function(right.rational_expectation_function),
state_equation_function(right.state_equation_function), 
measurement_equation_function(right.measurement_equation_function),
transition_prob_function(right.transition_prob_function), 
prior_distr_function(right.prior_distr_function),
current_x(right.current_x)
{}

CMSSM & CMSSM::operator=(const CMSSM &right)
{
        nL = right.nL;
        nTL = right.nTL;
        nS = right.nS;
        nNu = right.nNu;
        nXi = right.nXi;
        nZeta = right.nZeta;
        nZ = right.nZ;
        nY = right.nY;
        nE = right.nE;
        nU = right.nU;
	nExpectation = right.nExpectation; 
	A = right.A; 
	B = right.B; 
	C = right.C; 
	Psi = right.Psi; 
	Pi = right.Pi; 
        a = right.a;
        H = right.H;
        Phi_u = right.Phi_u;
        R = right.R;
        b = right.b;
        F = right.F;
        Phi_e = right.Phi_e;
        V = right.V;
	rational_expectation_function = right.rational_expectation_function; 
	state_equation_function = right.state_equation_function; 
	measurement_equation_function = right.measurement_equation_function; 
	transition_prob_function = right.transition_prob_function; 
	prior_distr_function = right.prior_distr_function; 
	current_x = right.current_x; 
        return *this;
}

int CMSSM::CheckModelFunctions() const
// 
// Returns:
// 	ERROR_OCCURRED: if rational_expectation_function, state_equation_function or measurement_equation_function or transition_prob_function 
// 		not properly set 
// 	SUCCESS: if all the above equations are properly set
{
	if (rational_expectation_function && state_equation_function && measurement_equation_function && transition_prob_function)
		return SUCCESS; 
	else 
		return ERROR_OCCURRED; 
}

int CMSSM::UpdateStateModelParameters(unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x)
{
	if (UpdateStateEquationParameter(t,y,x) == SUCCESS && UpdateMeasurementEquationParameter(t,y,x) == SUCCESS ) 
	{
		current_x = x;
		return SUCCESS; 
	}
	else 
		return ERROR_OCCURRED;  
}
