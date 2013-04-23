#include <cmath>
#include "CMSSM.hpp"

using namespace std; 

CMSSM::CMSSM() :
nL(0), nTL(0), 
nS(0), nNu(0), nXi(0), nZeta(0), 
nZ(0), nY(0), nE(0), nU(0), 
state_equation_parameter(), 
measurement_equation_parameter(), 
transition_prob_parameter(),
rational_expectation_function(NULL), 
state_equation_function(NULL),
measurement_equation_function(NULL),
transition_prob_function(NULL),
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


CMSSM::CMSSM(size_t _nL, size_t _nTL, size_t _nS, const TDenseVector &_sfp, const TDenseVector &_mfp, const TDenseVector &_tpp, RationalExpectationFunction *_ref, StateEquationFunction *_sef, MeasurementEquationFunction *_mef, TransitionProbMatrixFunction *_tpmf, const TDenseVector &_cx) :
nL(_nL), nTL(_nTL),
nS(_nS), 
nZ(0), nY(0), nE(0), nU(0),
state_equation_parameter(_sfp), 
measurement_equation_parameter(_mfp), 
transition_prob_parameter(_tpp),
rational_expectation_function(_ref), 
state_equation_function(_sef),
measurement_equation_function(_mef),
transition_prob_function(_tpmf), 
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
nZ(right.nZ), nY(right.nY),
nE(right.nE), nU(right.nU),
A(right.A), B(right.B), C(right.C), Psi(right.Psi), Pi(right.Pi), 
a(right.a), H(right.H), Phi_u(right.Phi_u), R(right.R),
b(right.b), F(right.F), Phi_e(right.Phi_e), V(right.V),
state_equation_parameter(right.state_equation_parameter), 
measurement_equation_parameter(right.measurement_equation_parameter),
transition_prob_parameter(right.transition_prob_parameter),
rational_expectation_function(right.rational_expectation_function),
state_equation_function(right.state_equation_function), 
measurement_equation_function(right.measurement_equation_function),
transition_prob_function(right.transition_prob_function), 
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
	state_equation_parameter = right.state_equation_parameter; 
	measurement_equation_parameter = right.measurement_equation_parameter; 
	transition_prob_parameter = right.transition_prob_parameter; 
	rational_expectation_function = right.rational_expectation_function; 
	state_equation_function = right.state_equation_function; 
	measurement_equation_function = right.measurement_equation_function; 
	transition_prob_function = right.transition_prob_function; 
	current_x = right.current_x; 
        return *this;
}

bool CMSSM::CheckModelFunctions() const
// 
// Returns:
// 	true: if rational_expectation_function, state_equation_function or measurement_equation_function or transition_prob_function 
// 		not properly set 
// 	false: if all the above equations are properly set
{
	if (rational_expectation_function && state_equation_function && measurement_equation_function && transition_prob_function)
		return false; 
	else 
		return true; 
}

int CMSSM::UpdateStateModelParameters(unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x)
{
	if (!UpdateStateEquationParameter(t,y,x) && !UpdateMeasurementEquationParameter(t,y,x) ) 
	{
		current_x = x;
		return 0; 
	}
	else 
		return 1;  
}
