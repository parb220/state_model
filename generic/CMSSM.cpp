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
state_equation_function(NULL),
measurement_equation_function(NULL),
transition_prob_function(NULL)
{
	A = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) );   
        B = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) );   
        C = vector<vector<TDenseVector> >(0,vector<TDenseVector>(0,TDenseVector() ) ); 
        Psi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) );   
        Pi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) );
	
	a = vector<TDenseVector>(0, TDenseVector());
        H = vector<TDenseMatrix>(0, TDenseMatrix());
        Phi_u = vector<TDenseMatrix>(0, TDenseMatrix());
        R = vector<TDenseMatrix>(0, TDenseMatrix());
        b = vector<TDenseVector>(0, TDenseVector());
        F = vector<TDenseMatrix>(0, TDenseMatrix());
        Phi_e = vector<TDenseMatrix>(0, TDenseMatrix());
        V = vector<TDenseMatrix>(0, TDenseMatrix());
}


CMSSM::CMSSM(size_t _nL, size_t _nTL, size_t _nS, const TDenseVector &_sfp, const TDenseVector &_mfp, const TDenseVector &_tpp, MakeABPsiPiC *_sef, MeasurementEquationFunction *_mef, TransitionProbMatrixFunction *_tpmf) :
nL(_nL), nTL(_nTL),
nS(_nS), nZ(0), nY(0), nE(0), nU(0),
state_equation_parameter(_sfp), 
measurement_equation_parameter(_mfp), 
transition_prob_parameter(_tpp), 
state_equation_function(_sef),
measurement_equation_function(_mef),
transition_prob_function(_tpmf)
{
	nNu = (size_t)pow(nS, nL+1); 
	nXi = (size_t)pow(nS, nTL+1); 
	nZeta = (size_t)pow(nS, nTL); 
	
	A = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) ); 
	B = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) ); 
	C = vector<vector<TDenseVector> >(0,vector<TDenseVector>(0,TDenseVector() ) ); 
	Psi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) ); 
	Pi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) );

        a = vector<TDenseVector>(nNu, TDenseVector());
        H = vector<TDenseMatrix>(nNu, TDenseMatrix());
        Phi_u = vector<TDenseMatrix>(nNu, TDenseMatrix());
        R = vector<TDenseMatrix>(nNu, TDenseMatrix());
        b = vector<TDenseVector>(nNu, TDenseVector());
        F = vector<TDenseMatrix>(nNu, TDenseMatrix());
        Phi_e = vector<TDenseMatrix>(nNu, TDenseMatrix());
        V = vector<TDenseMatrix>(nNu, TDenseMatrix());
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
state_equation_function(right.state_equation_function), 
measurement_equation_function(right.measurement_equation_function),
transition_prob_function(right.transition_prob_function)
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
	state_equation_function = right.state_equation_function; 
	measurement_equation_function = right.measurement_equation_function; 
	transition_prob_function = right.transition_prob_function; 
        return *this;
}

bool CMSSM::CheckStateMeasurementTransitionEquations() const
// 
// Returns:
// 	true: if state_equation_function or measurement_equation_function or transition_prob_function 
// 		not properly set 
// 	false: if all the above equations are properly set
{
	if (state_equation_function && measurement_equation_function && transition_prob_function)
		return false; 
	else 
		return true; 
}
