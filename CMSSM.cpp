#include <cmath>
#include "CMSSM.hpp"

using namespace std; 

CMSSM::CMSSM() :
nL(0), nTL(0), 
nS(0), nNu(0), nXi(0), nZeta(0), 
nZ(0), nY(0), nE(0), nU(0), 
x(), 
state_equation_parameter(), 
measurement_equation_parameter(), 
transition_prob_parameter()
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

        Q = TDenseMatrix();
}


CMSSM::CMSSM(size_t _nL, size_t _nTL, size_t _nS, size_t _nZ, size_t _nY, size_t _nE, size_t _nU, const TDenseVector &_x, const TDenseVector &_sfp, const TDenseVector &_mfp, const TDenseVector &_tpp) :
nL(_nL), nTL(_nTL),
nS(_nS), nZ(_nZ), nY(_nY), nE(_nE), nU(_nU),
x(_x),
state_equation_parameter(_sfp), 
measurement_equation_parameter(_mfp), 
transition_prob_parameter(_tpp) 
{
	nNu = (size_t)pow(nS, nL); 
	nXi = (size_t)pow(nS, nTL+1); 
	nZeta = (size_t)pow(nS, nTL); 
	
	A = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) ); 
	B = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) ); 
	C = vector<vector<TDenseVector> >(0,vector<TDenseVector>(0,TDenseVector() ) ); 
	Psi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) ); 
	Pi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix() ) );

        a = vector<TDenseVector>(nNu, TDenseVector(nY,0.0));
        H = vector<TDenseMatrix>(nNu, TDenseMatrix(nY,nZ,0.0));
        Phi_u = vector<TDenseMatrix>(nNu, TDenseMatrix(nY,nU,0.0));
        R = vector<TDenseMatrix>(nNu, TDenseMatrix(nY,nY,0.0));
        b = vector<TDenseVector>(nNu, TDenseVector(nZ,0.0));
        F = vector<TDenseMatrix>(nNu, TDenseMatrix(nZ,nZ,0.0));
        Phi_e = vector<TDenseMatrix>(nNu, TDenseMatrix(nZ,nE,0.0));
        V = vector<TDenseMatrix>(nNu, TDenseMatrix(nZ,nZ,0.0));

	Q = TDenseMatrix(nXi,nXi,0.0); 
}

CMSSM::CMSSM(const CMSSM &right) :
nL(right.nL), nTL(right.nTL),
nS(right.nS), nNu(right.nNu), nXi(right.nXi), nZeta(right.nZeta),
nZ(right.nZ), nY(right.nY),
nE(right.nE), nU(right.nU),
x(right.x),
A(right.A), B(right.B), C(right.C), Psi(right.Psi), Pi(right.Pi), 
a(right.a), H(right.H), Phi_u(right.Phi_u), R(right.R),
b(right.b), F(right.F), Phi_e(right.Phi_e), V(right.V),
Q(right.Q),
state_equation_parameter(right.state_equation_parameter), 
measurement_equation_parameter(right.measurement_equation_parameter),
transition_prob_parameter(right.transition_prob_parameter) 
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
        x = right.x;
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
	Q = right.Q; 
	state_equation_parameter = right.state_equation_parameter; 
	measurement_equation_parameter = right.measurement_equation_parameter; 
	transition_prob_parameter = right.transition_prob_parameter; 
        return *this;
}

