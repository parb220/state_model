#include <cmath>
#include "CMSSM.hpp"

using namespace std; 

CMSSM::CMSSM(size_t _nL, size_t _nTL, size_t _nS, size_t _nZ, size_t _nY, size_t _nE, size_t _nU, const TDenseVector &_x) :
nL(_nL), nTL(_nTL),
nS(_nS), nZ(_nZ), nY(_nY), nE(_nE), nU(_nU),
x(_x)
{
	nNu = (size_t)pow(nS, nL); 
	nXi = (size_t)pow(nS, nTL+1); 
	nZeta = (size_t)pow(nS, nTL); 
	
	A = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix(0,0) ) ); 
	B = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix(0,0) ) ); 
	C = vector<vector<TDenseVector> >(0,vector<TDenseVector>(0,TDenseVector(0) ) ); 
	Psi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix(0,0) ) ); 
	Pi = vector<vector<TDenseMatrix> >(0,vector<TDenseMatrix>(0,TDenseMatrix(0,0) ) );

        a = vector<TDenseVector>(nNu, TDenseVector(nY,0.0));
        H = vector<TDenseMatrix>(nNu, TDenseMatrix(nY,nZ,0.0));
        Phi_u = vector<TDenseMatrix>(nNu, TDenseMatrix(nY,nU,0.0));
        R = vector<TDenseMatrix>(nNu, TDenseMatrix(nY,nY,0.0));
        b = vector<TDenseVector>(nNu, TDenseVector(nZ,0.0));
        F = vector<TDenseMatrix>(nNu, TDenseMatrix(nZ,nZ,0.0));
        Phi_e = vector<TDenseMatrix>(nNu, TDenseMatrix(nZ,nE,0.0));
        V = vector<TDenseMatrix>(nNu, TDenseMatrix(nZ,nZ,0.0));
}

CMSSM::CMSSM(const CMSSM &right) :
nL(right.nL), nTL(right.nTL),
nS(right.nS), nNu(right.nNu), nXi(right.nXi), nZeta(right.nZeta),
nZ(right.nZ), nY(right.nY),
nE(right.nE), nU(right.nU),
x(right.x),
A(right.A), B(right.B), C(right.C), Psi(right.Psi), Pi(right.Pi), 
a(right.a), H(right.H), Phi_u(right.Phi_u), R(right.R),
b(right.b), F(right.F), Phi_e(right.Phi_e), V(right.V)
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
        return *this;
}

