#include <cmath>
#include "CMSSM.hpp"

using namespace std;

bool CMSSM::CheckSolution(MakeABPsiPiC *function)
{
	int error_code = StateEquationHelper(function); 
	if (error_code)
	{
		cerr << "Unable to solve system.\n"; 
		return true; 
	}

	// function->convert(A,B,Psi,Pi,C,x,state_equation_parameter) has 
	// been called in StateEquationHelper
	// It is therefore skipped here
	vector<double> p(2,0.0); 
	p[0] = x[x.dim-2]; 
	p[1] = x[x.dim-1]; 

	// Size
	size_t nOriginalRegime = 2;  

	// Stable manifold and annihilator regime
	vector<TDenseMatrix> Vs(nOriginalRegime);	// Schur vectors corresponding to the stable (||<1) eigen values
	vector<TDenseMatrix> Va(nOriginalRegime);  	// Annihilator vectors
	for (unsigned int i=0; i<nOriginalRegime; i++)
	{
		TDenseMatrix T, U;	// T: schur form of X; U: schur vector
		TDenseVector eR, eI; 	// eR and eI: real and imaginary parts of the eigen values
		TDenseMatrix XSchur = LeftSolve(A[i][i],B[i][i]); 
		Schur(T, eR, eI, U, XSchur, true); 
		
		int *select = new int[XSchur.rows]; 
		unsigned int nStable = 0; 
		for (unsigned int j=0; j<eR.dim; j++)
		{
			if (sqrt(eR[j]*eR[j] + eI[j]*eI[j]) < 1.0)
			{
				nStable ++; 
				select[j] =1; 
			}
			else 
				select[j] = 0; 
		}

		if (nStable == nZ)
		{
			Vs[i].Identity(nZ); 
			Va[i].Initialize(0.0, 0, nZ); 
		}
		else 
		{
			TDenseMatrix OrderT, OrderU; 
			TDenseVector OrderER, OrderEI; 
			OrderSchur(OrderT, OrderER, OrderEI, OrderU, T, U, select, true); 
			Vs[i].Resize(nZ, nStable); 
			Va[i].Resize(nZ, nZ-nStable); 
			for (unsigned int k=0; k<nZ; k++)
			{
				for (unsigned int l=0; l<nStable; l++)
					Vs[i].SetElement(OrderU(k,l), k, l); 
				for (unsigned int l=nStable; l<nZ; l++)
					Va[i].SetElement(OrderU(k,l), k, l-nStable); 
			}
		}
	}

	vector<vector<vector<TDenseVector> > >intermediate_b(2,vector<vector<TDenseVector> >(nOriginalRegime, vector<TDenseVector>(nOriginalRegime, TDenseVector(nZ,0.0) ) ) ); 
	vector<vector<vector<TDenseMatrix> > >intermediate_F(2,vector<vector<TDenseMatrix> >(nOriginalRegime, vector<TDenseMatrix>(nOriginalRegime, TDenseMatrix(nZ,nZ,0.0) ) ) ); 
	vector<vector<vector<TDenseMatrix> > >intermediate_Phi_e(2, vector<vector<TDenseMatrix> >(nOriginalRegime, vector<TDenseMatrix>(nOriginalRegime,TDenseMatrix(nZ,nE,0.0) ) ) ); 
	vector<TDenseVector>c_hat(nOriginalRegime,TDenseVector(nZ,0.0)); 

	intermediate_b[0][0][0] = b[0]; 
	intermediate_F[0][0][0] = F[0]; 
	intermediate_Phi_e[0][0][0] = Phi_e[0]; 

	intermediate_b[1][0][0] = b[5]; 
	intermediate_F[1][0][0] = F[5]; 
	intermediate_Phi_e[1][0][0] = Phi_e[5]; 

	intermediate_b[1][1][0] = b[6]; 
	intermediate_F[1][1][0] = F[6]; 
	intermediate_Phi_e[1][1][0 ] = Phi_e[6]; 

	intermediate_b[0][1][1] = b[10]; 
	intermediate_F[0][1][1] = F[10]; 
	intermediate_Phi_e[0][1][1] = Phi_e[10]; 

	intermediate_b[1][1][1] = b[15]; 
	intermediate_F[1][1][1] = F[15]; 
	intermediate_Phi_e[1][1][1] = Phi_e[15]; 

	intermediate_b[1][0][1] = b[12]; 
	intermediate_F[1][0][1] = F[12]; 
	intermediate_Phi_e[1][0][1] = Phi_e[12]; 

	for (unsigned int i=0; i<nOriginalRegime; i++)
		c_hat[i] = LeftSolve(A[i][i]-B[i][i], C[i][i]); 

	double total = 0.0; 
	TDenseMatrix InZ = Identity(nZ); 
	for (unsigned int i=0; i<nOriginalRegime; i++)
	{
		for (unsigned int j=0; j<nOriginalRegime; j++)
		{
			for (unsigned int k=0; k<2; k++)
			{
				if (i == j || k == 1)
					total += Norm( (InZ-Pi[i][j]*GeneralizedInverse(Pi[i][j]) ) * (A[i][j]-intermediate_Phi_e[k][i][j] - Psi[i][j]) ); 
			}
		}
	}

	for (unsigned int i=0; i<nOriginalRegime; i++)
	{
		total += Norm(Va[i]*(intermediate_b[0][i][i]-c_hat[i] + intermediate_F[0][i][i]*c_hat[i]) ); 
		total += Norm(Va[i]*intermediate_F[0][i][i]*Vs[i]); 
		total += Norm(Va[i]*intermediate_Phi_e[0][i][i]); 

		unsigned int j = (i == 0) ?  1 : 0; 

		total += Norm(Va[i]*(intermediate_b[1][i][j]-c_hat[i] + intermediate_F[0][i][j]*c_hat[j]) ); 
		total += Norm(Va[i]*intermediate_F[1][i][j]); 
		total += Norm(Va[i]*intermediate_Phi_e[1][i][j]);
	}

	for (unsigned int i=0; i<nOriginalRegime; i++)
	{
		for (unsigned int j=0; j<nOriginalRegime; j++)
		{
			for (unsigned int k=0; k<2; k++)
			{
				if (i==j || k==1)
					total += Norm( (InZ-Pi[i][j]*GeneralizedInverse(Pi[i][j]) ) * ( (A[i][j]*intermediate_F[k][i][j] - B[i][j])*c_hat[j] + A[i][j]*intermediate_b[k][i][j] -C[i][j] ) ); 
				if (k==1) 
					total += Norm( (InZ-Pi[i][j]*GeneralizedInverse(Pi[i][j]) ) * (A[i][i]*intermediate_F[k][i][j]-B[i][j]) ); 
				else if (i == j)
					total += Norm( (InZ-Pi[i][j]*GeneralizedInverse(Pi[i][j]) ) * (A[i][i]*intermediate_F[k][i][k]-B[i][j]) * Vs[j]);
			}
		}
	}

	for (unsigned int i=0; i<nOriginalRegime; i++)
	{
		unsigned j = (i==0) ? 1 : 0; 
		total += Norm(p[i]*GeneralizedInverse(Pi[i][i])*A[i][i]*intermediate_b[1][i][i] - C[i][i] + (A[i][i]*intermediate_F[1][i][i] - B[i][i])*c_hat[i] + (1.0-p[i])*GeneralizedInverse(Pi[j][i])*A[j][i]*intermediate_b[1][j][i]-C[j][i]+(A[j][i]*intermediate_F[1][j][i]-B[j][i])*c_hat[i]); 
		total += Norm(p[i]*GeneralizedInverse(Pi[i][i])*(A[i][i]*intermediate_F[1][i][i]-B[i][i]) + (1.0-p[i])*GeneralizedInverse(Pi[j][i])*(A[j][i]*intermediate_F[1][j][i]-B[j][i]) ); 
		total += Norm(GeneralizedInverse(Pi[i][i])*(A[i][i]*intermediate_b[0][i][i]-C[i][i]+(A[i][i]*intermediate_F[0][i][i]-B[i][i])*c_hat[i]) ); 
		total += Norm(GeneralizedInverse(Pi[i][i])*(A[i][i]*intermediate_F[0][i][i]-B[0][0])*Vs[i]); 
	}
} 
