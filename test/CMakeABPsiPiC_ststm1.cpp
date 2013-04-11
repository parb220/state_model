#include <vector>
#include "CMakeABPsiPiC_ststm1.hpp"

using namespace std; 

bool MakeABPsiPiC_ststm1::convert(vector<vector<TDenseMatrix> > &A, vector<vector<TDenseMatrix> > &B, vector<vector<TDenseMatrix> > &Psi, vector<vector<TDenseMatrix> >&Pi, vector<vector<TDenseVector> > &C, const TDenseVector &x, const TDenseVector &free_parameter)
// Microfounded NK model with Rotemburg's adjustment costs
//
// Return gensys form for state equation: 
//
// 	A(s(t))*x(t) = C(s(t)) + B(s(t))*x(t-1) + Psi(s(t))*epsilon(t) + Pi(s(t))*eta(t)
//
// err is one if some A{i} is not of full rank and zero otherwise.
{
	bool error_code = false; 
	unsigned i=0, j=x.dim; 
	size_t nRegime = 2; 
	
	double gsigma = x[i];	i++;
	TDenseVector gpsi2(nRegime,0.0); gpsi2.SetElement(x[i],0); gpsi2.SetElement(x[j],1);  i++; j++; 
	TDenseVector b(nRegime,0.0); b.SetElement(x[i],0); b.SetElement(x[j],1); i++; j++; 
	TDenseVector ggamma(nRegime,0.0); ggamma.SetElement(x[i],0); ggamma.SetElement(x[j],1); i++; j++; 
	double geta = x[i]; i++; 
	TDenseVector gphipi(nRegime,0.0); gphipi.SetElement(x[i],0); i++; 
	TDenseVector gphiy(nRegime,0.0); gphiy.SetElement(x[i],0); i++; 
	TDenseVector grhomu(nRegime,0.0); grhomu.SetElement(x[i],0); grhomu.SetElement(x[j],1); i++; j++; 
	TDenseVector grhorn(nRegime,0.0); grhorn.SetElement(x[i],0); grhorn.SetElement(x[j],1); i++; j++;	// Natural rate of interest
	TDenseVector grhoi(nRegime,0.0); grhoi.SetElement(x[i],0); grhoi.SetElement(x[j],1); i++; j++; 	// Monetary policy
	TDenseVector gsigmamu(nRegime,0.0); gsigmamu.SetElement(x[i],0); gsigmamu.SetElement(x[j],1); i++; j++; 
	TDenseVector gsigmarn(nRegime,0.0); gsigmarn.SetElement(x[i],0); gsigmarn.SetElement(x[j],1); i++; j++; // Natural rate of interest
	TDenseVector gsigmai(nRegime,0.0); gsigmai.SetElement(x[i],0); gsigmai.SetElement(x[j],1); i++; j++; 	// Monetary policy
	double giota = x[i]; 	// ZLB adjustment for level of fall in natural rate of ineterest	

	double glambdaz = free_parameter[GLAMBDAZ_STATE];	// glambdaz, Gross rate: (1+2%) annually per capita or (1+0.5%) quarterly per capita 
	double rn = free_parameter[RN_STATE];	// rn, Net rate: rn=4% annually or 1% quarterly (natural rate of interest or steady state real interest rate)
       	double gpistar = free_parameter[GPISTAR_STATE];	// gpistar, Net rate: log(1.02) -- 2% annually or 0.5% quarterly
        double sc = free_parameter[SC_STATE];	// sc, Steady state share of private consumption in C+G
        double Rlow = free_parameter[RLOW_STATE];	// Rlow, Net rate 10 basis points (0.1%) interest rate annually at zero bound or 0.025% quarte	

	double Rstar = gpistar + rn;	// corresponding to i in Zha's notes
	double gbeta = glambdaz / (1.0+rn);	// 0.995 because gbeta=glambdaz/(1+rn) where rn=4% annually (natural rate of interest or steady state real interest rate)
	double istar = -Rstar + Rlow;	// correspond to i^* in Zha's notes
	TDenseVector ci(nRegime,0.0); ci.SetElement((1.0-grhoi[1])*istar,1);
	TDenseVector crn(nRegime,0.0); crn.SetElement((1.0-grhorn[1])*(giota-Rstar),1);
	TDenseVector gchir(nRegime,0.0); gchir.SetElement(1.0,0); 

	size_t n_ = 7; 
	A = vector<vector<TDenseMatrix> >(nRegime, vector<TDenseMatrix>(nRegime, TDenseMatrix(n_,n_,0.0) ) ); 
	B = vector<vector<TDenseMatrix> >(nRegime, vector<TDenseMatrix>(nRegime, TDenseMatrix(n_,n_,0.0) ) ); 
	C = vector<vector<TDenseVector> >(nRegime, vector<TDenseVector>(nRegime, TDenseVector(n_,0.0) ) ); 
	Psi = vector<vector<TDenseMatrix> >(nRegime, vector<TDenseMatrix>(nRegime, TDenseMatrix(n_,3,0.0) ) ); 
	Pi = vector<vector<TDenseMatrix> >(nRegime, vector<TDenseMatrix>(nRegime, TDenseMatrix(n_,2,0.0) ) ); 

	for (unsigned int i=0; i<nRegime; i++)
	{
		for (unsigned int j=0; j<nRegime; j++)
		{
			A[i][j].SetElement(1.0+gbeta*ggamma[i],0,0); 
			A[i][j].SetElement(-gpsi2[i]*(geta+gsigma/sc),0,1);
			A[i][j].SetElement(-gpsi2(i),0,3);  
			A[i][j].SetElement(-gbeta,0,5);
			A[i][j].SetElement(gsigma+(gsigma-1.0)*b[i],1,1);
			A[i][j].SetElement(sc,1,2);
			A[i][j].SetElement(-sc,1,4);
			A[i][j].SetElement(-sc,1,5);
			A[i][j].SetElement(-gsigma,1,6);
			A[i][j].SetElement(-(1.0-grhoi[i])*gphipi[i],2,0);
			A[i][j].SetElement(-(1.0-grhoi[i])*gphiy[i],2,1);
			A[i][j].SetElement(1.0,2,2);
			A[i][j].SetElement(-(1.0-grhoi[i])*gchir[i],2,4);
			A[i][j].SetElement(1.0,3,4);
			A[i][j].SetElement(1.0,4,3);
			A[i][j].SetElement(1.0,5,0);
			A[i][j].SetElement(1.0,6,1);

			B[i][j].SetElement(ggamma[j],0,0);
			B[i][j].SetElement(-gpsi2[j]*(gsigma-1.0)*b[j]/sc,0,1);	
			B[i][j].SetElement((gsigma-1.0)*b[j],1,1);
			B[i][j].SetElement(grhoi[i],2,2);
			B[i][j].SetElement(grhorn[i],3,4);
			B[i][j].SetElement(grhomu[i],4,3);
			B[i][j].SetElement(1.0,5,5);
			B[i][j].SetElement(1.0,6,6);

			C[i][j].SetElement(ci[i],2); 
			C[i][j].SetElement(crn[i],3);

			Psi[i][j].SetElement(gsigmai[i],2,2);
			Psi[i][j].SetElement(gsigmarn[i],3,0);
			Psi[i][j].SetElement(gsigmamu[i],4,1);
			
			Pi[i][j].SetElement(1.0,5,0);
			Pi[i][j].SetElement(1.0,6,1);

			if (Rank(A[i][j]) < n_)
				error_code = true;  
		}
	}
	return error_code;
}
