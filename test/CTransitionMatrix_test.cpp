#include <cstdlib>
#include <cmath>
#include "CTransitionMatrix_test.hpp"

using namespace std; 

void CTransitionProbMatrixFunction_Test :: convert(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, const TDenseVector &x, size_t nS, size_t nTL, const TDenseVector &free_parameter) 
{
	if (nS != 4)
	{
		cerr << "nS in CTransitionProbMatrixFunction_Test::convert() not equal to 4\n"; 
		abort(); 
	}

	if (x.dim < 3)
	{
		cerr << "dimension of x in CTransitionProbMatrixFunction_Test::convert() less than 3, probably mixing old and new code.\n"; 
		abort(); 
	}

	// Transition matrix from path of regimes
	TDenseMatrix Qs(nS, nS, 0.0); 
	Qs.SetElement(0,3,1.0-x[x.dim-1]); 
	Qs.SetElement(1,1,x[x.dim-2]); 
	Qs.SetElement(2,1,1.0-x[x.dim-2]); 
	Qs.SetElement(2,2,1.0); 
	Qs.SetElement(3,3,x[x.dim-1]); 

	if (y[t][2] >= free_parameter[CUTOFF_TRANSITION])	// free_parameter[0]: cutoff
		Qs.SetElement(0,0,1.0);
	else 
		Qs.SetElement(1,0,1.0); 
	
	//
	// xi(t) = [s(t), s(t-1), ..., s(t-nTL)] encoded as 
	//
	// 	s(t)+nS*s(t-1)+...+nS^nTL*s(t-nTL)
	//
	// zeta(t) = [s(t), s(t-1), ...., s(t-nTL+1)]  encoded as
	//	
	//	s(t)+nS*s(t-1)+...+nS^(nTL-1)*s(t-nTL+1)
	//
	if (nTL >= 1)
	{
		size_t n = (size_t)pow(nS, nTL-1);  
		size_t nXi = (size_t)pow(nS, nTL+1); 
		Q.Zeros(nXi,nXi);
		for (unsigned int i=0; i<nS; i++)		// i = s(t+1)
			for (unsigned int j=0; j<nS; j++)	// j = s(t)
				for (unsigned int k=0; k<n; k++)	// k = s(t-1)+nS*s(t-2)+...nS^(nTL-1)*s(t-nTL+1)
					for (unsigned m=0; m<nS; m++)	// m = s(t-nTL)
						Q.SetElement(i+nS*(j+nS*k),j+nS*(k+n*m), Qs(i,j));	
	}
	else 
		Q = Qs; 
}
