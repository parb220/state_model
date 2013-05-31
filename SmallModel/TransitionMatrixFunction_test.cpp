#include <cstdlib>
#include <cmath>
#include "CMSSM_Error_Code.hpp"
#include "TransitionMatrixFunction_test.hpp" 

using namespace std; 

int TransitionProbMatrixFunction_test :: convert(TDenseMatrix &Q, unsigned int t, const vector<TDenseVector> &y, size_t nS, size_t nTL, const TDenseVector &x) 
{
	// Transition matrix from path of regimes
	TDenseMatrix Qs(nS, nS, 0.0); 
	Qs.SetElement(1.0-x[x.dim-1],0,3); 
	Qs.SetElement(x[x.dim-2],1,1); 
	Qs.SetElement(1.0-x[x.dim-2],2,1); 
	Qs.SetElement(1.0,2,2); 
	Qs.SetElement(x[x.dim-1],3,3); 

	if (y[t][2] >= fixed_parameter[2])	// fixed_parameter[2]: cutoff
		Qs.SetElement(1.0,0,0);
	else 
		Qs.SetElement(1.0,1,0); 
	
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
						Q.SetElement(Qs(i,j), i+nS*(j+nS*k),j+nS*(k+n*m));	
	}
	else 
		Q = Qs; 
	return SUCCESS; 
}
