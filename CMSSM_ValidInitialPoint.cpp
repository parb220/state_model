#include "dw_rand.h"
#include "CMSSM.hpp"

extern "C"
{
	void npsol_(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A, double *bl, double *bu, void *funcon(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate), void *funobj(int *mode, int *n, double *x, double *f, double *g, int *nstate), int *inform, int *iter, int *istate, double *c, double *cJac, double *clamda, double *f, double *g, double *R, double *x, int *iw, int *leniw, double *w, int *lenw);
}

using namespace std; 

MakeABPsiPiC* ObjectiveFunction_Validation::ABPsiPiC_function; 
TDenseVector ObjectiveFunction_Validation::state_equation_parameter; 
TDenseVector ObjectiveFunction_Validation::measurement_equation_parameter;
TDenseVector ObjectiveFunction_Validation::transition_prob_parameter; 

void *ObjectiveFunction_Validation::function(int *mode, int *n, double *x, double *f, double *g, int *nstate)
// A return value less than zero means that regime 1 is determinate.
{
	double error_return = 1.0e10;
	vector<vector<TDenseMatrix> > A, B, Psi, Pi; 
	vector<vector<TDenseVector> > C;  
	TDenseVector x_vector(x, *n); 

	bool gensys_err = ABPsiPiC_function->convert(A,B,Psi,Pi,C,x_vector,state_equation_parameter);
	if (gensys_err)
	{
		*f = error_return;
		return NULL; 
	} 
 
	size_t nZ = A[0][0].rows; 
	size_t nExpectation = Pi[0][0].cols; 
	TDenseMatrix V0a; 
	vector<double> E0a; 
	Annihilator(V0a, E0a, LeftSolve(A[0][0], B[0][0]) ); 
	if (V0a.rows != nExpectation) 
	{
		*f = error_return; 
		return NULL; 
	}
	if (Rank(V0a*LeftSolve(A[0][0],Pi[0][0]) ) != nExpectation)
	{
		*f = error_return; 
		return NULL; 
	}

	vector<double>a(2); 
	a[0] = E0a[nZ-nExpectation-1]-1.0; 	// a[0]=E0a[nZ-nExpectation]-1.0; 
	a[1] = 1.0-E0a[nZ-nExpectation]; 	// a[1]=1.0-E1a[nZ-nExpectation+1]; 
	vector<double>b(2); 
	b[0] = a[0] > 0.0 ? a[0] : 0.0; 
	b[1] = a[1] > 0.0 ? a[1] : 0.0; 

	*f=b[0]+b[1]; 
	if (*f <= 0)
	{
		*f = a[0] > a[1] ? a[0] : a[1]; 
		if (*f < -0.2)
			*f = -0.2; 
	}
	return NULL; 
}

bool CMSSM::ValidInitialPoint(TDenseVector &fval, const TDenseVector &x0, size_t max_count, TDenseVector &lb, TDenseVector &ub)
// Tries to find an x that is an valid initial point
//
// Return value 
// 	false:	success (no error)
// 	true:	error occurred
//
// Also returned:
// 	x:	valid initial point if successful
// 	fval:	vector of objective function values. A value <0 means that a determinate solution in regime 1 has been found
{
	bool error = true;
	const double INFINITE_BOUND = 10E10; 

	if (lb.dim != x0.dim)
		lb = TDenseVector(-INFINITE_BOUND, x0.dim); 
	if (ub.dim != x0.dim)
		ub = TDenseVector(INFINITE_BOUND, x0.dim); 

	for (unsigned int i=0; i<2; i++)
	{
		if (lb[lb.dim-i-1] < 1.0e-3)
			lb.SetElement(1.0e-3, lb.dim-i-1); 
		if (ub[ub.dim-i-1] > 1.0) 
			ub.SetElement(1.0, ub.dim-i-1); 
	}

	fval = TDenseVector(max_count, 0.0); 
	unsigned int count = 0; 
	
	// ====== variables for npsol ======
	int n = x.dim; 
	int nclin = 0; 	// number of linear constraints
	int ncnln = 0; 	// number of nonlinear constraints
	int nctotal = n + nclin + nctotal; 
	int ldA = 1;	// row dimension of A, because nclin=0
	int ldJ = 1; 	// row dimension of cJac, because ncnln = 0; 
	int ldR = n;	// row dimension of R 	
	double *A = new double[ldA*n];	// linear constraint matrix
	double *bl = new double[nctotal]; 	// lower bound of all constraints
	double *bu = new double[nctotal]; 	// upper bound of all constraints
	int *istate = new int[nctotal]; 	// indicate which constaints are used
	double *cJac = new double[ldJ*n]; 	// nonlinear constraint matrix
	double *clamda = new double[nctotal];	// Lagrangian coefficients
	double *R = new double[ldR*n]; 		// Hessian of Lagrangian
	double *x_raw = new double[n]; 		// solution
	int leniw = 3*n*nclin+2*ncnln; 
	int *iw = new int[leniw];  
	int lenw; 
	if (nclin==0 && ncnln ==0)
		lenw = 20*n; 
	else if (ncnln == 0)
		lenw = 2*n*n+20*n+11*nclin; 
	else 
		lenw = 2*n*n+n*nclin+2*n*ncnln+20*n+11*nclin+21*ncnln; 
	double *w = new double[lenw]; 
	int inform, iter; 
	double *c = new double [1];	// contains the values of the nonlinear constraint functions, because ncnln=0; 
	double f; 
	double *g = new double[n];	// contains objective gradient 
	// =================================
	
	memcpy(x_raw, x0.vector, sizeof(double)*x0.dim); 
	for (unsigned int i=0; i<n; i++)
	{
		bl[i] = lb[i];			// Setting lower bound  
		bu[i] = ub[i]; 			// Setting upper bound 
	}
	while (error && count < max_count)
	{
		npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, (ObjectiveFunction_Validation::function), &inform, &iter, istate, c, cJac, clamda, &f, g, R, x_raw, iw, &leniw, w, &lenw); 

		fval.SetElement(f, count); 
		if ( inform == 0 && fval[count] > 0)
		{
			for (unsigned int i=0; i<x.dim; i++)
			{
				if (lb[i] > -INFINITE_BOUND)
				{
					if (ub[i] < INFINITE_BOUND )
						x_raw[i] = (ub[i]-lb[i])*dw_uniform_rnd() + lb[i]; 
					else 
					{
						double grn_1 = dw_gaussian_rnd(), grn_2 = dw_gaussian_rnd(); 
						x_raw[i] = (grn_1*grn_1 + grn_2*grn_2) + lb[i]; 
					}
				}
				else 
				{
					if (ub[i] < INFINITE_BOUND) 
					{
						double grn_1 = dw_gaussian_rnd(), grn_2 = dw_gaussian_rnd(); 
						x_raw[i] = ub[i] - (grn_1*grn_1 + grn_2*grn_2); 
					}
					else 
						x_raw[i] = dw_gaussian_rnd(); 
				}
			}
			count ++; 
		}
		else 
			error = false; 
	}

	for (unsigned int i=0; i<x.dim; i++)
		x.SetElement(x_raw[i], i); 

	// ========= release npsol variables ==================
	delete [] A; 
	delete [] bl; 
	delete [] bu; 
	delete [] istate; 
	delete [] cJac; 
	delete [] clamda; 
	delete [] R; 
	delete [] x_raw; 
	delete [] iw; 
	delete [] w; 
	delete [] c; 
	delete [] g; 
	// ====================================================
}
