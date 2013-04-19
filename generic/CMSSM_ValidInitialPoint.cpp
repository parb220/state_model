#include <string>
#include "dw_rand.h"
#include "CMSSM.hpp"

extern "C"
{
	void npsol_(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A, double *bl, double *bu, void *funcon(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate), void *funobj(int *mode, int *n, double *x, double *f, double *g, int *nstate), int *inform, int *iter, int *istate, double *c, double *cJac, double *clamda, double *f, double *g, double *R, double *x, int *iw, int *leniw, double *w, int *lenw);
	void npoptn_(char *, int); 
}

using namespace std; 

CMSSM *ObjectiveFunction_Validation::model;

void *ObjectiveFunction_Validation::function(int *mode, int *n, double *x, double *f, double *g, int *nstate)
// A return value less than zero means that regime 1 is determinate.
{
	double error_return = 1.0e10;
	if (!model)
		*f = error_return; 
	else 
	{
		vector<vector<TDenseMatrix> > A, B, Psi, Pi; 
		vector<vector<TDenseVector> > C; 

		// Make a TDenseVector object, x_vector, to copy the content of x 
		TDenseVector x_vector(*n); 
		for (unsigned int i=0; i<*n; i++)
			x_vector.SetElement(x[i], i); 

		bool gensys_err = model->state_equation_function->convert(A,B,Psi,Pi,C,model->state_equation_parameter, x_vector);
		if (gensys_err)
			*f = error_return;
		else 
		{	 
			size_t nZ = A[0][0].rows; 
			size_t nExpectation = Pi[0][0].cols; 
			TDenseMatrix V0a; 
			vector<double> E0a; 
			Annihilator(V0a, E0a, LeftSolve(A[0][0], B[0][0]) ); 
			/*if (V0a.rows != nExpectation) 
				*f = error_return; 
			else if (Rank(V0a*LeftSolve(A[0][0],Pi[0][0]) ) != nExpectation)
				*f = error_return; 
			else 
			{*/
				vector<double>a(2); 
				/*
 * 				E0a: contains absolute eigen values sorted in an ascending order
 * 				nZ: number of eigen values (dimension of E0a)
 * 				nExpectation: number of explosible eigen values
 * 				nZ-nExpectation: number of stable eigen values
 * 				E0a[0 ... nZ-nExpectation-1]:	stable eigen values
 * 				E0a[nZ-nExpectation ... nZ-1]:	explosible eigen values
 * 				a[0]=E0a[nZ-nExpectation-1]-1.0: difference between the smallest stable eigen value and 1, should<0 if everthing else works 
 * 				a[1]=1.0-E0a[nZ-nExpectation]: difference between the smallest explosible eigen value and 1, should<0 if everthing else works 
 *				we would like both a[0] and a[1] to be negative, not not too negative (~0.2 should be perfect)
 *				However, if nExpectation is not 2 due to a particular x, then either a[0] or a[1] will be positive. We would like to avoid 
 *				this x by letting nposl keep driving *f down. 
 				*/
				a[0] = E0a[nZ-nExpectation-1]-1.0; 	// a[0]=E0a[nZ-nExpectation]-1.0; 
				a[1] = 1.0-E0a[nZ-nExpectation]; 	// a[1]=1.0-E1a[nZ-nExpectation+1]; 
				vector<double>b(2); 
				b[0] = a[0] > 0.0 ? a[0] : 0.0; 
				b[1] = a[1] > 0.0 ? a[1] : 0.0; 

				*f=b[0]+b[1]; 
				if (*f <= 0)
				{
					*f = a[0] > a[1] ? a[0] : a[1]; 
					/* The following seggments make the objective function not differentiable
 * 					if (*f < -0.2)
 *						*f = -0.2; 
					*/
					*f = 2.5*((*f+0.2)*(*f+0.2)-0.04); 
				}
			/*}*/
		}
	}	
}

int CMSSM::ValidInitialPoint(TDenseVector &fval, TDenseVector &x_optimal, const TDenseVector &x0, size_t max_count, TDenseVector &lb, TDenseVector &ub)
// Tries to find an x that is an valid initial point
//
// Return value 
// 	-1:     state_equation_function, measurement_equation_function or transition_prob_function not properly set
// 	0:	success (no error)
// 	1:	error occurred
//
// Also returned:
// 	x:	valid initial point if successful
// 	fval:	vector of objective function values. A value <0 means that a determinate solution in regime 1 has been found
{
	if (CheckStateMeasurementTransitionEquations()) 
		return -1; 

	int error = 1;
	const double INFINITE_BOUND = 1.0E20; 
	const double TOLERANCE = 0.0; 
	const string COLD_START = string("Cold Start"); 

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
	int n = x0.dim; 
	int nclin = 0; 	// number of linear constraints
	int ncnln = 0; 	// number of nonlinear constraints
	int nctotal = n + nclin + nctotal; 
	int ldA = nclin > 1 ? nclin : 1;	// row dimension of A, because nclin=0
	int ldJ = ncnln > 1 ? ncnln : 1; 	// row dimension of cJac, because ncnln = 0; 
	int ldR = n;	// row dimension of R 	
	double *A = new double[ldA*n];	// linear constraint matrix
	double *bl = new double[nctotal]; 	// lower bound of all constraints
	double *bu = new double[nctotal]; 	// upper bound of all constraints
	int *istate = new int[nctotal]; 	// indicate which constaints are used
	double *cJac = new double[ldJ*n]; 	// nonlinear constraint matrix
	double *clamda = new double[nctotal];	// Lagrangian coefficients
	double *R = new double[ldR*n]; 		// Hessian of Lagrangian
	double *x_raw = new double[n]; 		// solution
	int leniw = 3*n+nclin+2*ncnln; 
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
	
	memcpy(x_raw, x0.vector, sizeof(double)*n); 
	for (unsigned int i=0; i<n; i++)
	{
		bl[i] = lb[i];			// Setting lower bound  
		bu[i] = ub[i]; 			// Setting upper bound 
	}

 	ObjectiveFunction_Validation::model = this; 
	// Below is a revision based on Dan's code
	// We try a number, max_count, of times and finds the best solution
	double best_value = 1.0; 
	double *best_x = new double[n];
	while (count < max_count)
	{
		npoptn_((char*)COLD_START.c_str(), COLD_START.length());
		npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, ObjectiveFunction_Validation::function, &inform, &iter, istate, c, cJac, clamda, &f, g, R, x_raw, iw, &leniw, w, &lenw);
		fval.SetElement(f, count); 
		if (inform == 0 || inform == 1)
		{
			if (f < best_value)
			{
				best_value = f; 
				memcpy(best_x, x_raw, sizeof(double)*n); 
			}
		}
		for (unsigned int i=0; i<x0.dim; i++)
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
	if (best_value < TOLERANCE)
	{
		cout << "best_value "<< best_value << endl; 
		error = 0; 
		memcpy(x_raw, best_x, sizeof(double)*n); 
	}
	else
		error = 1; 
	delete [] best_x; 
	
 	/* Below is almost the exact copy of Dan's code
	while (error && count < max_count)
	{
		npoptn_((char*)COLD_START.c_str(), COLD_START.length()); 
		npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, ObjectiveFunction_Validation::function, &inform, &iter, istate, c, cJac, clamda, &f, g, R, x_raw, iw, &leniw, w, &lenw); 

		fval.SetElement(f, count); 
		if (fval[count] < TOLERANCE)
		// if (inform == 0 && fval[count] < TOLERANCE )
			error = 0; 
		else 
		{
			for (unsigned int i=0; i<x0.dim; i++)
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
	} */
	if (error == 0)
	{	
		x_optimal.Resize(x0.dim); 
		for (unsigned int i=0; i<x0.dim; i++)
			x_optimal.SetElement(x_raw[i], i); 
	}

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
	return error; 
}
