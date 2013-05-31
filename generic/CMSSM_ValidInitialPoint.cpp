#include <string>
#include "dw_rand.h"
#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"
#include "optimization.hpp"

using namespace std; 

double CMSSM::MINUS_INFINITY_LOCAL; 

class ObjectiveFunction_Validation
{
public:
        static CMSSM *model;
        static void *function(int *mode, int *n, double *x, double *f, double *g, int *nstate);
};

CMSSM *ObjectiveFunction_Validation::model; 

void *ObjectiveFunction_Validation::function(int *mode, int *n, double *x, double *f, double *g, int *nstate)
// A return value less than zero means that regime 1 is determinate.
{
	double error_return = -CMSSM::MINUS_INFINITY_LOCAL;  
	if (!model)
		*f = error_return; 
	else 
	{
		vector<vector<TDenseMatrix> > A, B, Psi, Pi; 
		vector<vector<TDenseVector> > C; 

		size_t nZ = model->nZ;  
		size_t nY = model->nY; 
		size_t nU = model->nU; 
		size_t nE = model->nE; 
		size_t nExpectation = model->nExpectation; 

		// Make a TDenseVector object, x_vector, to copy the content of x 
		TDenseVector x_vector(*n); 
		for (unsigned int i=0; i<*n; i++)
			x_vector.SetElement(x[i], i); 

		int gensys_err = model->rational_expectation_function->convert(A,B,Psi,Pi,C,x_vector,nZ,nY,nU,nE,nExpectation);
		if (gensys_err)
			*f = error_return;
		else 
		{	 
			TDenseMatrix V0a; 
			TDenseVector E0a; 
			Annihilator(V0a, E0a, LeftSolve(A[0][0], B[0][0])); 
			vector<double>a(2); 
			/*
 * 			E0a: contains absolute eigen values sorted in an ascending order
 * 			nZ: number of eigen values (dimension of E0a)
 * 			nExpectation: number of explosible eigen values
 * 			nZ-nExpectation: number of stable eigen values
 * 			E0a[0 ... nZ-nExpectation-1]:	stable eigen values
 * 			E0a[nZ-nExpectation ... nZ-1]:	explosible eigen values
 * 			a[0]=E0a[nZ-nExpectation-1]-1.0: difference between the smallest stable eigen value and 1, should<0 if everthing else works 
 * 			a[1]=1.0-E0a[nZ-nExpectation]: difference between the smallest explosible eigen value and 1, should<0 if everthing else works 
 *			we would like both a[0] and a[1] to be negative, not not too negative (~0.2 should be perfect)
 *			However, if nExpectation is not 2 due to a particular x, then either a[0] or a[1] will be positive. We would like to avoid 
 *			this x by letting nposl keep driving *f down. 
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
 * 				if (*f < -0.2)
 *					*f = -0.2; 
				*/
				*f = 2.5*((*f+0.2)*(*f+0.2)-0.04); 
			}
		}
	}	
}

int CMSSM::ValidInitialPoint(TDenseVector &fval, TDenseVector &x_optimal, const TDenseVector &x0, size_t max_count, const TDenseVector &lower_bound, const TDenseVector &upper_bound)
// Tries to find an x that is an valid initial point
//
// Return value 
// 	MODEL_NOT_PROPERLY_SET:     state_equation_function, measurement_equation_function or transition_prob_function not properly set
// 	SUCCESS:	success (no error)
// 	ERRO_OCCURRED:	error occurred
//
// Also returned:
// 	x:	valid initial point if successful
// 	fval:	vector of objective function values. A value <0 means that a determinate solution in regime 1 has been found
{
	if (CheckModelFunctions() != SUCCESS) 
		return MODEL_NOT_PROPERLY_SET; 

	int error;
	const double INFINITE_BOUND = -MINUS_INFINITY_LOCAL; 
	const double TOLERANCE = 0.0; 
	const string COLD_START = string("Cold Start"); 
	const string NO_PRINT_OUT = string("Major print level = 0");
	const string DERIVATIVE_LEVEL = string("Derivative level = 0"); 

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
	// lower and upper bounds
	for (unsigned int i=0; i<n; i++)
	{
		bl[i] = lower_bound.dim > i ? lower_bound[i] : -INFINITE_BOUND; 
		bu[i] = upper_bound.dim > i ? upper_bound[i] : INFINITE_BOUND;  
	}
	// since the last 2 parameters correspond to probabilites
	for (unsigned int i=n-1; i>=n-2; i--)
        {
                if (bl[i] < 1.0e-3)
			bl[i] = 1.0e-3; 
		if (bu[i] > 1.0)
                        bu[i] = 1.0; 
        }


 	ObjectiveFunction_Validation::model = this; 
	// Below is a revision based on Dan's code
	// We try a number, max_count, of times and finds the best solution
	double best_value = -MINUS_INFINITY_LOCAL; 
	double *best_x = new double[n];
	while (count < max_count)
	{
		npoptn_((char*)DERIVATIVE_LEVEL.c_str(), DERIVATIVE_LEVEL.length());
	        npoptn_((char*)COLD_START.c_str(), COLD_START.length());
        	npoptn_((char*)NO_PRINT_OUT.c_str(), NO_PRINT_OUT.length());
		npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, ObjectiveFunction_Validation::function, &inform, &iter, istate, c, cJac, clamda, &f, g, R, x_raw, iw, &leniw, w, &lenw);
		fval.SetElement(f, count); 
		if (f < best_value)
		{
			best_value = f; 
			memcpy(best_x, x_raw, sizeof(double)*n); 
		}
		for (unsigned int i=0; i<x0.dim; i++)
		{
                	if (bl[i] > -INFINITE_BOUND)
                	{
                		if (bu[i] < INFINITE_BOUND )
                			x_raw[i] = (bu[i]-bl[i])*dw_uniform_rnd() + bl[i]; 
                		else 
                		{
                			double grn_1 = dw_gaussian_rnd(), grn_2 = dw_gaussian_rnd(); 
                			x_raw[i] = (grn_1*grn_1 + grn_2*grn_2) + bl[i]; 
                		}
                	}
                	else 
                	{
                		if (bu[i] < INFINITE_BOUND) 
                		{
                			double grn_1 = dw_gaussian_rnd(), grn_2 = dw_gaussian_rnd(); 
                			x_raw[i] = bu[i] - (grn_1*grn_1 + grn_2*grn_2); 
                		}
                		else 
                			x_raw[i] = dw_gaussian_rnd(); 
                	}
                }
                count ++; 
	}
	if (best_value < TOLERANCE)
	{
		error = SUCCESS; 
		memcpy(x_raw, best_x, sizeof(double)*n); 
	}
	else
		error = ERROR_OCCURRED; 
	delete [] best_x; 
	
 	/* Below is almost the exact copy of Dan's code
	while (error && count < max_count)
	{
		npoptn_((char*)COLD_START.c_str(), COLD_START.length()); 
		npoptn_((char*)NO_PRINT_OUT.c_str(), NO_PRINT_OUT.length()); 
		npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, ObjectiveFunction_Validation::function, &inform, &iter, istate, c, cJac, clamda, &f, g, R, x_raw, iw, &leniw, w, &lenw); 
		fval.SetElement(f, count); 
		if (fval[count] < TOLERANCE)
			error = SUCCESS;  
		else 
		{
			for (unsigned int i=0; i<x0.dim; i++)
			{
				if (bl[i] > -INFINITE_BOUND)
				{
					if (bu[i] < INFINITE_BOUND )
						x_raw[i] = (bu[i]-bl[i])*dw_uniform_rnd() + bl[i]; 
					else 
					{
						double grn_1 = dw_gaussian_rnd(), grn_2 = dw_gaussian_rnd(); 
						x_raw[i] = (grn_1*grn_1 + grn_2*grn_2) + bl[i]; 
					}
				}
				else 
				{
					if (bu[i] < INFINITE_BOUND) 
					{
						double grn_1 = dw_gaussian_rnd(), grn_2 = dw_gaussian_rnd(); 
						x_raw[i] = bu[i] - (grn_1*grn_1 + grn_2*grn_2); 
					}
					else 
						x_raw[i] = dw_gaussian_rnd(); 
				}
			}
			count ++; 
		}
	} */
	if (error == SUCCESS)
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
