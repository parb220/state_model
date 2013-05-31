#include <string>
#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"
#include "optimization.hpp"

using namespace std; 
class MinusLogPosterior_NPSOL
{
public:
        static CMSSM *model;
        static vector<TDenseVector> y;
        static vector<TDenseVector> z_0;
        static vector<TDenseMatrix> P_0;
        static TDenseVector initial_prob;

        static void *function(int *mode, int *n, double *x, double *f, double *g, int *nstate);
};

class MinusLogPosterior_CSMINWEL
{
public: 
	static CMSSM *model; 
	static vector<TDenseVector> y; 
	static vector<TDenseVector> z_0; 
	static vector<TDenseMatrix> P_0; 
	static TDenseVector initial_prob; 

	static double function(double *x, int n, double **args, int *dims); 
}; 

CMSSM * MinusLogPosterior_NPSOL::model; 
vector<TDenseVector> MinusLogPosterior_NPSOL::y; 
vector<TDenseVector> MinusLogPosterior_NPSOL::z_0; 
vector<TDenseMatrix> MinusLogPosterior_NPSOL::P_0; 
TDenseVector MinusLogPosterior_NPSOL::initial_prob; 

void * MinusLogPosterior_NPSOL::function(int *mode, int *n, double *x_array, double *f, double *g, int *nstate)
{
	// Make TDenseVector out of x_array
        TDenseVector x(*n); 
	for (unsigned int i=0; i<*n; i++)
		x.SetElement(x_array[i],i); 

        double minus_log_posterior= -CMSSM::MINUS_INFINITY_LOCAL;
	double log_posterior;
	if (model->LogPosterior(log_posterior, x, y, z_0, P_0, initial_prob) == SUCCESS)
		minus_log_posterior = -log_posterior; 

	*f = minus_log_posterior; 
}

CMSSM * MinusLogPosterior_CSMINWEL::model; 
vector<TDenseVector> MinusLogPosterior_CSMINWEL::y; 
vector<TDenseVector> MinusLogPosterior_CSMINWEL::z_0; 
vector<TDenseMatrix> MinusLogPosterior_CSMINWEL::P_0; 
TDenseVector MinusLogPosterior_CSMINWEL::initial_prob; 

double MinusLogPosterior_CSMINWEL::function(double *x_array, int n, double **args, int *dim)
{
	TDenseVector x(n); 
	for (unsigned int i=0; i<n; i++)
		x.SetElement(x_array[i],i); 

	double minus_log_posterior = -CMSSM::MINUS_INFINITY_LOCAL, log_posterior; 
	if (model->LogPosterior(log_posterior, x, y, z_0, P_0, initial_prob) == SUCCESS)
		return -log_posterior; 
	else 
		return minus_log_posterior; 
}

int CMSSM::Maximize_LogPosterior_NPSOL(double &log_posterior_optimal, TDenseVector &x_optimal, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob, const TDenseVector &x0) 
// Returns
// 	MODEL_NOT_PROPERLY_SET:	if state_equation_function, measurement_equation_function or transition_prob_function not properly set
// 	>=0:	inform code returned by npsol_
{
	if (prior_distr_function == NULL || CheckModelFunctions() != SUCCESS)
	{
		log_posterior_optimal = MINUS_INFINITY_LOCAL; 
		return MODEL_NOT_PROPERLY_SET; 
	}
 
	int error; 

	// Setting MinusLogPosterior static variables
	MinusLogPosterior_NPSOL::model = this; 
	MinusLogPosterior_NPSOL::y = y; 
	MinusLogPosterior_NPSOL::z_0 = z_0; 
	MinusLogPosterior_NPSOL::P_0 = P_0; 
	MinusLogPosterior_NPSOL::initial_prob = initial_prob; 

	const double INFINITE_BOUND = -MINUS_INFINITY_LOCAL ;
	const string COLD_START = string("Cold Start");
	const string NO_PRINT_OUT = string("Major print level = 0"); 
	const string DERIVATIVE_LEVEL = string("Derivative level = 0"); 
	// npsol unconstrained 
	int n = x0.dim; 
	int nclin = 0; 
	int ncnln = 0;
	int nctotl = n + nclin + ncnln;  
	int ldA = nclin > 1 ? nclin : 1;	// ldA >= 1 and ldA >= nclin, row dimension of A
	int ldJ = ncnln > 1 ? ncnln : 1; 	// ldJ >= 1 and ldJ >= ncnln, row dimenision of cJac
	int ldR = n; 				// row dimension of R
	double *A = new double[ldA*n];		// linear constraint matrix
	double *bl = new double[nctotl];	// lower bound of all constraints
	double *bu = new double[nctotl];	// upper bound of all constraints
	int *istate = new int[nctotl]; 	// indicate whether bl or bu or both should be used
	double *cJac = new double[ldJ*n]; 	// Jacobian matrix for the nonlinear constraint
	double *clambda = new double[nctotl];	// Lagrangin coefficients (not relevant here because our problem is unconstrained)
	double *R = new double[ldR*n]; 		// Hessian of Lagrangian
	double *x_raw = new double[n]; 		
	int leniw = 3*n+nclin+2*ncnln; 
	int *iw = new int[leniw];		// integer working space
	int lenw; 
	if (nclin == 0 && ncnln == 0)
		lenw = 20*n; 
	else if (ncnln == 0)
		lenw = 2*n*n+20*n+11*nclin; 
	else 
		lenw = 2*n*n+n*nclin+2*n*ncnln+20*n+11*nclin+21*ncnln; 
	double *w = new double[lenw]; 			// double working space
	int inform, iter; 
	double *c = new double[1]; 		// value of nonlinear constraint function, not used because ncnln = 0; 
	double f; 				// value of the objective function
	double *g = new double[n]; 			// gradident 
	
	//=====================================================
	memcpy(x_raw, x0.vector, n*sizeof(double) ); 
	for (unsigned int i=0; i<n; i++)
	{
		bl[i] = -INFINITE_BOUND; 
		bu[i] = INFINITE_BOUND; 
	}

	npoptn_((char*)DERIVATIVE_LEVEL.c_str(), DERIVATIVE_LEVEL.length()); 
	npoptn_((char*)COLD_START.c_str(), COLD_START.length());
	npoptn_((char*)NO_PRINT_OUT.c_str(), NO_PRINT_OUT.length());
	npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, MinusLogPosterior_NPSOL::function, &inform, &iter, istate, c, cJac, clambda, &f, g, R, x_raw, iw, &leniw, w, &lenw); 

	x_optimal.Resize(n); 
	for (unsigned int i=0; i<n; i++)
		x_optimal.SetElement(x_raw[i], i); 
	log_posterior_optimal = -f; 
	error = inform; 
	// inform returned by npsol_ has the following meanings:
	// 	<0:	Either funcon or funobj has set mode to this negative value
	// 	0:	The iterates have converged to x that satisfies all optimal conditions.
	// 	1:	x satisfies all optimal conditions, but the iterates have not converged
	// 	2:	No feasible solution because the linear constraints and bounds cannot be met
	// 	3:	No feasible solution because the nonlinear constraints cannot be met
	// 	4:	Major iteration limit is reached
	// 	6:	x does not satisfy the first-order optimality conditions, and no improved point for the merit function could be found during the final linesearch
	// 	7:	The function derivates returned by funcon and funobj appear to be incorrect
	// 	9:	An input parameter is invalid  

	// release memory
	delete []A; 
	delete []bl; 
	delete []bu; 
	delete []istate; 
	delete []cJac; 
	delete []clambda; 
	delete []R; 
	delete []x_raw; 
	delete []iw; 
	delete []w; 
	delete []c; 
	delete []g; 
	return error; 
} 

int CMSSM::Maximize_LogPosterior_CSMINWEL(double &log_posterior_optimal, TDenseVector &x_optimal, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob, const TDenseVector &x0) 
// Returns
// 	MODEL_NOT_PROPERLY_SET:	if state_equation_function, measurement_equation_function or transition_prob_function not properly set
// 	>=0:	inform code returned by npsol_
{
	if (prior_distr_function == NULL || CheckModelFunctions() != SUCCESS)
	{
		log_posterior_optimal = MINUS_INFINITY_LOCAL; 
		return MODEL_NOT_PROPERLY_SET; 
	}

	// Setting MinusLogLikelihood static variables
	MinusLogPosterior_CSMINWEL::model = this;
        MinusLogPosterior_CSMINWEL::y = y;
        MinusLogPosterior_CSMINWEL::z_0 = z_0;
        MinusLogPosterior_CSMINWEL::P_0 = P_0;
        MinusLogPosterior_CSMINWEL::initial_prob = initial_prob;

	// CSMINWEL
	size_t n = x0.dim;
        TDenseMatrix H(n,n,0.0);
        TDenseVector g(n,0.0);
        double *x_array, fh, crit;
        int itct, nit, fcount, retcodeh;
        x_array = new double[n];
        memcpy(x_array, x0.vector, sizeof(double)*n);
        crit = 1.0e-3;
        nit = 50;

	H.Identity(n);
        g.Zeros(n);
        dw_csminwel(&MinusLogPosterior_CSMINWEL::function,x_array,n,H.matrix,g.vector,NULL,&fh,crit,&itct,nit, &fcount,&retcodeh,NULL,NULL);
	
	x_optimal.Resize(n);
        for (unsigned int i=0; i<n; i++)
                x_optimal.SetElement(x_array[i],i);
        log_posterior_optimal = -fh;
        delete x_array;
        return retcodeh;
}
