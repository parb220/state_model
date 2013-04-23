#include <string>
#include "CMSSM.hpp"

extern "C"
{
        void npsol_(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A, double *bl, double *bu, void *funcon(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate), void *funobj(int *mode, int *n, double *x, double *f, double *g, int *nstate), int *inform, int *iter, int *istate, double *c, double *cJac, double *clamda, double *f, double *g, double *R, double *x, int *iw, int *leniw, double *w, int *lenw);
        void npoptn_(char *, int);
}


using namespace std; 

class MinusLogLikelihood
{
public:
        static CMSSM *model;
        static vector<TDenseVector> y;
        static vector<TDenseVector> z_0;
        static vector<TDenseMatrix> P_0;
        static TDenseVector initial_prob;

        static void *function(int *mode, int*n, double *x, double *f, double *g, int *nstate);
};

CMSSM * MinusLogLikelihood::model; 
vector<TDenseVector> MinusLogLikelihood::y; 
vector<TDenseVector> MinusLogLikelihood::z_0; 
vector<TDenseMatrix> MinusLogLikelihood::P_0; 
TDenseVector MinusLogLikelihood::initial_prob; 

void * MinusLogLikelihood::function(int *mode, int *n, double *x_array, double *f, double *g, int *nstate)
{
        double minus_log_likelihood=1.0e30;

	// Make TDenseVector out of x_array
        TDenseVector x(*n); 
	for (unsigned int i=0; i<*n; i++)
		x.SetElement(x_array[i],i); 

	// Modify model according to x
	if ( !model->UpdateStateModelParameters(0,y,x) )
	{
		vector<vector<TDenseVector> > z_tm1;
        	vector<vector<TDenseMatrix> > P_tm1;
        	vector<TDenseVector> p_tm1;
       
		// Get initial values fo z_0 and P_0 from first few observations
        	size_t initial_period = 4;
        	vector<TDenseVector> sub_y(y.begin(), y.begin()+initial_period);
        	double log_likelihood;

		int kalman_error = model->KalmanFilter(log_likelihood, z_tm1, P_tm1, p_tm1, sub_y, z_0, P_0, initial_prob, x);
        	if (!kalman_error)
        	{
        		vector<TDenseVector> new_z_0 = z_tm1.back();
        		vector<TDenseMatrix> new_P_0 = P_tm1.back();
       
        		// Kalman filter 
        		sub_y = vector<TDenseVector>(y.begin()+initial_period, y.end());
			kalman_error = model->KalmanFilter(log_likelihood, z_tm1, P_tm1, p_tm1, sub_y, new_z_0, new_P_0, initial_prob, x);
                       if (!kalman_error)
        			minus_log_likelihood = -log_likelihood;
        	}
	*f = minus_log_likelihood; 
	}
}

int CMSSM::Minimize_MinusLogLikelihood(double &minus_log_likelihood_optimal, TDenseVector &x_optimal, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob, const TDenseVector &x0) 
// Returns
// 	-1:	if state_equation_function, measurement_equation_function or transition_prob_function not properly set
// 	>=0:	inform code returned by npsol_
{
	if (CheckModelFunctions() )
		return -1;
 
	int error =1; 

	// Setting MinumsLogLikelihood static variables
	MinusLogLikelihood::model = this; 
	MinusLogLikelihood::y = y; 
	MinusLogLikelihood::z_0 = z_0; 
	MinusLogLikelihood::P_0 = P_0; 
	MinusLogLikelihood::initial_prob = initial_prob; 

	const double INFINITE_BOUND = 1.0E20;
	const string COLD_START = string("Cold Start");
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
	
	npoptn_((char*)COLD_START.c_str(), COLD_START.length());
	npsol_(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, NULL, MinusLogLikelihood::function, &inform, &iter, istate, c, cJac, clambda, &f, g, R, x_raw, iw, &leniw, w, &lenw); 

	x_optimal.Resize(n); 
	for (unsigned int i=0; i<n; i++)
		x_optimal.SetElement(x_raw[i], i); 
	minus_log_likelihood_optimal = f; 
	error = inform; 

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
