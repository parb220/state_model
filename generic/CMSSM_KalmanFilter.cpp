#include <cmath>
#include "CMSSM_Error_Code.hpp"
#include "CMSSM.hpp"
#include "dw_math.h"

using namespace std; 

int CMSSM::KalmanFilter(double &log_likelihood, vector<vector<TDenseVector> > &z_tm1, vector<vector<TDenseMatrix> > &P_tm1, vector<TDenseVector > &p_tm1, const vector<TDenseVector> &y, const vector<TDenseVector> &z_0, const vector<TDenseMatrix> &P_0, const TDenseVector &initial_prob, const TDenseVector &x)
// Calculates
// 	log_likelihood = p (Y(T) | parameters)
//
// Returned values
// 	error_code:
// 		MODEL_NOT_PROPERLY_SET:	state_equation_function, measurement_equation_function or transition_prob_function not properly set
// 		SUCCESS:	success
//		ERROR_OCCURRED:	error occurred
// 
// Results:
//
// 	for k=0,1,...,nXi-1 and t=0,1,...,T-1
// 	z_tm1(t,k) = E( z(t) | Y(t-1), xi(t)=k ), 
// 	P_tm1(t,k) = variance( z(t) | Y(t-1), xi(t)=k )
// 	p_tm1(t,k) = p( xi(t)=k | Y(t-1) )
//
// Auxiliary values (not returned but used in calculations)
// 	for k=0,1,...,nXi-1, j=0,1,...,nZeta-1, and t=0,1,...,T-1
//	z_t(t,k) = E( z(t) | Y(t), xi(t)=k )
//	P_t(t,k) = variance( z(t) | Y(t), xi(t)=k )
//	N_tm1(t,k) = variance( y(t) | Y(t-1), xi(t) = k )
//	ey_tm1(t,k) = y(t) - E( y(t) | Y(t-1), xi(t) = k )
//	Iz_t(t,j) = E ( z(t) | Y(t), zeta(t) = j )
//	IP_t(t,j) = variance( z(t) | Y(t-1), zeta(t) = j ) ??? variance( z(t) | Y(t), zeta(t)=j )
//	p_t(t,k) = p( xi(t)=k | Y(t) )
//	log_conditional_likelihood(t,k) = p( y(t) | Y(t-1), xi(t) )
//
// Input values
// 	y(t):	nY-by-1 TDenseVector, measurment variables at time t, t=0, 1, ..., T-1
// 	z_0:	nZ-by-1 TDenseVector, E( z(0) | xi(0)=k ) k=0, 1, ..., nXi-1
// 	P_0:	nZ-by-nZ TDenseMatrix, variance( z(0) | xi(0)=k ) k=0, 1, ..., nXi-1
// 	initial_prob[k]:	p( xi(0) = k), k=0, 1, ..., nXi-1
//
// Valules inherent in CMSSM model
// 	nL:	number of lagged regimes that enter the measurement or state equations
// 	nTL:	total number of lagged regimes to track
// 	nS:	number of regimes
// 	nZ:	number of state variables
// 	nY:	number of measurement variables
// 	nU:	number of measurement noise variables
// 	nE:	number of state noise variables
// 	nNu:	nS^nL, number of regimes affecting coefficients
// 	nXi:	nS^(nTL+1), total number of regimes tracked
// 	nZeta:	nS^nTL
//
// Values inferred from input parameters and CMSSM inherent parameters
// 	T:	y.size(), number of time points for which measurments are available
{
	const double MINUS_INFINITY_LOCAL = -1.0e30; 
	const double CONSTANT = -0.918938533204673*nY; 
	
	if (CheckModelFunctions() != SUCCESS) 
	{
		log_likelihood = MINUS_INFINITY_LOCAL; 	
		return MODEL_NOT_PROPERLY_SET; 
	}

	if (UpdateStateModelParameters(0,y,x) != SUCCESS)
	{
		log_likelihood = MINUS_INFINITY_LOCAL; 
		return ERROR_OCCURRED; 
	}
	
	// Sizes and constants
	size_t T = y.size(); 	
	size_t loop_2 = (size_t)pow(nS,nL); 
	size_t loop_1 = (size_t)pow(nS,nTL-nL); 
	
	// output values
	z_tm1 = vector<vector<TDenseVector> >(T,vector<TDenseVector>(nXi) );  
	P_tm1 = vector<vector<TDenseMatrix> >(T,vector<TDenseMatrix>(nXi) );
	p_tm1 = vector<TDenseVector> (T); 

	// intermediate values
	vector<vector<TDenseVector> > z_t = vector<vector<TDenseVector> >(T,vector<TDenseVector>(nXi) ); 
	vector<vector<TDenseMatrix> > P_t = vector<vector<TDenseMatrix> >(T,vector<TDenseMatrix>(nXi) );  
	vector<vector<TDenseMatrix> > N_tm1 = vector<vector<TDenseMatrix> >(T,vector<TDenseMatrix>(nXi ) ); 
	vector<vector<TDenseVector> > ey_tm1 = vector<vector<TDenseVector> >(T,vector<TDenseVector>(nXi ) ); 
	vector<vector<TDenseVector> > Iz_t = vector<vector<TDenseVector> >(T,vector<TDenseVector>(nZeta ) ); 
	vector<vector<TDenseMatrix> > IP_t = vector<vector<TDenseMatrix> >(T,vector<TDenseMatrix>(nZeta) ); 
	vector<TDenseVector > p_t = vector<TDenseVector >(T); 
	vector<TDenseVector > log_conditional_likelihood = vector<TDenseVector>(T); 	


	for (unsigned int t=0; t<T; t++)
	{
		for (unsigned int i=0; i<nXi; i++)
		{
			z_tm1[t][i] = TDenseVector(nZ,0.0); 
			P_tm1[t][i] = TDenseMatrix(nZ,nZ,0.0); 
			z_t[t][i] = TDenseVector(nZ,0.0); 
			P_t[t][i] = TDenseMatrix(nZ,nZ,0.0); 
			N_tm1[t][i] = TDenseMatrix(nY,nY,0.0); 
			ey_tm1[t][i] = TDenseVector(nY,0.0); 
		}
		for (unsigned int i=0; i<nZeta; i++)
		{
			Iz_t[t][i] = TDenseVector(nZ,0.0); 
			IP_t[t][i] = TDenseMatrix(nZ,nZ,0.0); 
		}
		p_tm1[t] = TDenseVector(nXi,0.0); 
		p_t[t] = TDenseVector(nXi,0.0); 
		log_conditional_likelihood[t] = TDenseVector(nXi,0.0); 
	}

	// initialization
	for (unsigned int i=0; i<nXi; i++)
	{
		z_tm1[0][i] = z_0[i];
		P_tm1[0][i] = P_0[i]; 
	}
	
	// debugging OK tables
	vector<vector<bool> > OK_t(T,vector<bool>(nXi,false)); // OK_t(t,xi)=true if z_t(t,xi) and P_t(t,xi) have been properly computed.
	vector<vector<bool> > OK_tm1(T,vector<bool>(nXi,false)); // OK_tm1(t,xi)=true if z_tm1(t,xi), P_tm1(t,xi), ey_tm1(t,xi) and N_tm1(t,xi) have been properly computed
	vector<vector<bool> > OK_zeta(T,vector<bool>(nZeta,false)); // OK_tm1(t,xi)=true if Iz_t(t,zeta) and IP_t(t,zeta) have been properly computed

	OK_tm1[0]=vector<bool>(nXi,true); 
	log_likelihood=0.0; 
	TDenseMatrix Q; 
	for (unsigned int t=0; t<T; t++)
	{
		// p_tm1(t) = p( xi|Y(t-1) )
		if (t == 0)
			p_tm1[t] = initial_prob;
		else
		{
			if (GetTranstionProbMatrix(Q, t-1,y,x) != SUCCESS)
			{
				log_likelihood = MINUS_INFINITY_LOCAL;
                                return ERROR_OCCURRED;
			}
			p_tm1[t].Multiply(Q,p_t[t-1]);  
		}

		// === Kalman filter starts here ===========================================
		// integrate z_t(t-1,zeta) and P_t(t-1,zeta) if t>0
		if (t > 0)
		{
			for (unsigned int zeta=0; zeta<nZeta; zeta++)
			{
				double conditional_prob = 0.0; 
				for (unsigned int s=0; s<nS; s++)
				{
					unsigned int xi = zeta+nZeta*s; 
					if(p_t[t-1][xi] > 0)
					{
						if(!OK_t[t-1][xi])
						{
							log_likelihood = MINUS_INFINITY_LOCAL; 
							return ERROR_OCCURRED; 
						}
						conditional_prob += p_t[t-1][xi]; 
						Iz_t[t-1][zeta] = Iz_t[t-1][zeta] + z_t[t-1][xi]*p_t[t-1][xi]; 
						IP_t[t-1][zeta] = IP_t[t-1][zeta] + P_t[t-1][xi]*p_t[t-1][xi]; 
					}	
				}
				if (conditional_prob > 0.0)
				{
					Iz_t[t-1][zeta] = Iz_t[t-1][zeta] * (1.0/conditional_prob); 
					IP_t[t-1][zeta] = IP_t[t-1][zeta] * (1.0/conditional_prob);
					// force symmetric
					IP_t[t-1][zeta] = (IP_t[t-1][zeta] + Transpose(IP_t[t-1][zeta]) )*0.5; 
					OK_zeta[t-1][zeta] = true; 
				}
			}
		}

		for (unsigned int jj=0; jj<loop_1; jj++)
		{
			for (unsigned int ii=0; ii<loop_2; ii++)
			{
				unsigned int zeta = ii+loop_2*jj; 
				if (t==0 || OK_zeta[t-1][zeta])
				{
					for (unsigned s=0; s<nS; s++)
					{
						unsigned nu = s+nS*ii; 
						unsigned xi = s+nS*zeta; 

						// computes z_tm1[t][xi] and P_tm1[t][xi] if t>0
						if (t>0)
						{
							z_tm1[t][xi].Add(b[nu],F[nu]*Iz_t[t-1][zeta]);
							if (nE > 0)
								P_tm1[t][xi].Add(F[nu]*MultiplyTranspose(IP_t[t-1][zeta], F[nu]), V[nu]); 
							else 
								P_tm1[t][xi]=F[nu]*MultiplyTranspose(IP_t[t-1][zeta], F[nu]); 
							// force symmetric
							P_tm1[t][xi] = (P_tm1[t][xi]+Transpose(P_tm1[t][xi]) )*0.5;  
						}
						OK_tm1[t][xi] = true; 
						
						// temporary matrix, K=cov( z(t), y(t)| Y(t-1), xi(t) = xi)
						TDenseMatrix K = MultiplyTranspose(P_tm1[t][xi], H[nu]); 

						// compute y_tm1[t][xi], ey_tm1[t][xi] and N_tm1[t][xi]
						ey_tm1[t][xi] =  y[t] - a[nu] - H[nu]*z_tm1[t][xi];
						if (nU > 0)
							N_tm1[t][xi].Add(H[nu]*K, R[nu]); 
						else 
							N_tm1[t][xi] = H[nu]*K; 

						// force symmetry
						N_tm1[t][xi] = (N_tm1[t][xi]+Transpose(N_tm1[t][xi]) ) *0.5; 

						// Theoretically, N_tm1[t][xi], P_tm1[t][xi], IP_t[t-1][zeta] and P_t[t][xi] all
						// represent covariance matrices, and are therefore all semi positive. Consequently,
						// ey_tm1[t][xi]*Inverse(N_tm1[t][xi])*ey_tm1[t][xi] should >=0;
						// However, computationally, because 
						//	P_t[t][xi] = P_tm1[t][xi]-K* Inverse(N_tm1[t][xi]) 
						//	ey_tm1[t][xi]*Inverse(N_tm1[t][xi])*ey_tm1[t][xi]	
						// where the inverse are both based on LU decomposition and is not numically stable.
						// So if ey_tm1[t][xi]*Inverse(N_tm1[t][xi])*ey_tm1[t][xi] renders a negative result,
						// it should be abondoned. 

						// Cholesky decomposition of N_tm1[t][xi]
						TDenseMatrix Cholesky_T; 
						try 
						{	Cholesky_T = Cholesky(N_tm1[t][xi]); }
						catch (...)
						{
							log_conditional_likelihood[t].SetElement(MINUS_INFINITY_LOCAL, xi); 
                                                        OK_t[t][xi] = false;
                                                        continue;
						}
						double log_abs_determinant = N_tm1[t][xi].LogAbsDeterminant_Cholesky(Cholesky_T);
						// double log_abs_determinant = N_tm1[t][xi].LogAbsDeterminant(); 
						TDenseVector tempV; 
						try 
						{ 	tempV=LeftSolve_Cholesky(Cholesky_T, ey_tm1[t][xi]); }
							// tempV = LeftSolve(N_tm1[t][xi], ey_tm1[t][xi]);  }
						catch (...)
						{
							log_conditional_likelihood[t].SetElement(MINUS_INFINITY_LOCAL, xi); 
							OK_t[t][xi] = false;
							continue;  
						}
						double inner_product = InnerProduct(ey_tm1[t][xi],tempV); 
						if (inner_product >= 0)
						{
							log_conditional_likelihood[t].SetElement(CONSTANT - 0.5*(log_abs_determinant+inner_product), xi);
							// computes z_t[t][xi] and P[t][xi]
							z_t[t][xi].Add(z_tm1[t][xi], K*tempV);
							// If the program can reach here, then N_tm1[t][xi] is inversible
							// Therefore tempM should not incurr errors
							TDenseMatrix tempM = LeftSolve_Cholesky(Cholesky_T,Transpose(K)); 
							// TDenseMatrix tempM = LeftSolve(N_tm1[t][xi], Transpose(K)); 
							P_t[t][xi].Subtract(P_tm1[t][xi], K*tempM ); 
							// force symmetry
							P_t[t][xi] = 0.5*(P_t[t][xi] + Transpose(P_t[t][xi]) ); 
							OK_t[t][xi] = true; 
						}
						else
						{
							log_conditional_likelihood[t].SetElement(MINUS_INFINITY_LOCAL, xi);
							OK_t[t][xi] = false; 	
						}
					}
				}
				else 
				{
					for (unsigned int s=0; s<nS; s++)
					{
						unsigned int xi = s + nS*zeta; 
						log_conditional_likelihood[t].SetElement(MINUS_INFINITY_LOCAL, xi); 
					}
				}
			}
		}

		// == Hamilton filter =====
		double scale = MINUS_INFINITY_LOCAL; 
		for (unsigned int xi=0; xi<nXi; xi++)
		{
			if (p_tm1[t][xi] > 0.0)
			{
				scale = AddScaledLogs(1.0,scale,p_tm1[t][xi],log_conditional_likelihood[t][xi]); 
				p_t[t].SetElement(log( p_tm1[t][xi] ) + log_conditional_likelihood[t][xi], xi); 
			}
			else 
				p_t[t].SetElement(MINUS_INFINITY_LOCAL, xi); 
		}

		// update log liklihoood
		log_likelihood += scale; 

		// update p_t
		for (unsigned int xi=0; xi<nXi; xi++)
			p_t[t].SetElement(exp( p_t[t][xi] - scale), xi);  
		// End Hamilton filter

		// Get State space parameter
		if (UpdateStateModelParameters(t+1, y, x) != SUCCESS)
		{
			log_likelihood == MINUS_INFINITY_LOCAL; 
			return ERROR_OCCURRED; 
		} 
	}
	return SUCCESS;
}
