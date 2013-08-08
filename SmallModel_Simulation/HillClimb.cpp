#include <vector>
#include "dw_dense_matrix.hpp"
#include "CMSSM_Error_Code.hpp"
#include "CMSSM_test_1st.hpp"
#include "CMSSM_test_2nd.hpp"
#include "CSampleIDWeight.h"

using namespace std; 

int MakeLbUb_test(TDenseVector &lb, TDenseVector &ub, size_t bDim, int choice=1);

int RandomInit_test_1st(TDenseVector &x, const vector<int> &pos); 
int RandomInit_test_2nd(TDenseVector &x, const vector<int> &pos);

int HillClimb(CSampleIDWeight &optimal, const TDenseVector &x0, CMSSM_test_1st &model_1st, CMSSM_test_2nd &model_2nd, CMSSM_test_1st & model_all, const vector<int> &locs_x1, const vector<int> &locs_x2, const vector<int> &locs_xall, const vector<TDenseVector> &y1st, const vector<TDenseVector> &y2nd, const vector<TDenseVector> &y)
{
	// 1st regime (s_t*, s_{t-1}*) = (1*, 1*)
	// Starts in state one with probability one
	TDenseVector initial_prob_1st(model_1st.nXi, 0.0);
	initial_prob_1st(0) = 1.0;
	// Initial values for initial Kalman filter
	vector<TDenseVector> z0_1st(model_1st.nXi);
	vector<TDenseMatrix> P0_1st(model_1st.nXi);
	for (unsigned int i=0; i<model_1st.nXi; i++)
	{
		z0_1st[i].Zeros(model_1st.nZ);
        	P0_1st[i]=Identity(model_1st.nZ)*100.0;
	}

	// 2nd regime
	TDenseVector initial_prob_2nd(model_2nd.nXi, 0.0);
	initial_prob_2nd(10) = 1.0;     // 10: the position for (3*, 3*) (staying in the ZLB regime)
	vector<TDenseVector> z0_2nd(model_2nd.nXi);
	vector<TDenseMatrix> P0_2nd(model_2nd.nXi);

	// variables
	double log_likelihood_1st;
	double log_posterior_1st, log_posterior_2nd, log_posterior_all;
	vector<TDenseVector> z_tm1_last_1st;
	vector<TDenseMatrix> P_tm1_last_1st;
	TDenseVector p_tm1_last_1st;

	TDenseVector lb, ub; 
	TDenseVector x_input, initial_x1, initial_x2, initial_xall; 
	TDenseVector xOptimal_1st, xOptimal_2nd, xOptimal_all; 
	double lpOptimal_1st, lpOptimal_2nd, lpOptimal_all; 
	optimal.calculated = false; 
	optimal.weight = CMSSM::MINUS_INFINITY_LOCAL; 
	
	size_t n_try = 20; 
	bool bad = false; 
	
	x_input = x0; 
	for (unsigned int i_try=0; i_try<n_try; i_try++)
	{
		bad = false; 
		
		// 1st regime
		// initialization 
		initial_x1 = x_input; 
		if(i_try>0 && RandomInit_test_1st(initial_x1, locs_x1) != SUCCESS)
			bad = true; 
		// bound
		if(!bad && MakeLbUb_test(lb, ub, initial_x1.dim) != SUCCESS)
			bad = true; 
		// log_posterior for initial_x1
		if(!bad && (model_1st.LogPosterior(log_posterior_1st, initial_x1, y1st, z0_1st, P0_1st, initial_prob_1st) != SUCCESS || log_posterior_1st<=CMSSM::MINUS_INFINITY_LOCAL))
			bad = true; 
		// optimizing
		// if(!bad && model_1st.Maximize_LogPosterior_NPSOL(lpOptimal_1st,xOptimal_1st,y1st,z0_1st,P0_1st,initial_prob_1st,initial_x1, ub,lb) == MODEL_NOT_PROPERLY_SET)
		if (!bad && model_1st.Maximize_LogPosterior_CSMINWEL(lpOptimal_1st,xOptimal_1st,y1st,z0_1st,P0_1st,initial_prob_1st,initial_x1) == MODEL_NOT_PROPERLY_SET)
			bad = true; 
		// update x_input(locs_x1)
		if(!bad)
		{
			for (unsigned int j=0; j<locs_x1.size(); j++)
				x_input(locs_x1[j]) = xOptimal_1st(locs_x1[j]);
		}

		// 2nd regime
		// obtain z_tm1_last_1st and P_tm1_last_1st to intialize z_02nd and P0_2nd
		if(!bad && model_1st.LogLikelihood(log_likelihood_1st, z_tm1_last_1st, P_tm1_last_1st, p_tm1_last_1st, x0, y1st, z0_1st, P0_1st, initial_prob_1st) != SUCCESS)
			bad = true; 
		// initialize
		if(!bad)
		{
			for (unsigned int j=0; j<model_2nd.nXi; j++)
                	{
                		z0_2nd[j].CopyContent(z_tm1_last_1st[0]);
                		P0_2nd[j].CopyContent(P_tm1_last_1st[0]);
                	}
                	initial_x2 = x_input;
			if(i_try>0 && RandomInit_test_2nd(initial_x2, locs_x2) != SUCCESS)
				bad = true; 
		}
		// bound
		if(!bad && MakeLbUb_test(lb,ub,initial_x2.dim) != SUCCESS )
			bad = true; 
		// log-posterior for initial_x2
		if(!bad && (model_2nd.LogPosterior(log_posterior_2nd, initial_x2, y2nd, z0_2nd, P0_2nd, initial_prob_2nd) != SUCCESS || log_posterior_2nd <= CMSSM::MINUS_INFINITY_LOCAL))
			bad = true; 
		// optimizing
		// if(!bad && model_2nd.Maximize_LogPosterior_NPSOL(lpOptimal_2nd,xOptimal_2nd,y2nd,z0_2nd,P0_2nd,initial_prob_2nd,initial_x2, ub,lb) == MODEL_NOT_PROPERLY_SET)
		if (!bad && model_2nd.Maximize_LogPosterior_CSMINWEL(lpOptimal_2nd,xOptimal_2nd,y2nd,z0_2nd,P0_2nd,initial_prob_2nd,initial_x2) == MODEL_NOT_PROPERLY_SET) 
			bad = true; 
		// Update x_input(locs_x2)
		if(!bad)
                {
                	for (unsigned int j=0; j<locs_x2.size(); j++)
                		x_input(locs_x2[j]) = xOptimal_2nd(locs_x2[j]);                  
                }
		
		// all parameters
		// initialize and bound
		if(!bad)
		{
			initial_xall = x_input;
			if(!bad && MakeLbUb_test(lb,ub,initial_xall.dim) != SUCCESS )
				bad = true; 
		}
		// log_posterior for initial_xall
		if(!bad && (model_all.LogPosterior(log_posterior_all, initial_xall, y, z0_1st, P0_1st, initial_prob_1st) != SUCCESS || log_posterior_all <= CMSSM::MINUS_INFINITY_LOCAL) )
			bad = true; 

		// optimize
		// if(!bad && model_all.Maximize_LogPosterior_NPSOL(lpOptimal_all,xOptimal_all,y,z0_1st,P0_1st,initial_prob_1st,initial_xall,ub, lb) == MODEL_NOT_PROPERLY_SET)
		if(!bad && model_all.Maximize_LogPosterior_CSMINWEL(lpOptimal_all,xOptimal_all,y,z0_1st,P0_1st,initial_prob_1st,initial_xall) == MODEL_NOT_PROPERLY_SET)
			bad = true; 
		if (!bad)
		{
			// update x_input(locs_xall)
			for (unsigned int j=0; j<locs_xall.size(); j++)
				x_input(locs_xall[j]) = xOptimal_all(locs_xall[j]);
			// log_posterior foa final solution
			if (model_all.LogPosterior(log_posterior_all, x_input, y, z0_1st, P0_1st, initial_prob_1st) != SUCCESS)
				bad = true; 
		}
		if (!bad)
		{
			if (log_posterior_all > optimal.weight)
			{
				optimal.calculated = true; 
				optimal.weight = log_posterior_all; 
				optimal.data.CopyContent(x_input); 
			}
		}		
	}

	if (optimal.weight <= CMSSM::MINUS_INFINITY_LOCAL)
		return ERROR_OCCURRED; 
	else 
		return SUCCESS; 
}

