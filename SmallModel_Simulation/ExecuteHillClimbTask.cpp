#include <vector>
#include <mpi.h>
#include "dw_dense_matrix.hpp"
#include "CMSSM.hpp" 
#include "CMSSM_Error_Code.hpp"
#include "CMSSM_test_1st.hpp"
#include "CMSSM_test_2nd.hpp"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"

using namespace std; 

TDenseVector InitializeParameter(size_t n, const TDenseVector &fixed_parameter);

int HillClimb(CSampleIDWeight &optimal, const TDenseVector &x0, CMSSM_test_1st &model_1st, CMSSM_test_2nd &model_2nd, CMSSM_test_1st & model_all, const vector<int> &locs_x1, const vector<int> &locs_x2, const vector<int> &locs_all, const vector<TDenseVector> &y1st, const vector<TDenseVector> &y2nd, const vector<TDenseVector> &y); 

double ExecuteHillClimbTask(size_t nFree, const TDenseVector &fixed_parameter, const vector<int> &locs_x1, const vector<int> &locs_x2, const vector<int> &locs_xall, CMSSM_test_1st &model_1st, CMSSM_test_2nd &model_2nd, CMSSM_test_1st &model_all, const vector<TDenseVector> &y1st, const vector<TDenseVector> &y2nd, const vector<TDenseVector> &y, const CEESParameter &parameter, CStorageHead &storage)
{
	double *rPackage = new double[N_MESSAGE]; 
	MPI_Status status; 

	MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status); 
	unsigned int level = (unsigned int)rPackage[LEVEL_INDEX]; 
	storage.restore(parameter.BinIndex_Start(level), parameter.BinIndex_End(level)); 
	size_t nSolution = (size_t)rPackage[LENGTH_INDEX]; 
	
	// Hill Climb
	TDenseVector x0;   
	CSampleIDWeight sample; 
	double max_log_posterior = CMSSM::MINUS_INFINITY_LOCAL; 
	for (unsigned int i=0; i<nSolution; i++)
	{
		x0= InitializeParameter(nFree, fixed_parameter);

		if (HillClimb(sample, x0, model_1st, model_2nd, model_all, locs_x1, locs_x2, locs_xall, y1st, y2nd, y) == SUCCESS)
		{
			int bIndex = parameter.BinIndex(-sample.weight, level); 
			storage.DepositSample(bIndex, sample); 
			max_log_posterior = max_log_posterior > sample.weight ? max_log_posterior : sample.weight; 
		}
	}
	storage.finalize(parameter.BinIndex_Start(level), parameter.BinIndex_End(level)); 
	storage.ClearDepositDrawHistory(parameter.BinIndex_Start(level), parameter.BinIndex_End(level)); 

	double *sPackage = new double[N_MESSAGE]; 
	sPackage[H0_INDEX] = max_log_posterior; 
	MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, 0, HILL_CLIMB_TAG, MPI_COMM_WORLD); 
	delete [] rPackage; 
	delete [] sPackage; 
	return max_log_posterior; 
}
