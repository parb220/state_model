#include <vector>
#include <cmath>
#include <mpi.h>
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "mpi_parameter.h"

using namespace std; 

double DispatchHillClimbTask(const vector<int> &node_pool, size_t number_hill_climb, const CEESParameter &parameter, CStorageHead &storage)
{
	size_t n_solution_per_node = (size_t)ceil((double)(number_hill_climb)/(double)node_pool.size());
	double *sPackage = new double[N_MESSAGE]; 
	sPackage[LENGTH_INDEX] = n_solution_per_node; 
	sPackage[LEVEL_INDEX] = parameter.number_energy_level; 

	for (unsigned int i=0; i<node_pool.size(); i++)
		MPI_Send(sPackage, N_MESSAGE, MPI_DOUBLE, node_pool[i], HILL_CLIMB_TAG, MPI_COMM_WORLD); 
	delete [] sPackage; 

	int success=0; 
	MPI_Status status; 
	double *rPackage = new double[N_MESSAGE]; 
	double max_log_posterior; 
	for (unsigned int i=0; i<node_pool.size(); i++)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, MPI_ANY_SOURCE, HILL_CLIMB_TAG, MPI_COMM_WORLD, &status); 
		if (i == 0)
			max_log_posterior = rPackage[H0_INDEX];
		else 
			max_log_posterior = max_log_posterior > rPackage[H0_INDEX] ? max_log_posterior : rPackage[H0_INDEX]; 
	}
	delete [] rPackage; 

	storage.consolidate(parameter.BinIndex_Start(parameter.number_energy_level),parameter.BinIndex_End(parameter.number_energy_level)); 
	return max_log_posterior; 
}
