#include <iostream>
#include <vector>
#include <mpi.h>
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "master_deploying.hpp"
#include "mpi_parameter.h"

using namespace std; 

int master_deploying(size_t nNode, bool if_tuning_done, size_t number_hill_climb, size_t n_initial, const CEESParameter &parameter, CStorageHead &storage)
{
	if (number_hill_climb > 0)
	{
		vector<int> node_pool(nNode-1); 
		for (int i=1; i<nNode; i++)
			node_pool[i-1] =i; 
		DispatchHillClimbTask(node_pool, number_hill_climb, parameter, storage); 
	}
	
	// node_group
	vector<vector<int> >node_group(n_initial); 
	int i=1; 
	while (i<nNode)
	{
		for (int j=0; j<node_group.size() && i<nNode; j++)
		{
			node_group[j].push_back(i); 
			i++; 
		}
	}

	// tuning
	if(!if_tuning_done)
		DispatchTuneSimulation(node_group, parameter, storage); 	
	else 
	{
		for (int level=parameter.highest_level; level>=parameter.lowest_level; level--)
			DispatchSimulation(node_group, parameter, storage, parameter.simulation_length, level, SIMULATION_TAG);
	}
	
	double *sMessage= new double [N_MESSAGE];
        for (int i=0; i<node_group.size(); i++)
                for (int j=0; j<node_group[i].size(); j++)
                        MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, node_group[i][j], END_TAG, MPI_COMM_WORLD);
        delete [] sMessage;
	return 0; 
}
