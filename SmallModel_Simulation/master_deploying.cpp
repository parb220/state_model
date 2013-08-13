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
	double max_log_posterior, received_log_posterior; 
	if (number_hill_climb > 0)
	{
		vector<int> node_pool(nNode-1); 
		for (unsigned int i=1; i<nNode; i++)
			node_pool[i-1] =i; 
		max_log_posterior = DispatchHillClimbTask(node_pool, number_hill_climb, parameter, storage); 
		cout << "HillClimb: max_log_posterior: " << max_log_posterior << endl; 
	}
	
	// node_group
	vector<vector<int> >node_group(n_initial); 
	unsigned int i=1; 
	while (i<nNode)
	{
		for (unsigned int j=0; j<node_group.size() && i<nNode; j++)
		{
			node_group[j].push_back(i); 
			i++; 
		}
	}

	// tuning
	if(!if_tuning_done)
		max_log_posterior=DispatchTuneSimulation(node_group, parameter, storage); 	
	else 
	{
		for (int level=parameter.highest_level; level>=parameter.lowest_level; level--)
		{
			received_log_posterior=DispatchSimulation(node_group, parameter, storage, parameter.simulation_length, level, SIMULATION_TAG);
			if (level == parameter.highest_level)
				max_log_posterior = received_log_posterior; 
			else 
				max_log_posterior = max_log_posterior > received_log_posterior ? max_log_posterior : received_log_posterior; 
		}
	}
	cout << "Simulation: max_log_posterior: " << max_log_posterior << endl;  
	
	double *sMessage= new double [N_MESSAGE];
        for (unsigned int i=0; i<node_group.size(); i++)
                for (unsigned int j=0; j<node_group[i].size(); j++)
                        MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, node_group[i][j], END_TAG, MPI_COMM_WORLD);
        delete [] sMessage;
	return 0; 
}
