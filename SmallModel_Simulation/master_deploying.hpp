#ifndef _MASTER_DEPLOYING_HEADER_
#define _MASTER_DEPLOYING_HEADER_

int DispatchHillClimbTask(const vector<int> &node_pool, size_t number_hill_climb, const CEESParameter &parameter, CStorageHead &storage); 

double DispatchTuneSimulation(const vector<vector<int> > &node_group, const CEESParameter &parameter, CStorageHead &storage);

double DispatchSimulation(const vector<vector<int> > &nodeGroup, const CEESParameter &parameter, CStorageHead &storage, size_t simulation_length, int level, int tag);

#endif
