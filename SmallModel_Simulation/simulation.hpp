#ifndef _SIMULATION_HEADER_
#define _SIMULATION_HEADER_

int master_deploying(size_t nNode, bool if_tuning_done, size_t number_hill_climb, size_t n_initial, const CEESParameter &parameter, CStorageHead &storage);

int slave_computing(bool if_original, size_t number_hill_climb, size_t n_initial, const string &data_file_name, double minus_infinity, CEESParameter &parameter, CStorageHead &storage);

#endif
