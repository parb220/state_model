#ifndef _SLAVE_COMPUTING_HEADER
#define _SLAVE_COMPUTING_HEADER

bool ReadOutData(const string &data_file_name, vector<TDenseVector> &y1st, vector<TDenseVector> &y2nd, vector<TDenseVector> &y); 

void SetUpModel(size_t nY, TDenseVector &fixed_parameter, vector<int> &locs_x1, vector<int> &locs_x2, vector<int> &locs_xall, CMSSM_test_1st &model_1st, CMSSM_test_2nd &model_2nd, CMSSM_test_1st &model_all); 

double ExecuteHillClimbTask(size_t nFree, const TDenseVector &fixed_parameter, const vector<int> &locs_1st, const vector<int> &locs_2nd, const vector<int> &locs_all, CMSSM_test_1st &model_1st, CMSSM_test_2nd &model_2nd, CMSSM_test_1st &model_all, const vector<TDenseVector> &y1st, const vector<TDenseVector> &y2nd, const vector<TDenseVector> &y, const CEESParameter &parameter, CStorageHead &storage); 

bool ExecuteTuningTask_BeforeSimulation(size_t period, size_t max_period, CEquiEnergy_CMSSM &model, CStorageHead &storage, const CEESParameter &parameter, unsigned int group_index, size_t pool_size); 

bool ExecuteTuningTask_AfterSimulation(size_t period, size_t max_period, CEquiEnergy_CMSSM &model, const CEESParameter &parameter, unsigned int group_index); 

bool ExecuteSimulationTask(double &max_log_posterior, bool if_within, bool if_write_sample_file, bool if_storage, CEquiEnergy_CMSSM &model, CStorageHead &storage, const CEESParameter &parameter, unsigned int my_rank, unsigned int group_index, size_t pool_size, int message_tag); 

#endif

