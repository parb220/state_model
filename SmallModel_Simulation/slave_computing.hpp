#ifndef _SLAVE_COMPUTING_HEADER
#define _SLAVE_COMPUTING_HEADER

bool ReadOutData(const string &data_file_name, vector<TDenseVector> &y1st, vector<TDenseVector> &y2nd, vector<TDenseVector> &y); 

void SetUpModel(size_t nFree, size_t nY, TDenseVector &fixed_parameter, CMSSM_test &model_1st, CMSSM_test_2nd &model_2nd, CMSSM_test &model_all); 

void ExecuteHillClimbTask(size_t nFree, const TDenseVector &fixed_parameter, CMSSM_test &model_1st, CMSSM_test_2nd &model_2nd, CMSSM_test &model_all, const vector<TDenseVector> &y1st, const vector<TDenseVector> &y2nd, const vector<TDenseVector> &y, const CEESParameter &parameter, CStorageHead &storage); 

bool ExecuteTuningTask_BeforeSimulation(size_t period, size_t max_period, CEquiEnergy_CMSSM_test &model, CStorageHead &storage, const CEESParameter &parameter, int group_index, size_t pool_size); 

bool ExecuteTuningTask_AfterSimulation(size_t period, size_t max_period, CEquiEnergy_CMSSM_test &model, const CEESParameter &parameter, int group_index); 

bool ExecuteSimulationTask(bool if_within, bool if_write_sample_file, bool if_storage, CEquiEnergy_CMSSM_test &model, CStorageHead &storage, const CEESParameter &parameter, int my_rank, int group_index, size_t pool_size, int message_tag); 

#endif

