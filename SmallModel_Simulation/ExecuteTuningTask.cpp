#include <sstream>
#include <cstdlib>
#include <fstream>
#include "dw_dense_matrix.hpp"
#include "CEquiEnergy_CMSSM_test.hpp"
#include "CStorageHead.h"
#include "CEESParameter.h"
#include "CMetropolis.h"
#include "CSampleIDWeight.h"
#include "storage_parameter.h"

using namespace std; 
bool ExecuteTuningTask_BeforeSimulation(size_t period, size_t max_period, CEquiEnergy_CMSSM_test &model, int group_index, size_t pool_size)
{
	// start point
	model.storage->RestoreForFetch(model.energy_level+1);
	if (model.storage->empty(model.energy_level+1) || !model.Initialize_RandomlyPickFrom_K_BestSample(pool_size, model.energy_level+1)) 
		return false; 
	model.storage->ClearDepositDrawHistory(model.energy_level+1);

	// save the start point
	stringstream convert;
	convert.str(string()); 
	convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << model.energy_level << "." << group_index; 
	string start_point_file = model.parameter->storage_dir + convert.str();
	ofstream output_file;
	output_file.open(start_point_file.c_str(), ios::binary|ios::out); 
	if (!output_file)
		return false; 
	else 
	{
		CSampleIDWeight x_complete; 
		vector<int> locs_variable=model.target_model->locs_variable(); 
		vector<int> locs_constant=model.target_model->locs_constant(); 
		x_complete.data.Resize(locs_variable.size()+locs_constant.size()); 
		x_complete.data.SetSubVector(locs_variable, model.current_sample.data); 
		x_complete.data.SetSubVector(locs_constant, model.target_model->GetConstantPart()); 	
	
		write(output_file, &x_complete); 
	}
	output_file.close(); 
	
	// tuning 
        convert.str(string());
       	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_1ST << model.energy_level << "." << group_index;
	string block_file_name = model.parameter->storage_dir + convert.str();	
	
	if (model.energy_level == model.parameter->number_energy_level-1)
	{
		if (!model.metropolis->AdaptiveBeforeSimulation_OnePass(model.current_sample, period, max_period, block_file_name))
		{
                	cerr << "CMetropolis::AdaptiveBeforeSimulation() : Error in writing " << block_file_name << endl;
                       	abort();
                }
	}
	else 
	{
		// based on samples of the previous level (not group specific) to estimate covariance matrix
		convert.str(string());
                convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level+1;
                string sample_file_name = model.parameter->storage_dir + convert.str();
                if (!model.metropolis->AdaptiveAfterSimulation_OnePass(model.current_sample, period, max_period, sample_file_name, block_file_name))
                {
                	cerr << "CMetroplis::AdaptiveAfterSimulation() : Error in reading " << sample_file_name << " or writing " << block_file_name << endl;
               	abort();
                }
	}
	return true; 
}


bool ExecuteTuningTask_AfterSimulation(size_t period, size_t max_period, CEquiEnergy_CMSSM_test &model, int group_index)
{
	// start point
	stringstream convert; 
	convert << model.parameter->run_id << "/" << model.parameter->run_id << START_POINT << model.energy_level << "." << group_index; 
	string start_point_file = model.parameter->storage_dir + convert.str(); 
	if (!model.InitializeFromFile(start_point_file))
		return false;  

	convert.str(string());
       	convert << model.parameter->run_id << "/" << model.parameter->run_id << BLOCK_2ND << model.energy_level << "." << group_index; 
        string block_file_name = model.parameter->storage_dir + convert.str();
	// sample file: based on the samples of the current level (group specific)
	convert.str(string());
        convert << model.parameter->run_id << "/" << model.parameter->run_id << VARIANCE_SAMPLE_FILE_TAG << model.energy_level << "." << group_index;	
	string sample_file_name = model.parameter->storage_dir + convert.str();

	if (!model.metropolis->AdaptiveAfterSimulation_OnePass(model.current_sample, period, max_period, sample_file_name, block_file_name) )
	{
		cerr << "CMetroplis::AdaptiveAfterSimulation() : Error in reading " << sample_file_name << " or writing " << block_file_name << endl;
		abort();
	}
	return true; 
}
