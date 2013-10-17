#include <vector>
#include <ctime>
#include <mpi.h>
#include <cstdlib>
#include "dw_dense_matrix.hpp"
#include "CMSSM_test.hpp"
#include "CMSSM_test_2nd.hpp"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "CEquiEnergy_CMSSM_test.hpp"
#include "CMetropolis.h"
#include "mpi_parameter.h"
#include "slave_computing.hpp"

using namespace std; 
bool GetCommunicationParameter(const double *rPackage, size_t package_size, CEESParameter *parameter)
{
        parameter->simulation_length = (size_t)(rPackage[LENGTH_INDEX]);
        parameter->burn_in_length = (size_t)(rPackage[BURN_INDEX]);
        parameter->thin = (size_t)(rPackage[thin_INDEX]);
        parameter->THIN = (size_t)(rPackage[THIN_INDEX]);

        parameter->SetTemperature();
        return true;
}

int slave_computing(bool if_original, size_t number_hill_climb, size_t n_initial, const string &data_file_name, double minus_infinity, CEESParameter &parameter, CStorageHead &storage)
{
	// reading data
	vector<TDenseVector> y1st, y2nd, y; 
	if (!ReadOutData(data_file_name, y1st, y2nd, y))
	{
		cerr << "slave_computing : ReadOutData() : error occurred in LoadData()\n"; 
		abort(); 
	}
	size_t nY = y[0].dim; 

	// CMSSM models 
	CMSSM::MINUS_INFINITY_LOCAL = minus_infinity;
	size_t nFree = 21+8+2;  // 21 (15+6): model parameters; 6=2X3 (Delta in Dan's sunspont notation) + 2=2x1 (gamma in Dan's sunspot notation): sunspot parameters (with 2 endogeneous errors and 3 fundamental shocks); 2 probabilites of staying in ZLB

	TDenseVector fixed_parameter; 	//fixed_parameters for (RationalExpectationFunction, StateEquationFunction, MeasurementEquationFunction) and TransitionProbMatrix
	CMSSM_test model_1st;
	CMSSM_test_2nd model_2nd; 
	CMSSM_test model_all; 
	SetUpModel(nFree, nY, fixed_parameter, model_1st, model_2nd, model_all); 
	
	// HillClimb
	if (number_hill_climb > 0)
	{
		ExecuteHillClimbTask(nFree, fixed_parameter, model_1st, model_2nd, model_all, y1st, y2nd, y, parameter, storage); 
	}

	// EquiEnergyModel 
	CEquiEnergy_CMSSM_test simulation_model; 
	simulation_model.target_model = &model_all; 
	simulation_model.timer_when_started = time(NULL); 
	if (if_original)
		simulation_model.if_bounded = false; 
	CMetropolis metropolis(&simulation_model); 
	simulation_model.metropolis = &metropolis; 
	simulation_model.parameter = &parameter;
	simulation_model.storage = &storage; 
	simulation_model.y = y; 
	simulation_model.initial_prob=TDenseVector(model_1st.nXi,0.0); 
	simulation_model.initial_prob[0] = 1.0; 
	simulation_model.z_0 = vector<TDenseVector>(model_1st.nXi); 
	simulation_model.P_0 = vector<TDenseMatrix>(model_1st.nXi); 
	for (int i=0; i<model_1st.nXi; i++)
	{
		simulation_model.z_0[i].Zeros(model_1st.nZ); 
		simulation_model.P_0[i] = Identity(model_1st.nZ)*100.0; 
	}
	
	// tuning, simulation	
	int my_rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	MPI_Status status; 

	double *rPackage=new double[N_MESSAGE], *sPackage=new double[N_MESSAGE]; 
	int group_index; 
	bool if_within, if_write_file, if_storage; 
	while (1)
	{
		MPI_Recv(rPackage, N_MESSAGE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
		if (status.MPI_TAG == END_TAG)
		{
			delete []sPackage; 
			delete []rPackage; 
			exit(0); 
		}
		simulation_model.energy_level = (int)rPackage[LEVEL_INDEX]; 
		group_index = (int)rPackage[GROUP_INDEX]; 
		if (!GetCommunicationParameter(rPackage, N_MESSAGE, simulation_model.parameter) )
		{
			cerr << "GetCommunicationParameter() : Error occurred.\n"; 
			abort(); 
		}
		simulation_model.t_bound = simulation_model.parameter->t[simulation_model.energy_level]; 		
		if (status.MPI_TAG == TUNE_TAG_BEFORE_SIMULATION || status.MPI_TAG == TUNE_TAG_AFTER_SIMULATION)
		{
			size_t period = 20; 
			size_t max_period = 16*period; 
			if (status.MPI_TAG == TUNE_TAG_BEFORE_SIMULATION)
			{
				if (!ExecuteTuningTask_BeforeSimulation(period, max_period, simulation_model, group_index, n_initial) )
				{
					cerr << "ExecuteTuningTask_BeforeSimulation() : Error occurred :: sample file reading or block_file writing or start_tune_point writing error.\n";
                                	abort();
				}
			}
			else 
			{
				if (!ExecuteTuningTask_AfterSimulation(period, max_period, simulation_model, group_index) )
				{
					cerr << "ExecuteTuningTask_AfterSimulation() : Error occurred :: start_tune_point file reading or sample file reading or block_file writing error.\n";
                                        abort();
				}
			}	
		}
		else if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST || status.MPI_TAG == TUNE_TAG_SIMULATION_SECOND || status.MPI_TAG == SIMULATION_TAG)
		{
			if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST)
                                if_within = true;
                        else
                                if_within = false;

                        if (status.MPI_TAG == TUNE_TAG_SIMULATION_FIRST || status.MPI_TAG == TUNE_TAG_SIMULATION_SECOND)
                                if_write_file = true;
                        else
                                if_write_file = false;

                        if (status.MPI_TAG == TUNE_TAG_SIMULATION_SECOND || status.MPI_TAG == SIMULATION_TAG)
                                if_storage = true;
                        else
                                if_storage = false;

                        if (!ExecuteSimulationTask(if_within, if_write_file, if_storage, simulation_model,  my_rank, group_index, n_initial, status.MPI_TAG) )
			{
				cerr << "ExecuteSimulationTask : Error in simulation.\n"; 
				abort(); 
			}
		}
		MPI_Send(sPackage,N_MESSAGE,MPI_DOUBLE, 0, status.MPI_TAG, MPI_COMM_WORLD); 
	}
}

