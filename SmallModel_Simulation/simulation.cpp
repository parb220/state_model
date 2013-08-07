#include <iostream>
#include <ctime>
#include <getopt.h>
#include <cstdlib>
#include <mpi.h>
#include "dw_rand.h"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "CStorageHead.h"
#include "simulation.hpp"

using namespace std; 

int main(int argc, char **argv)
{
	// random number generator
	dw_initialize_generator(time(NULL));
	
	// Initialize MPI
	MPI_Init(&argc, &argv); 
	int my_rank, nNode; 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &nNode); 	

	static struct option long_options[] = 
	{
		{"data_file", required_argument, 0, 'd'}, 
		{"minus_infinity", required_argument, 0, 'm'},
		{"RunID", required_argument, 0, 'r'}, 
		{"Original", no_argument, 0, 'o'}, 
		{"nLevel", required_argument, 0, 'n'}, 
		{"HK", required_argument, 0, 'H'}, 	
		{"lHigh", required_argument, 0, 'L'}, 
		{"lLow", required_argument, 0, 'l'},
		{"H0", required_argument, 0, 'h'}, 
		{"T0", required_argument, 0, 't'}, 
		{"cFactor", required_argument, 0, 'c'}, 
		{"thin", required_argument, 0, 'i'}, 
		{"ndraws", required_argument, 0, 'w'},
		{"nInitial", required_argument, 0, 'I'}, 
		{"HillClimb", required_argument, 0, 'C'},
		{"TuningDone", no_argument, 0, 'T'},
		{0, 0, 0, 0}
	}; 
	int option_index = 0; 
	
	string data_file_name; 
	double minus_infinity = -1.0e30; 
	
	CEESParameter parameter; 
	parameter.storage_marker = 10000; 
	parameter.run_id = time(NULL); 
	parameter.storage_dir = getenv("HOME")+string("/state_model/mssm/result/");
	bool if_original = false;
	parameter.number_energy_level = 10; 
 	parameter.pee =0.3; 
	parameter.hk_1 = 200.0; 
	parameter.highest_level = parameter.number_energy_level -1; 
	parameter.lowest_level = 0; 
	parameter.h0 = 0.0; 
	parameter.t0 = 1.0; 
	parameter.c_factor = 1.0; 
	parameter.deposit_frequency = 50; 
	parameter.simulation_length = 10000; 
	size_t n_initial = (size_t)nNode; 
	size_t number_hill_climb = 0; 
	bool if_tuning_done = false; 

	while (1)
	{
		int c = getopt_long(argc, argv, "d:m:r:on:H:L:l:h:t:c:i:w:I:C:T", long_options, &option_index); 
		if (c == -1)
			break; 
		switch(c)
		{
			case 'd':
				data_file_name = string(optarg); break; 
			case 'm':
				minus_infinity = atof(optarg); break; 
			case 'r':
				parameter.run_id = atoi(optarg); break; 
			case 'o':
				if_original = true; break;
			case 'n':
				parameter.number_energy_level = atoi(optarg); break;  	
			case 'H':
				parameter.hk_1 = atof(optarg); break; 
			case 'L':
				parameter.highest_level = atoi(optarg); break; 
			case 'l':
				parameter.lowest_level = atoi(optarg); break; 			
			case 'h':
				parameter.h0 = atof(optarg); break; 
			case 't':
				parameter.t0 = atof(optarg); break; 
			case 'c':
				parameter.c_factor = atof(optarg); break;
			case 'i':
				parameter.deposit_frequency = atoi(optarg); break; 
	 		case 'w':
				parameter.simulation_length = atoi(optarg); break; 
			case 'I':
				n_initial = atoi(optarg); break; 
			case 'C':
				number_hill_climb = atoi(optarg); break; 
			case 'T':
				if_tuning_done = true; break; 
			default:
				break; 
		}
	}

	if (data_file_name.empty() ) 
	{
		cerr << "Usage: " << argv[0] << " -d data file -i minus infinity.\n"; 
		abort(); 
	}
	if (if_original)
	{
		parameter.number_energy_level =1; 
		parameter.highest_level = parameter.lowest_level = 0; 
	}
	parameter.SetEnergyBound(); 
	parameter.SetTemperature(); 
	int number_bin = (parameter.number_energy_level+1) * (parameter.number_energy_level+1); 
	CStorageHead storage(parameter.run_id, parameter.storage_marker, number_bin, parameter.storage_dir, my_rank); 

	if (my_rank == 0)
		master_deploying(nNode, if_tuning_done, number_hill_climb, n_initial, parameter, storage); 
	else 
		slave_computing(if_original, number_hill_climb, n_initial, data_file_name, minus_infinity, parameter, storage); 
}
