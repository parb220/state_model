#include <iostream>
#include <sstream>
#include <ctime>
#include <getopt.h>
#include <cstdlib>
#include <mpi.h>
#include "dw_rand.h"
#include "CSampleIDWeight.h"
#include "CEESParameter.h"
#include "mpi_parameter.h"
#include "CStorageHead.h"
#include "simulation.hpp"

using namespace std; 

int main(int argc, char **argv)
{
	// random number generator
	// dw_initialize_generator(time(NULL));
	dw_initialize_generator(0); 
	
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
		{"lHigh", required_argument, 0, 'L'}, 
		{"lLow", required_argument, 0, 'l'},
		{"H0", required_argument, 0, 'h'}, 
		{"T0", required_argument, 0, 't'}, 
		{"TK", required_argument, 0, 'T'},
		{"thin", required_argument, 0, 'i'}, 
		{"ndraws", required_argument, 0, 'w'},
		{"nInitial", required_argument, 0, 'I'}, 
		{"HillClimb", required_argument, 0, 'C'},
		{"TuningDone", no_argument, 0, 'D'},
		{0, 0, 0, 0}
	}; 
	int option_index = 0; 
	
	string data_file_name; 
	double minus_infinity = -1.0e30; 
	
	CEESParameter parameter; 
	parameter.storage_marker = 10000; 
	stringstream convert; 
	convert.str(string()); 
	convert << time(NULL); 
	parameter.run_id = convert.str(); 
	parameter.storage_dir = getenv("HOME")+string("/state_model/mssm/result/");
	bool if_original = false;
	parameter.number_energy_level = 10; 
 	parameter.pee =0.3; 
	parameter.highest_level = parameter.number_energy_level -1; 
	parameter.lowest_level = 0; 
	parameter.h0 = 0.0; 
	parameter.t0 = 1.0; 
	parameter.tk_1 = 1000.0; 
	parameter.deposit_frequency = 50; 
	parameter.simulation_length = 10000; 
	size_t n_initial = 1; 
	size_t number_hill_climb = 0; 
	bool if_tuning_done = false; 

	while (1)
	{
		int c = getopt_long(argc, argv, "d:m:r:on:H:L:l:h:t:c:i:w:I:C:D", long_options, &option_index); 
		if (c == -1)
			break; 
		switch(c)
		{
			case 'd':
				data_file_name = string(optarg); break; 
			case 'm':
				minus_infinity = atof(optarg); break; 
			case 'r':
				parameter.run_id = string(optarg); break; 
			case 'o':
				if_original = true; break;
			case 'n':
				parameter.number_energy_level = atoi(optarg); break;  	
			case 'L':
				parameter.highest_level = atoi(optarg); break; 
			case 'l':
				parameter.lowest_level = atoi(optarg); break; 
			case 'h':
				parameter.h0 = atof(optarg); break; 
			case 't':
				parameter.t0 = atof(optarg); break; 
			case 'T': 
				parameter.tk_1 = atof(optarg); break; 
			case 'i':
				parameter.deposit_frequency = atoi(optarg); break; 
	 		case 'w':
				parameter.simulation_length = atoi(optarg); break; 
			case 'I':
				n_initial = atoi(optarg); break; 
			case 'C':
				number_hill_climb = atoi(optarg); break; 
			case 'D':
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
	parameter.SetTemperature(); 
	CStorageHead storage(my_rank, parameter.run_id, parameter.storage_marker,parameter.storage_dir,parameter.number_energy_level); 

	if (my_rank == 0)
	{
		if(!storage.makedir())
		{
			cerr << "Error in making directory for " << parameter.run_id << endl;
			double *sMessage= new double [N_MESSAGE];
			for (int i=1; i<nNode; i++)
				MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD); 
			delete [] sMessage; 
			exit(1);  
		}
		master_deploying(nNode, if_tuning_done, number_hill_climb, n_initial, parameter, storage); 
	}
	else 
		slave_computing(if_original, number_hill_climb, n_initial, data_file_name, minus_infinity, parameter, storage); 
}
