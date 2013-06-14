#include <iomanip>
#include <cstdio>
#include <ctime>
#include <mpi.h>
#include <getopt.h>
#include "CMSSM_Error_Code.hpp"
#include "ReadWriteFile.hpp"
#include "dw_rand.h"
#include "dw_dense_matrix.hpp"
#include "est_all_mpi.hpp"

using namespace std;

int main(int argc, char **argv)
{
	static struct option long_options[] =
        {
		{"number_iterations", required_argument, 0, 'N'},
		{"number_seeds", required_argument, 0, 'S'},
                {"data_file", required_argument, 0, 'd'},
                {"initial_value_file", required_argument, 0, 'v'},
                {"max_tries", required_argument, 0, 't'},
                {"minus_infinity", required_argument, 0, 'i'},
                {"output_file", required_argument, 0, 'f'},
                {0, 0, 0, 0}
        };
	int option_index = 0; 
	string data_file_name, initial_value_file_name, output_file_name; 
	size_t n_iterations = 10, n_seeds = 1, n_tries = 30; 
	double minus_infinity = -1.0e30; 

	while (1)
	{
		int c=getopt_long(argc, argv, "N:S:d:v:t:i:f:", long_options, &option_index); 
		if (c == -1)
			break; 
		switch(c)
		{
			case 'N':
				n_iterations = atoi(optarg); break; 
			case 'S':
				n_seeds = atoi(optarg); break; 
			case 'd':
				data_file_name = string(optarg); break; 
			case 'v':
				initial_value_file_name = string(optarg); break; 
			case 't':
				n_tries = atoi(optarg); break; 
			case 'i': 
				minus_infinity = atof(optarg); break; 
			case 'f': 
				output_file_name = string(optarg); break; 
		}
	}

	if (data_file_name.length()<=0 || output_file_name.length() <=0)
	{
		cerr << "Usage: " << argv[0] << " -N iteration number -S seed number -d data file -v initial value file -t maximum tries -i minus infinity -f output file.\n";
                abort();
	}

	// test if data file exists
	TDenseVector qm_date;           // dates
        vector<TDenseVector> qdata;     // variables
        if (LoadData(qm_date, qdata, data_file_name) != SUCCESS)
        {
                cerr << "------LoadData(): error occurred ------\n";
                abort();
        }

	// random number generator
	dw_initialize_generator(time(NULL));

	// Initialize MPI 
	MPI_Init(&argc, &argv);
	int my_rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (my_rank == 0)
		est_all_master(n_iterations, n_seeds, output_file_name); 
	else 
		est_all_slave(data_file_name, initial_value_file_name, n_tries, minus_infinity, output_file_name);  	

	// end
	MPI_Finalize; 
} 
