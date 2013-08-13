#include <mpi.h>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "dw_dense_matrix.hpp"
#include "est_all_mpi.hpp"

using namespace std;

class greater_operator
{
public: 
	bool operator() (const TDenseVector &x, const TDenseVector &y)
	{
		if (x.dim ==0 || y.dim == 0 || x.dim != y.dim)
		{
			cerr << "solution_comp_operator(): dimension must be non-zero and consistent.\n"; 
			abort(); 
		}
		return (x[0] > y[0]); 
	}
} solution_greater;  

void est_all_master(size_t n_iterations, size_t n_seeds, const string &output_file_name ) 
{
	int n_nodes, sMessage, rMessage;
	MPI_Status status; 
	MPI_Comm_size(MPI_COMM_WORLD, &n_nodes); 	
	
	vector<TDenseVector> solutions, best_solutions, solutions_merged; 
	string file_name; 
	stringstream convert; 

	for (unsigned int i_iteration=0; i_iteration<n_iterations; i_iteration++)
	{
		best_solutions.clear(); 
		sMessage = i_iteration; 
		for (unsigned int i_node =1; i_node<(int)n_nodes; i_node ++)
			MPI_Send(&sMessage, 1, MPI_INT, i_node, WORK_TAG, MPI_COMM_WORLD); 
		
		for (unsigned int i_node =1; i_node<(int)n_nodes; i_node ++)
		{
			MPI_Recv(&rMessage, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);				if (status.MPI_TAG == WORK_TAG)
			{
				// load solutions 
				convert.str(string()); 
				convert << "." << i_iteration << "." << status.MPI_SOURCE; 
				file_name = output_file_name + convert.str(); 
				if (LoadSolution(file_name, solutions) == IO_SUCCESS)
				{
					// descending order
					sort(solutions.begin(), solutions.end(), solution_greater); 
					// merge solutions and best_solutions
					solutions_merged.resize(solutions.size()+best_solutions.size()); 
					merge(solutions.begin(), solutions.end(), best_solutions.begin(), best_solutions.end(), solutions_merged.begin(), solution_greater);
					// pick the top n_seeds merged_solutions
					best_solutions.resize(n_seeds < solutions_merged.size() ? n_seeds: solutions_merged.size()); 
					for (unsigned int i_seed=0; i_seed<n_seeds; i_seed++)
						best_solutions[i_seed].CopyContent(solutions_merged[i_seed]); 
				}
				else 
				{
					cerr << "Error in loading solutions from " << file_name << endl; 
					abort(); 
				}
			}
		}
		// save seeds
		convert.str(string()); 
		convert << "." << i_iteration+1 << ".ini"; 
		file_name = output_file_name + convert.str(); 
		if (SaveSeed(file_name, best_solutions) != IO_SUCCESS)
		{
			cerr << "Error in saving seeds to " << file_name << endl; 
			abort(); 
		}
	}

	// end 
	sMessage = 0; 
	for (unsigned int i_node =1; i_node<(int)n_nodes; i_node ++)
		MPI_Send(&sMessage, 1, MPI_INT, i_node, FINISH_TAG, MPI_COMM_WORLD);
}
