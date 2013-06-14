#include "est_all_mpi.hpp"
#include <iomanip>
#include <vector>
#include <fstream>

using namespace std; 

int SaveSolution(const string &file_name, const std::vector<TDenseVector>& solutions)
{
	if (solutions.empty() || solutions[0].dim<=0)
		return IO_ERROR; 
	
	ofstream output_file(file_name.c_str()); 
	if (!output_file)
		return IO_ERROR; 
	output_file << solutions.size() << "\t" << solutions[0].dim << endl; 
	for (unsigned int i =0; i<solutions.size(); i++)
	{
		for (unsigned int j=0; j<solutions[i].dim; j++)
			output_file << setprecision(20) << solutions[i][j] << " "; 
		output_file << endl; 
	}
	output_file.close(); 
	return IO_SUCCESS; 
}

int LoadSolution(const string &file_name, vector<TDenseVector> &solutions)
{
	ifstream input_file(file_name.c_str()); 
	if (!input_file)
		return IO_ERROR; 
	size_t number, dim; 
	input_file >> number >> dim; 
	if (number <=0 || dim <=0)
		return IO_ERROR;
	solutions.resize(number); 
	for (unsigned int i=0; i<number; i++)
	{
		solutions[i].Resize(dim); 
		for (unsigned int j=0; j<dim; j++)
			input_file >> solutions[i][j]; 
	}
	input_file.close(); 
	return IO_SUCCESS; 
}

int SaveSeed(const string &file_name, const vector<TDenseVector> &best_solutions)
{
	if (best_solutions.empty() || best_solutions[0].dim<=0)
		return IO_ERROR; 
	ofstream output_file(file_name.c_str()); 
	if (!output_file)
		return IO_ERROR; 
	output_file << best_solutions.size() << "\t" << best_solutions[0].dim << endl; 	
	for (unsigned int i=0; i<best_solutions.size(); i++)
	{
		for (unsigned int j=2; j<best_solutions[i].dim; j++)
			output_file <<setprecision(20) << best_solutions[i][j] << " "; 
		output_file << endl; 
	}
	output_file.close(); 
	return IO_SUCCESS; 
}
