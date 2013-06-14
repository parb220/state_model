#ifndef _EST_ALL_MPI_
#define _EST_ALL_MPI_

#include <vector>
#include <string>
#include "dw_dense_matrix.hpp"

using namespace std; 

const int WORK_TAG =1; 
const int FINISH_TAG = 0; 
const int IO_ERROR = 1; 
const int IO_SUCCESS = 0; 

void est_all_master(size_t, size_t, const string &);
void est_all_slave(const string & data_file_name, const string & initial_value_file_name, size_t n_tries, double minus_infinity, const string &output_file_name);
int LoadSolution(const string &file_name, vector<TDenseVector> &solutions);
int SaveSeed(const string &file_name, const vector<TDenseVector> &best_solutions);
int SaveSolution(const string &file_name, const std::vector<TDenseVector>& solutions);

#endif
