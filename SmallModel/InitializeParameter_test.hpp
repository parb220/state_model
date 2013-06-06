#ifndef _INITIALIZE_PARAMETER_
#define _INITIALIZE_PARAMETER_
#include "dw_dense_matrix.hpp"

TDenseVector InitializeParameter(size_t n, const TDenseVector &fixed_parameter); 
int RandomInit_test_1st(TDenseVector &x, const std::vector<int> &pos); 
int RandomInit_test_2nd(TDenseVector &x, const std::vector<int> &pos); 

#endif
