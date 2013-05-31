#ifndef _READ_WRITE_FILE_
#define _READ_WRITE_FILE_ 

#include <vector>
#include "dw_dense_matrix.hpp"

using namespace std; 
int LoadData(TDenseVector&, vector<TDenseVector>& , const string &filename);
int LoadInitialValue(vector<TDenseVector> &, const string &filename); 

#endif
