#include <fstream>
#include "ReadWriteFile.hpp"

using namespace std; 

bool LoadData(vector<double> &qm_date, vector<TDenseVector> &qdata, const string &file_name)
{
	ifstream input_file(file_name.c_str()); 
	if (!input_file)
		return true;	// error occurred
	size_t nSample, nDim; 	// number of records, and dimension of each record
	input_file >> nSample >> nDim; 
	if (nSample <=0 || nDim <= 1)	// should at least contain 1 record
		return true; 		// and each record should as least has 2 fields
					// (date and data)
	qm_date.resize(nSample); 
	qdata = vector<TDenseVector>(nSample,TDenseVector(nDim-1,0.0)); 	
	double field; 
	for (unsigned int i=0; i<nSample; i++)
	{
		input_file >> qm_date[i]; 
		for (unsigned int j=0; j<nDim-1; j++)
		{
			input_file >> field; 
			qdata[i].SetElement(field, j); 
		}
	}
	input_file.close(); 
	return false; 
} 
