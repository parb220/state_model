#include <fstream>
#include "CMSSM_Error_Code.hpp"
#include "ReadWriteFile.hpp"

using namespace std; 
int LoadInitialValue(vector<TDenseVector> &initialX, const string &file_name)
{
	ifstream input_file(file_name.c_str()); 
	if (!input_file)
		return ERROR_OCCURRED; 
	size_t nData, nDim; 
	input_file >> nData >> nDim; 
	if (nData <= 0 || nDim <=1 )
	{
		input_file.close(); 
		return ERROR_OCCURRED; 
	}
	initialX.resize(nData); 
	double field; 
	for (unsigned int i=0; i<nData; i++)
	{
		initialX[i].Resize(nDim); 
		for (unsigned int j=0; j<nDim; j++)
		{
			input_file >> field; 
			initialX[i].SetElement(field, j); 
		}
	}
	input_file.close(); 
	return SUCCESS; 
}

int LoadData(TDenseVector &qm_date, vector<TDenseVector> &qdata, const string &file_name)
{
	ifstream input_file(file_name.c_str()); 
	if (!input_file)
		return ERROR_OCCURRED;	// error occurred
	size_t nSample, nDim; 	// number of records, and dimension of each record
	input_file >> nSample >> nDim; 
	if (nSample <=0 || nDim <= 1)	// should at least contain 1 record
	{
		input_file.close(); 
		return ERROR_OCCURRED; 	// and each record should as least has 2 fields
					// (date and data)
	}
	qm_date.Resize(nSample); 
	qdata.resize(nSample); 
	double field; 
	for (unsigned int i=0; i<nSample; i++)
	{
		input_file >> field; 
		qm_date.SetElement(field, i);
		qdata[i].Resize(nDim-1); 
		for (unsigned int j=0; j<nDim-1; j++)
		{
			input_file >> field; 
			qdata[i].SetElement(field, j); 
		}
	}
	input_file.close(); 
	return SUCCESS; 
} 
