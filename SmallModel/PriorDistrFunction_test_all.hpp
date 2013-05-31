#ifndef _PRIOR_DISTR_FUNCTION_TEST_ALL_
#define _PRIOR_DISTR_FUNCTION_TEST_ALL_

#include "CMSSM.hpp"

class PriorDistrFunction_test_all : public PriorDistributionFunction
{
public:
	virtual double log_pdf(const TDenseVector &x) ;
        PriorDistrFunction_test_all() : PriorDistributionFunction() {}
        PriorDistrFunction_test_all(const TDenseVector &p) : PriorDistributionFunction(p) {}
        virtual ~PriorDistrFunction_test_all() {}
};

#endif
