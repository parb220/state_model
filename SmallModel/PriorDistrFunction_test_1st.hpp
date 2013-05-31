#ifndef _PRIOR_DISTR_FUNCTION_TEST_1ST_
#define _PRIOR_DISTR_FUNCTION_TEST_1ST_

#include "CMSSM.hpp"

class PriorDistrFunction_test_1st : public PriorDistributionFunction
{
public:
	virtual double log_pdf(const TDenseVector &x) ;
        PriorDistrFunction_test_1st() : PriorDistributionFunction() {}
        PriorDistrFunction_test_1st(const TDenseVector &p) : PriorDistributionFunction(p) {}
        virtual ~PriorDistrFunction_test_1st() {}
};

#endif
