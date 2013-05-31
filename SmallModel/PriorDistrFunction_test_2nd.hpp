#ifndef _PRIOR_DISTR_FUNCTION_TEST_2ND_
#define _PRIOR_DISTR_FUNCTION_TEST_2ND_

#include "CMSSM.hpp"

class PriorDistrFunction_test_2nd : public PriorDistributionFunction
{
public:
	virtual double log_pdf(const TDenseVector &x) ;
        PriorDistrFunction_test_2nd() : PriorDistributionFunction() {}
        PriorDistrFunction_test_2nd(const TDenseVector &p) : PriorDistributionFunction(p) {}
        virtual ~PriorDistrFunction_test_2nd() {}
};

#endif
