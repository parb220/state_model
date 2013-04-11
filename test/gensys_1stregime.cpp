/*
 * Not necessary, because the original matlab file does this for testing
 */

#include "dw_dense_matrix.hpp"
#include "findequilibrium.hpp"

using namespace std; 

bool gensys_1stregime(const TDenseVector &parameter)
{
	bool error_code = false; 
	size_t nEquation = 7; 	// number of equations
	double div = 1.0+1.0e09; 

	double gsigma = 2.0;
	double gpsi2 = 0.5;
	double b = 0.7;
	double ggamma = 0.3;
	double geta = 0.5;
	double gphipi = 2.5;
	double gphiy = 0.5;
	double grhomu = 0.9;
	double grhorn = 0.9;	// Natural rate of interest.
	double grhoi = 0.9;  	// Monetary policy.
	double gsigmamu = 0.1;	
	double gsigmarn = 0.1;	// Natural rate of interest.
	double gsigmai = 0.1;	// Monetary policy.
	double giota = 0.0;	// ZLB adjustment for level of fall in natural rate of interest.

	double glambdaz = parameter[GLAMBDAZ_STATE]; 
	double rn = parameter[RN_STATE]; 
	double gpistar = parameter[GPISTAR_STATE]; 
	double sc = parameter[SC_STATE]; 
	double Rlow = parameter[RLOW_STATE]; 

	double Rstar = gpistar + rn; 
	double gbeta = glambdaz / (1.0+rn); 
	double istar = -Rstar + Rlow; 
	double ci = 0.0; 
	double crn = 0.0; 
	double gchir = 1.0; 

	TDenseMatrix A00(nEquation,nEquation,0.0);
	A00.SetElement(1.0+gbeta*ggamma,0,0); 
	A00.SetElement(-gpsi2*(geta+gsigma/sc),0,1); 	
	A00.SetElement(-gpsi2,0,3); 
	A00.SetElement(-gbeta,0,5); 
	A00.SetElement(gsigma+(gsigma-1.0)*b,1,1); 
	A00.SetElement(sc,1,2); 
	A00.SetElement(-sc,1,4); 
	A00.SetElement(-sc,1,5); 
	A00.SetElement(-gsigma,1,6); 
	A00.SetElement(-(1.0-grhoi)*gphipi,2,0);
	A00.SetElement(-(1.0-grhoi)*gphiy,2,1);
	A00.SetElement(1.0,2,2); 
	A00.SetElement(-(1.0-grhoi)*gchir,2,4);
	A00.SetElement(1.0,3,4); 
	A00.SetElement(1.0,4,3); 
	A00.SetElement(1.0,5,0); 
	A00.SetElement(1.0,6,1);  

	TDenseMatrix B00(nEquation, nEquation,0.0); 
	B00.SetElement(ggamma,0,0);
	B00.SetElement(-gpsi2*(gsigma-1.0)*b/sc,0,1);
        B00.SetElement((gsigma-1.0)*b,1,1);
        B00.SetElement(grhoi,2,2);
        B00.SetElement(grhorn,3,4);
        B00.SetElement(grhomu,4,3);
        B00.SetElement(1.0,5,5);
        B00.SetElement(1.0,6,6);

	TDenseVector C00(nEquation,0.0); 
	C00.SetElement(ci,2); 
	C00.SetElement(crn,3); 

	TDenseMatrix Psi00(nEquation,3,0.0); 
	Psi00.SetElement(gsigmai,2,2);
        Psi00.SetElement(gsigmarn,3,0);
        Psi00.SetElement(gsigmamu,4,1);

	TDenseMatrix Pi00(nEquation,2,0.0);
        Pi00.SetElement(1.0,5,0);
        Pi00.SetElement(1.0,6,1);

	if (Rank(A00) < nEquation)
		error_code = true; 

	// Stops at Line 97 of gensys_1stregime
}
