#include "MakeLbUb_ststm1.hpp"

void MakeLbUb_ststm1(TDenseVector &lb, TDenseVector &ub, int choice)
{
	const double INFINITE_BOUND = 10E20; 
	size_t nX = 14;
	size_t nRegime = 2;  

	double gsigma = 2.0;
	double gpsi2 = 0.5;
	double b = 0.7;
	double ggamma = 0.3;
	double geta = 0.5;
	double gphipi = 2.5;
	double gphiy = 0.5;
	double grhomu = 0.9;
	double grhorn = 0.9;	// Natural rate of interest.
	double grhor = 0.9;	// Monetary policy.
	double gsigmamu = 0.1;	
	double gsigmarn = 0.1;	// Natural rate of interest.
	double gsigmai = 0.1;  	// Monetary policy.
	double giota = 0.0;  	// ZLB adjustment for level of fall in natural rate of interest.	

	lb.Zeros(nX*nRegime+nRegime); 	// nX * nRegime: nX parameters per regime
	ub.Zeros(nX*nRegime+nRegime);	// nRegime:	probability per regime
	unsigned int i=0, j=nX; 
	if (choice == 1)	// Wide search	
	{
		lb.SetElement(1.0,i); ub.SetElement(5.0,i); i++;	// gsigma (risk aversion)
		ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;   // gpsi2 (slope of Phillips curve)
		ub.SetElement(0.9999,i);   ub.SetElement(0.9999,j); i++; j++;  // b (habit)
		ub.SetElement(0.9999,i);   ub.SetElement(0.9999,j); i++; j++;   // ggamma (inflation indexation)
		ub.SetElement(3.0,i); i++;					// geta (inverse labor elasticity)
		lb.SetElement(1.0,i); ub.SetElement(5.0,i); i++;		// gphipi (policy response to inflation)
		ub.SetElement(1.0,i); i++;					// gphiy (plicy response to output gap)
		ub.SetElement(0.9999,i); ub.SetElement(0.9999,j); i++; j++;  	// grhomu (persistence to markup)
		ub.SetElement(0.9999,i); ub.SetElement(0.9999,j); i++; j++;  	// grhorn (persistence to natural rate of interest)
		ub.SetElement(0.9999,i); ub.SetElement(0.9999,j); i++; j++; 	// grhoi (persistence to monetary policy)
		ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;	// gsigmamu (s.d. of markup shock)
		ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;     // gsigmarn (s.d. of demand shock)
		ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;   	// gsigmai (s.d. of policy shock)
		lb.SetElement(-1.0,i); ub.SetElement(1.0,i); 			// giota (level adjustment for ZLB drop in the interest rate)
		lb.SetElement(-INFINITE_BOUND,23); ub.SetElement(INFINITE_BOUND,23);    // Sunspot component
		lb.SetElement(-INFINITE_BOUND,24); ub.SetElement(INFINITE_BOUND,24); 	// Sunspot component
		lb.SetElement(-INFINITE_BOUND,25); ub.SetElement(INFINITE_BOUND,25);    // Sunspot component
		lb.SetElement(-INFINITE_BOUND,26); ub.SetElement(INFINITE_BOUND,26);    // Sunspot component
		lb.SetElement(-INFINITE_BOUND,27); ub.SetElement(INFINITE_BOUND,27);    // Sunspot component
		lb.SetElement(-INFINITE_BOUND,28); ub.SetElement(INFINITE_BOUND,28);    // Sunspot component
		ub.SetElement(1.0, 29); 					// probability of staying in regime 2
	}
	else if (choice == 2)	// Unrestricted search on gpsi2
	{
		lb.SetElement(1.0,i); ub.SetElement(5.0,i); i++;        // gsigma (risk aversion)
		lb.SetElement(-0.1,j); ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;	// gpsi2 (slope of Phillips curve)
		ub.SetElement(0.9999,i);   ub.SetElement(0.9999,j); i++; j++;  // b (habit)
		ub.SetElement(0.9999,i);   ub.SetElement(0.9999,j); i++; j++;   // ggamma (inflation indexation)
		ub.SetElement(2,i); i++; 					// geta (inverse labor elasticity)
		lb.SetElement(1.0,i); ub.SetElement(5.0,i); i++;                // gphipi (policy response to inflation)
		ub.SetElement(1.0,i); i++;                                      // gphiy (plicy response to output gap)
		ub.SetElement(0.9999,i); ub.SetElement(0.9999,j); i++; j++;     // grhomu (persistence to markup)
                ub.SetElement(0.9999,i); ub.SetElement(0.9999,j); i++; j++;     // grhorn (persistence to natural rate of interest)
                ub.SetElement(0.9999,i); ub.SetElement(0.9999,j); i++; j++;     // grhoi (persistence to monetary policy)
		ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;     // gsigmamu (s.d. of markup shock)
                ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;     // gsigmarn (s.d. of demand shock)
                ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;     // gsigmai (s.d. of policy shock)
		lb.SetElement(-INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,i); 		// giota (level adjustment for ZLB drop in the interest rate)
		lb.SetElement(-INFINITE_BOUND,23); ub.SetElement(INFINITE_BOUND,23);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,24); ub.SetElement(INFINITE_BOUND,24);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,25); ub.SetElement(INFINITE_BOUND,25);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,26); ub.SetElement(INFINITE_BOUND,26);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,27); ub.SetElement(INFINITE_BOUND,27);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,28); ub.SetElement(INFINITE_BOUND,28);    // Sunspot component
		ub.SetElement(1.0, 29);                                         // probability of staying in regime 2
	}
	else 	// Unrestricted search
	{
		lb.SetElement(1.0,i); ub.SetElement(5.0,i); i++;        // gsigma (risk aversion)
		lb.SetElement(-INFINITE_BOUND,j); ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;      // gpsi2 (slope of Phillips curve)
		ub.SetElement(0.9999,i);   ub.SetElement(0.9999,j); i++; j++;  // b (habit)
                ub.SetElement(0.9999,i);   ub.SetElement(0.9999,j); i++; j++;   // ggamma (inflation indexation)
		ub.SetElement(INFINITE_BOUND,i); i++;                           // geta (inverse labor elasticity)
		lb.SetElement(1.0,i); ub.SetElement(INFINITE_BOUND,i); i++; 	// gphipi (policy response to inflation)
		ub.SetElement(INFINITE_BOUND,i); i++; 				// gphiy (plicy response to output gap)
		ub.SetElement(0.9999,i); ub.SetElement(0.9999,j); i++; j++;     // grhomu (persistence to markup)
                ub.SetElement(0.9999,i); ub.SetElement(0.9999,j); i++; j++;     // grhorn (persistence to natural rate of interest)
                ub.SetElement(0.9999,i); ub.SetElement(0.9999,j); i++; j++;     // grhoi (persistence to monetary policy)
		ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;     // gsigmamu (s.d. of markup shock)
                ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;     // gsigmarn (s.d. of demand shock)
                ub.SetElement(INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,j); i++; j++;     // gsigmai (s.d. of policy shock)
		lb.SetElement(-INFINITE_BOUND,i); ub.SetElement(INFINITE_BOUND,i);              // giota (level adjustment for ZLB drop in the interest rate)
		lb.SetElement(-INFINITE_BOUND,23); ub.SetElement(INFINITE_BOUND,23);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,24); ub.SetElement(INFINITE_BOUND,24);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,25); ub.SetElement(INFINITE_BOUND,25);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,26); ub.SetElement(INFINITE_BOUND,26);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,27); ub.SetElement(INFINITE_BOUND,27);    // Sunspot component
                lb.SetElement(-INFINITE_BOUND,28); ub.SetElement(INFINITE_BOUND,28);    // Sunspot component
                ub.SetElement(1.0, 29);	
	}	
}
