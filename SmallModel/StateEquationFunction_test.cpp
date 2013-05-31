#include <cmath>
#include "CMSSM_Error_Code.hpp"
#include "StateEquationFunction_test.hpp"

using namespace std; 

int StateEquationFunction_test::convert(vector<TDenseVector> &b, vector<TDenseMatrix> &F, vector<TDenseMatrix> &Phi_e, vector<TDenseMatrix> &V, const vector<vector<TDenseMatrix> > &A, const vector<vector<TDenseMatrix> > &B, const vector<vector<TDenseMatrix> > &Psi, const vector<vector<TDenseMatrix> >&Pi, const vector<vector<TDenseVector> > &C, const TDenseVector &x, size_t nZ, size_t nE, size_t nExpectation, size_t nNu)
// Return value of error has the following meaning
// 	SUCCESS: success
// 	1: some A returned by MakeABPsiPiC is not invertible
// 	2: p0<=0 or p0>1 or p1<=0 or p1>1
// 	3: number of explosive roots is not equal to the number of expectational error. Thus
// 		regime 0 is not determinate
// 	4: V0a*LeftSolve(A[0][0],Pi[0][0]) is not of full row rank
// 	5: V0a*LeftSolve(A[0][1],Pi[0][1]) is not of full row rank
// 	6: V1a*LeftSolve(A[1][1],Pi[1][1]) is not of full row rank
// 	7: V1a*LeftSolve(A[1][0],Pi[1][0]) is not of full row rank
// 	8: LeftSolve(A[i][i]-B[i][i],C[i][i]) error for some i 
// 	9: Annihilator error
{
	// sizes 
	size_t originalNS = 2; 
	size_t nR = 6; 

	// probability of staying in zero lower bound
	double p0 = x(x.dim-2);
	double p1 = x(x.dim-1);
	if ( (p0<=0) || (p0>1) || (p1<=0) || (p1>1) )
		return 2; 

	//intermediate_b, intermediate_F, intermediate_Phi_e, and intermediate_V
	vector<TDenseVector> intermediate_b(nR); 
	vector<TDenseMatrix> intermediate_F(nR);
	vector<TDenseMatrix> intermediate_Phi_e(nR); 
	vector<TDenseMatrix> intermediate_V(nR);
	for (unsigned int i=0; i<nR; i++)
	{
		intermediate_b[i] = TDenseVector(nZ,0.0); 
		intermediate_F[i] = TDenseMatrix(nZ,nZ,0.0); 
		intermediate_Phi_e[i] = TDenseMatrix(nZ,nE,0.0); 
		intermediate_V[i] = TDenseMatrix(nZ,nZ,0.0); 
	}

	// c_hat, d, and I
	vector<TDenseVector>c_hat(originalNS); 
	vector<TDenseVector>d(originalNS);
	for (unsigned int i=0; i<originalNS; i++)
	{
		c_hat[i] = TDenseVector(nZ,0.0); 
		d[i] = TDenseVector(nZ,0.0);
	}

	try {	c_hat[0] = LeftSolve( (A[0][0]-B[0][0]), C[0][0] ); }
	catch(...)
	{	return 8; }
	try {	c_hat[1] = LeftSolve( (A[1][1]-B[1][1]), C[1][1] ); }
	catch (...)
	{	return 8; }

	d[0] = C[0][1] - A[0][1]*c_hat[0] + B[0][1]*c_hat[1]; 
	d[1] = C[1][0] - A[1][0]*c_hat[1] + B[1][0]*c_hat[0];
	TDenseMatrix I; I.Identity(nZ); 

	// Multiply by inverse of A
	TDenseMatrix A00InvPi00 = LeftSolve(A[0][0], Pi[0][0]); 
	TDenseMatrix A10InvPi10 = LeftSolve(A[1][0], Pi[1][0]);	// ??? Pi[1][0] or Pi[0][0] 
	TDenseMatrix A11InvPi11 = LeftSolve(A[1][1], Pi[1][1]); 
	TDenseMatrix A01InvPi01 = LeftSolve(A[0][1], Pi[0][1]); 

	TDenseMatrix A00InvB00 = LeftSolve(A[0][0], B[0][0]); 
	TDenseMatrix A10InvB10 = LeftSolve(A[1][0], B[1][0]); 
	TDenseMatrix A11InvB11 = LeftSolve(A[1][1], B[1][1]); 
	TDenseMatrix A01InvB01 = LeftSolve(A[0][1], B[0][1]); 

	TDenseMatrix A00InvPsi00 = LeftSolve(A[0][0], Psi[0][0]); 
	TDenseMatrix A01InvPsi01 = LeftSolve(A[0][1], Psi[0][1]); 
	TDenseMatrix A11InvPsi11 = LeftSolve(A[1][1], Psi[1][1]); 
	TDenseMatrix A10InvPsi10 = LeftSolve(A[1][0], Psi[1][0]);

	// Annihilators
	TDenseMatrix V0a, V1a; 
	TDenseVector E0a, E1a;
	try { Annihilator(V0a, E0a, A00InvB00); }
	catch(...)
	{	return 9; }
	try { Annihilator(V1a, E1a, A11InvB11); }
	catch(...)
	{	return 9; }

	// Regime 0 determinate?
	if (V0a.rows != nExpectation)
		return 3; 
	
	// Set Mij, Nij, Delta22, gamma12
	TDenseMatrix X = V0a * A00InvPi00; 
	if (Rank(X) != nExpectation)
		return 4; 
	TDenseMatrix M00 = LeftSolve(X, V0a); 

	X = V0a*A01InvPi01;
	if (Rank(X) != nExpectation)
		return 5; 
	TDenseMatrix M01 = LeftSolve(X, V0a); 	

	TDenseMatrix N11, M11, N10, M10; 
	if (V1a.rows == 0)	
	{
		N11.Identity(nExpectation); 
		M11.Zeros(nExpectation, nZ); 
		N10.Identity(nExpectation); 
		M10.Zeros(nExpectation, nZ); 
	}
	else 
	{
		X = V1a*A11InvPi11; 
		if (Rank(X) != V1a.rows)
			return 6;		
		NullSpace(N11, X); 	// N11 = NullSpace(X); 
		M11.Multiply(X.GeneralizedInverse(), V1a);	// M11 = pinv(X)*V1a
		
		X = V1a * A10InvPi10; 
		if (Rank(X) != V1a.rows)
			return 7;
		NullSpace(N10, X); 	// N10 = NullSpace(X);
		M10.Multiply(X.GeneralizedInverse(), V1a); 	// M10 = pinv(X)*V1a
	
		if (N11.cols != N10.cols)
			return 7;
	}

	// Delta11 = reshape(x(end-1-N11.cols*nE:end-2),N11.cols, nE); 
	TDenseMatrix Delta11(N11.cols, nE); 
	unsigned int counter_reshape = x.dim-2-Delta11.rows *Delta11.cols; 
	for (unsigned int j=0; j<Delta11.cols; j++)
	{
		for (unsigned int i=0; i<Delta11.rows; i++)
		{
			Delta11.SetElement(x(counter_reshape), i,j);
			counter_reshape++;
		}
	}	
	// gamma10 = x(end-1-nExpectations*nE-N10.cols:end-2-nExpectations*nE); 
	TDenseVector gamma10(N10.cols); 
	counter_reshape = x.dim-2-nExpectation*nE-N10.cols; 
	for (unsigned int i=0; i<gamma10.dim; i++)
	{
		gamma10.SetElement(x(counter_reshape), i); 
		counter_reshape ++; 
	}
	// G0 = -M10*A10InvPsi10 + N10*Delta11
	// G1 = -M01*A01InvPsi01
	TDenseMatrix G0 = -M10*A10InvPsi10 + N10*Delta11; 
	TDenseMatrix G1 = -M01*A01InvPsi01; 

	// ======= s(t+1) =0, s(t) = 0, p(t) =1 =========================
	X = I - A00InvPi00*M00; 
	intermediate_F[0] = X * A00InvB00; 
	intermediate_Phi_e[0] = X * A00InvPsi00; 
	intermediate_b[0] = c_hat[0]-intermediate_F[0]*c_hat[0];  

	// ======= s(t+1)=0 s(t)=0 p(t) =p0 =============================
	X=A00InvPi00*M10*((1.0-p0)/p0);
	intermediate_F[1] = A00InvB00 + X*A10InvB10; 
	intermediate_Phi_e[1] = A00InvPsi00 + A00InvPi00*G0; 
	intermediate_b[1] = c_hat[0] - intermediate_F[1]*c_hat[0] + X*LeftSolve(A[1][0], d[1]) +  A00InvPi00*N10*gamma10*((1.0-p0)/p0) ; 

	// ====== s(t+1)=1 s(t)=0 p(t)=p0 ===============================
	X = I - A10InvPi10*M10; 
	intermediate_F[2] = X * A10InvB10; 
	intermediate_Phi_e[2] = X*A10InvPsi10 + A10InvPi10*N10*Delta11;
	intermediate_b[2] = c_hat[1] - intermediate_F[2]*c_hat[0] + X*( LeftSolve(A[1][0], d[1]) - A10InvPi10*N10*gamma10 );  

	// ======= s(t+1)=1 s(t)=1 p(t)=1 ===============================
	X = I - A11InvPi11*M11; 
	intermediate_F[3] = X*A11InvB11; 
	intermediate_Phi_e[3] = X*A11InvPsi11 + A11InvPi11*N11*Delta11; 
	intermediate_b[3] = c_hat[1] - intermediate_F[3]*c_hat[1]; 

	// ======= s(t+1)=1 s(t)=1 p(t)=p1 ==============================
	X = A11InvPi11*M01 * ((1.0-p1)/p1); 
	intermediate_F[4] = A11InvB11 + X*A01InvB01; 
	intermediate_Phi_e[4] = A11InvPsi11 + A11InvPi11*G1; 
	intermediate_b[4] = c_hat[1] - intermediate_F[4]*c_hat[1] + X*LeftSolve(A[0][1],d[0]); 
	
	// ======= s(t+1)=0 s(t)=1 p(t)=p1 ==============================
	X = I - A01InvPi01*M01; 
	intermediate_F[5] = X*A01InvB01; 
	intermediate_Phi_e[5] = X* A01InvPsi01; 
	intermediate_b[5] = c_hat[0]-intermediate_F[5]*c_hat[1] + X *LeftSolve(A[0][1], d[0]); 

	// Compute V
	if (nE)
	{
		for (unsigned int i=0; i<nR; i++)
			intermediate_V[i].MultiplyTranspose(intermediate_Phi_e[i], intermediate_Phi_e[i]);  
	}

	// Form solution 
	b = vector<TDenseVector>(nNu); 
	F = vector<TDenseMatrix>(nNu); 
	Phi_e = vector<TDenseMatrix>(nNu); 
	V = vector<TDenseMatrix>(nNu); 
	for (unsigned int i=0; i<nNu; i++)
	{
		b[i] = TDenseVector(nZ,0.0); 
		F[i] = TDenseMatrix(nZ,nZ,0.0); 
		Phi_e[i] = TDenseMatrix(nZ,nE,0.0); 
		V[i] = TDenseMatrix(nZ,nZ,0.0); 
	}

	// 0 <==> 0
	b[0] = intermediate_b[0]; 
	F[0] = intermediate_F[0]; 
	Phi_e[0] = intermediate_Phi_e[0]; 
	V[0] = intermediate_V[0]; 

	// 1 <==> 0
	b[1] = intermediate_b[0]; 
	F[1] = intermediate_F[0]; 
	Phi_e[1] = intermediate_Phi_e[0]; 
	V[1] = intermediate_V[0]; 

	// 5 <==> 1
	b[5] = intermediate_b[1]; 
	F[5] = intermediate_F[1]; 
	Phi_e[5] = intermediate_Phi_e[1];
	V[5] = intermediate_V[1];

	// 6 <==> 2
	b[6] = intermediate_b[2]; 
	F[6] = intermediate_F[2]; 
	Phi_e[6] = intermediate_Phi_e[2]; 
	V[6] = intermediate_V[2]; 

	// 10 <==> 3
	b[10] = intermediate_b[3]; 
	F[10] = intermediate_F[3]; 
	Phi_e[10] = intermediate_Phi_e[3]; 
	V[10] = intermediate_V[3]; 

	// 11 <==> 3
	b[11] = intermediate_b[3]; 
	F[11] = intermediate_F[3]; 
	Phi_e[11] = intermediate_Phi_e[3]; 
	V[11] = intermediate_V[3]; 

	// 12 <==> 5
	b[12] = intermediate_b[5]; 
	F[12] = intermediate_F[5]; 
	Phi_e[12] = intermediate_Phi_e[5]; 
	V[12] = intermediate_V[5]; 

	// 15 <==> 4
	b[15] = intermediate_b[4]; 
	F[15] = intermediate_F[4]; 
	Phi_e[15] = intermediate_Phi_e[4]; 
	V[15] = intermediate_V[4]; 

	return SUCCESS; 
}
