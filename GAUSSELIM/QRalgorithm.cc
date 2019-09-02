#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "irram-random.h"
#include "strassen.h"
#include <vector>
#include <iostream>
#include <ctime>
using namespace iRRAM;

INTEGER pt (int p)
{
    if (p == 1)
	return INTEGER(2);

    if (p % 2 == 0)
    {
	return pt (p/2)* pt (p/2);
    }
    else
	return pt(p/2) * pt(p/2) * 2;
}

RATIONAL err (int p)
{
    return RATIONAL(1, pt(p));
}

void print(REALMATRIX M)
{
    cout << "print matrix:\n";
    for (int i = 0; i < (int) M.maxrow; i++)
    {
	for (int j = 0; j< (int) M.maxcolumn; j++)
	{
	    cout << M(i,j) << " ";
	}
	cout << "\n";
    }
}

/*
 * Various matrix norms where current max-norm does not have
 * Sub-multiplicativity.
 * Fnorm : Frobenius norm
 * Cnorm : Column norm
 * Rnorm : Row norm
 */
REAL Fnorm(REALMATRIX M)
{
    REAL norm = 0;
    for (unsigned int i = 0; i < M.maxrow; i++)
	for (unsigned int j = 0; j < M.maxcolumn; j++)
	    norm += M(i,j)*M(i,j);
    return sqrt(norm);
}

REAL Cnorm(REALMATRIX M)
{
    return 0;
}

REAL Rnorm(REALMATRIX M)
{
    return 0;
}



/*
 * Some matrix operation utilities
 * Some of these should be in 
 * REALMATRIXs definiton
 */
REALMATRIX neg(REALMATRIX M)

{
    return zeroes(M.maxrow, M.maxcolumn) - M;
}

REALMATRIX transpose(REALMATRIX M)
{
    REALMATRIX N = REALMATRIX(M.maxcolumn, M.maxrow);
    for (unsigned int i = 0; i < M.maxcolumn; i++)
	for (unsigned int j = 0; j < M.maxrow; j++)
	    N(i,j) = M(j,i);
    return N;
}

// return colum matrix below some index (i,j)
REALMATRIX colum(REALMATRIX M, int i, int j)
{
    REALMATRIX A = REALMATRIX(M.maxrow - i, 1);
    for (int k = i; k < (int) M.maxrow; k++)
	A(k-i,0) = M(k, j);
    return A;
}

// [A B]
// [C D]
REALMATRIX block_matrix(REALMATRIX A, REALMATRIX B, REALMATRIX C, REALMATRIX D)
{
    int row = (int) (A.maxrow + C.maxrow);
    int col = (int) (A.maxcolumn + B.maxcolumn);
    REALMATRIX S = REALMATRIX(row, col);
    for(int i = 0; i < row; i ++)
    {
	for(int j = 0; j < col; j++)
	{
	    if (i < (int) A.maxrow && j < (int) A.maxcolumn)
		S(i,j) = A(i,j);
	    if (i < (int) A.maxrow && j >= (int) A.maxcolumn)
		S(i,j) = B(i,j - A.maxcolumn);
	    if (i >= (int) A.maxrow && j < (int) A.maxcolumn)
		S(i,j) = C(i - A.maxrow,j);
	    if (i >= (int) A.maxrow && j >= (int) A.maxcolumn)
		S(i,j) = D(i - A.maxrow, j - A.maxcolumn);
	}
    }
    return S;
}

// [A 0]
// [0 B]
REALMATRIX block_diag_matrix(REALMATRIX A, REALMATRIX B)
{
    int n = (int) (A.maxrow + B.maxrow);
    REALMATRIX S = REALMATRIX(n, n);
    for(int i = 0; i < n; i ++)
    {
	for(int j = 0; j < n; j++)
	{
	    if (i < (int) A.maxrow && j < (int) A.maxcolumn)
		S(i,j) = A(i,j);
	    else 
	    {   
		if (i >= (int) A.maxrow && j >= (int) A.maxcolumn)
		    S(i,j) = B(i - A.maxrow, j - A.maxcolumn);
		else
		    S(i,j) = 0;
	    }
	}
    }
    return S;
}

// [A B C]
// [E F G]
// [H I J] => F (M(ii) ~ M(j-1,j-1))
REALMATRIX block_decompose(REALMATRIX A, unsigned int i, unsigned int j)
{
    REALMATRIX M = REALMATRIX(j-i, j-i);
    for (unsigned int n = i; n < j; n++)
	for (unsigned int m = i; m < j; m++)
	    M(n-i,m-i) = A(n,m);
    return M;
}

void block_decompose(REALMATRIX A, unsigned int i, unsigned int j, REALMATRIX* L)
{
    int r = A.maxrow;
    int c = A.maxcolumn;
    L[0] = REALMATRIX(i,j);
    L[1] = REALMATRIX(i,c-j);
    L[2] = REALMATRIX(r-i,j);
    L[3] = REALMATRIX(r-i,c-j);
    for(int n=0; n<r; n++)
	for(int m=0; m<c; m++)
	    if(n<i)
		if(m<j)
		    L[0](n,m)=A(n,m);
		else
		    L[1](n,m-j)=A(n,m);
	    else
		if(m<j)
		    L[2](n-i,m)=A(n,m);
		else
		    L[3](n-i,m-j)=A(n,m);
}




/* Reducing a symmetric matrix into a Tridiagonal matrix.
 * First one requires the matrix to be special. 
 * Second one works for every symmetric matrices however,
 * it returns a tridiagonal matrix which approximates
 * eigenvalues of the original one. it is not a similarity transformation.
 * The natural question arises: if a hessenberg reduction computable...
 */
REALMATRIX hessenberg_reduction(REALMATRIX M)
{
    REALMATRIX c, u, P;
    REAL s;

    int n = (int) M.maxrow;
    for(int k = 0; k < n-2; k++)
    {

	c = colum(M, k+1, k);
	s = Fnorm(c);
	c(0, 0) = c(0, 0) - s;
	s = Fnorm(c);
	u = c / s;
	P = eye(n - k - 1) - 2 * u * transpose(u);
	P = block_diag_matrix(eye(k+1), P);
	M = mult(P, mult(M, P));
    }
    return M;
}

// output: Hessenberg matrix H whose eigenvalues are perturbed
// at most by 2^p (p is considered as a negative number)
REALMATRIX hessenberg_reduction(REALMATRIX M, int p)
{
    REALMATRIX c, u, P;
    REALMATRIX L[4];
    REAL s, q;
    REAL approx = power(2, p - 1);
    int n = (int) M.maxrow;
    REAL err;
    for(int k = 0; k < n-2; k++)
    {
	std::cout<<k;
	c = colum(M, k+1, k);
	s = Fnorm(c);
	//    q = sqrt(s*s - c(0,0)*c(0,0));
	err = 0;
	err = abs(c(1,0));
	for (int i=1; i< c.maxrow; i++)
	    err = maximum(err, abs(c(i,0)));

	if (choose(err < approx / n / n, err > approx / 2 / n / n ) == 2)
	{
	    c(0, 0) = c(0, 0) - s;
	    s = Fnorm(c);
	    u = c / s;
	    P = eye(n - k - 1) - 2 * u * transpose(u);
	    //    P = block_diag_matrix(eye(k+1), P);
	    block_decompose (M, k+1, k+1, L);
	    L[1] = L[1]*P;
	    L[2] = P*L[2];
	    L[3] = mult(P, mult(L[3], P));
	    for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
		    if(i < k+1)
			if(j > k)
			    M(i,j) = L[1](i,j-k-1);
			else
			    continue;
		    else
			if(j < k+1)
			    M(i,j) = L[2](i-k-1, j);
			else
			    M(i,j) = L[3](i-k-1, j-k-1);    
	}
	else
	{
	    std::cout<<"b";
	    for(int i = k+2; i < n; i++)
	    {
		M(k, i) -= M(i, k);
		M(i, k) = 0;
	    }
	
	}
    }
    return M;
}

// 2^{-p} reduce the hessenberg matrix M
std::vector<REALMATRIX> split(REALMATRIX M, int p)
{
    REAL abs_value;
    int n = (int) M.maxrow;
    unsigned int prev = 0;
    REAL approx = power(2, p -1) / n;
    std::vector<REALMATRIX> matrices;
    for(int i=1; i<(int) M.maxrow; i++)
    {
	abs_value = abs(M(i, i-1));
	if(choose(abs_value >  approx / 2, abs_value < approx) == 2)
	{
	    matrices.push_back(block_decompose(M, prev, i));
	    prev = i;
	}
    }
    matrices.push_back(block_decompose(M, prev, n));
    return matrices;
} 



/* QR decomposition and QR iterations
 * As matrices we deal with is now tridiagonal matrices,
 * we use Givens rotations to decompose the matrixes.
 */
// Constructs Givens matrix
// i < j <= n, c^2 + s^2 = 1
REALMATRIX Givens (unsigned int n, unsigned int i, unsigned int j, REAL  c, REAL s)
{
    REALMATRIX G = eye(n);
    G(i,i) = c;
    G(j,j) = c;
    G(j,i) = s;
    G(i,j) = - s;
    return G;
}

// QR Decomposition on Hessenberg matrix via Givens rotation
std::pair<REALMATRIX, REALMATRIX> Hessenberg_QR_decomposition(REALMATRIX H)
{
    int n = H.maxrow;
    REALMATRIX Q = eye(n);
    REAL r, c, s;
    REAL tmp1, tmp2;
    for (int i = 0; i < n-1; i++)
    {
	r = sqrt(H(i,i)*H(i,i) + H(i+1, i)*H(i+1, i));
	c = H(i,i) / r;
	s = H(i+1,i) / r;
	for(int l=0; l<n; l++)
	{
	    tmp1 = Q(l,i); tmp2 = Q(l,i+1);
	    Q(l,i) = tmp1*c + tmp2*s;
	    Q(l,i+1) = - tmp1 * s + tmp2 * c;

	    tmp1 = H(i,l); tmp2 = H(i+1,l);
	    H(i,l) = tmp1*c +tmp2*s;
	    H(i+1,l) = - tmp1*s + tmp2*c;
	}    
    }
    return std::pair<REALMATRIX, REALMATRIX>(Q, H);
}

// QR step
REALMATRIX QR(REALMATRIX H)
{
    int n = H.maxrow;
    REALMATRIX Q = eye(n);
    REAL r, c, s;
    REAL tmp1;
    REAL Givens[2*n-2];
    int tmp;
    for (int i = 0; i < n-1; i++)
    {
	r = sqrt(H(i,i)*H(i,i) + H(i+1, i)*H(i+1, i));
	c = H(i,i) / r;
	s = H(i+1,i) / r;
	Givens[i*2] = c;
	Givens[i*2+1] = s;
	tmp = i + 3;
	if(i == n-2)
	    tmp = i + 2;
	for(int l=i; l<tmp; l++)
	{
	    tmp1 = H(i,l); tmp2 = H(i+1,l);
	    H(i,l) = tmp1*c +tmp2*s;
	    H(i+1,l) = - tmp1*s + tmp2*c;
	}    
    }

    for (int i=0; i<n-1; i++)
    {
	for (int j=i-1; j<i+2; j++)
	{
	    if(j == -1 || j == n+1)
		continue;
	    tmp1 = H(j, i) * Givens[2*i] + H(j,i+1) * Givens[2*i + 1];
	    H(j, i+1) = - H(j, i) * Givens[2*i+1] + H(j,i+1) *Givens[2*i]; 
	    H(j,i) = tmp1;
	}
    }

    return H;
}
    

int sign(REAL x, REAL p)
{
    if(choose(x < p, x > - p) == 1)
	return - 1;
    return 1;
}

REALMATRIX Hessenberg_QR_step(REALMATRIX M)
{
    int n = M.maxrow;

    // Soft Wilkinson shift with k = 1/8
    REAL delta = (M(n-1,n-1) - M(n,n))/2;
    REAL s = M(n,n) - sign(delta, M(n,n-1) / 16) * M(n,n-1) *  M(n,n-1) / (abs(delta) + sqrt(delta*delta + M(n,n-1)*M(n,n-1)));

    for(int i=0; i<n; i++)
    {
	M(i,i) = M(i,i) - s;
    }

    std::pair<REALMATRIX, REALMATRIX> QR = Hessenberg_QR_decomposition(M);
    
    REALMATRIX H = mult(QR.second, QR.first);
    
    for(int i=0; i<n; i++)
	H(i,i) = H(i,i) + s;
    
    return H;
}

std::vector<REAL> tridiagonal_QR_algorithm (REALMATRIX M, int p)
{
    int n = M.maxrow;
    std::vector<REAL> eigens = std::vector<REAL>();
    eigens.reserve(n);
    for (int j = 0; j < n; j++)
	eigens.push_back(0);

    REAL error = power(2, p - 1) / n;
    // RATIONAL error = err (-p+1) / n;
    REALMATRIX H = M;
    int m = n - 1;
    for (int i = 0; i < n-2; i++)
    {
	while(choose(abs(H(m,m-1)) > error / 2, abs(H(m,m-1)) < error) == 1)
	{
	    H = Hessenberg_QR_step(H);
	}
	cout <<i<<" ";
	eigens[i] = H(m,m);
	H = block_decompose(H, 0, m);
	m -= 1;
    }
    cout << "\n";

    REAL a1 = H(0,0);
    REAL a2 = H(1,1);
    REAL b = H(1,0);
    REAL tmp = sqrt((a1-a2)*(a1-a2) + 4*b*b);
    REAL l1 = (a1 + a2 + tmp)/2;
    REAL l2 = (a1 + a2 - tmp)/2;
    eigens[n-2] = l1;
    eigens[n-1] = l2;
    return eigens;
}

std::vector<REAL> real_symmetric_eigenvalues(REALMATRIX M, int p)
{

    std::vector<REAL> eigen, eigens;
    int n = M.maxrow;
    int ind = 0;
    eigens.reserve(n);
    for(int i=0; i<n; i++)
	eigens.push_back(0);

    std::cout<<"  Hessenberg reduction start... ";  
    M = hessenberg_reduction(M, p - 3);

    std::cout<<"\n  reducing Hessenberg matrix...\n";
    std::vector<REALMATRIX> T = split(M, p - 3);
  
    std::cout<<"  QR iteration...\n";
    for(unsigned int k = 0; k < T.size(); k++)
    {
	eigen = tridiagonal_QR_algorithm(T[k], p - 3);
	for(unsigned int i = 0; i < T[k].maxrow; i++)
	{
	    eigens[ind] = eigen[i];
	    ind++;
	}
    }
    return eigens;
}


void compute()
{
    std::cout<<"\nrestarting iteration...\n";
    int n = 30;
    int p = - 50;

    REALMATRIX M = gaussian_assymetric(n);
    M = mult(M, M);
    std::cout<<"  random matrix generated...\n";
    
    // measuring runtime start
    std::clock_t c_start = std::clock();
  
    std::vector<REAL> eigens = real_symmetric_eigenvalues(M, p);
    for (int i=0; i<n; i++)
	cout << eigens[i] << "\n";

    // measuring runtime end
    std::clock_t c_end = std::clock();
    long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
}




