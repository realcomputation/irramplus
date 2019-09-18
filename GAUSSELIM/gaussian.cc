#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"
#include "random-matrix.h"
using namespace iRRAM;

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


// apply gaussian elimination to a regular square matrix A
REALMATRIX gaussian_elimination(REALMATRIX A)
{

  int n = A.maxrow;
  int pi, pj;
  REAL m;
  REAL temp;
  REAL div[n];
  for(int i=0; i<n-1; i++)
    {
      m = 0;
      for(int j=i; j<n; j++)
	for(int k=i; k<n; k++)
	  m = maximum(m, abs(A(j,k)));
      for(int j=i;j<n;j++)
	for(int k=i; k<n; k++)
	  if (choose(abs(A(j,k)) > m / 2, abs(A(j,k)) < m) == 1)
	    { pi = j; pj = k;}

      for(int j=0; j<n; j++)
	{
	  temp = A(j,i);
	  A(j,i) = A(j, pj);
	  A(j,pj) = temp;
	}
      for(int j=0; j<n; j++)
	{
	  temp = A(i,j);
	  A(i,j) = A(pi, j);
	  A(pi, j) = temp;
	}
      for(int j=i+1; j<n; j++)
	div[j] = A(i,j)/A(i,i);
      div[i] = 1;

      for(int j=i+1; j<n; j++)
	{
	  for(int k=i+1; k<n; k++)
	    A(j,k) = A(j,k) - div[k] * A(j,i);
	  A(j,i) = 0;
	}

    }
  return A;
}

REAL gaussian_minor(REALMATRIX A, int k, double alpha)
{

  	int n = A.maxrow;


  	int pi, pj, tp;
  	REAL m;
  	REAL temp, tmp;
  	REAL div[n];

  	REAL sort[n*n];
  	int sort_order[n*n];
  	int sortx[n*n];
  	int sorty[n*n];
  	int sort_tmp;

	for(int i=0; i<k; i++)
	{	
		tp=0;
	    for(int j=i; j<n; j++)
	    {
			for(int k=i; k<n; k++)
		  	{
		  		sort[tp] = abs(A(j,k));
		  		sort_order[tp] = tp;
		  		sortx[tp] = j; 
		  		sorty[tp] = k;
		  		tp ++;
		  	}
		}
		for(int ix=0; ix<tp-1; ix++)
			for(int iy=ix+1; iy<tp; iy++)
			{
				if(sort[ix] < sort[iy])
				{
					tmp = sort[ix];
					sort[ix] = sort[iy];
					sort[ix] = tmp;
					sort_tmp = sortx[ix];
					sortx[ix] = sorty[iy];
					sorty[iy] = sort_tmp;

					sort_tmp = sorty[ix];
					sorty[ix] = sorty[iy];
					sorty[iy] = sort_tmp;
				}
			}


		pi = sortx[(int) (alpha * double((n-i)*(n-i)))];
		pj = sorty[(int) (alpha * double((n-i)*(n-i)))];
	      for(int j=0; j<n; j++)
		{
		  temp = A(j,i);
		  A(j,i) = A(j, pj);
		  A(j,pj) = temp;
		}
	      for(int j=0; j<n; j++)
		{
		  temp = A(i,j);
		  A(i,j) = A(pi, j);
		  A(pi, j) = temp;
		}
	      for(int j=i+1; j<n; j++)
		div[j] = A(i,j)/A(i,i);
	      div[i] = 1;

	      for(int j=i+1; j<n; j++)
		{
		  for(int k=i+1; k<n; k++)
		    A(j,k) = A(j,k) - div[k] * A(j,i);
		  A(j,i) = 0;
		}

	    }
	    m = 1;
	    for(int i=0; i<k; i++)
	    	m = m * abs(A(i,i));
	  return m;


}

REAL experiment()
{
	REALMATRIX A = gaussian_matrix(10);
	cout << gaussian_minor(A, 10, 0) << ", " << gaussian_minor(A,10,0.9) << "\n";
}


void compute()
{
	experiment();
	exit(0);
  	sizetype error, error_f;   
  	int errorexp, errorexp_f; 
  	REALMATRIX A, B;
  	int n = 50;
	A = 100*gaussian_symmetric_matrix(n);
  	A.geterror(error);
	errorexp =  (int)std::floor(std::log2(error.mantissa)) + error.exponent;
  	B = A*A*A; //= gaussian_elimination(A);
  	B.geterror(error_f);
	errorexp_f =  (int)std::floor(std::log2(error_f.mantissa)) + error_f.exponent;
	std::cout << errorexp << ", " << errorexp_f << ", "<< errorexp + (int)std::floor(std::log2(11)*n) << "\n";

	if(REAL(1) > REAL(1))
		exit(1);

}
