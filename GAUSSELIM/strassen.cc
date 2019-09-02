#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <vector>
using namespace iRRAM;

void square_block_decompose(const REALMATRIX& x, REALMATRIX* blocks)
{
  int n = x.maxrow;
  REALMATRIX A = REALMATRIX(n/2,n/2);
  REALMATRIX B = REALMATRIX(n/2,n/2);
  REALMATRIX C = REALMATRIX(n/2,n/2);
  REALMATRIX D = REALMATRIX(n/2,n/2);
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<n; j++)
    {
      if(i < n/2)
	if(j < n/2)
	  A(i,j) = x(i,j);
        else
	  B(i, j - n/2) = x(i,j);
      else
	if(j < n/2)
	  C(i - n/2, j) = x(i,j);
	else
	  D(i-n/2,j-n/2) = x(i,j);
    }
  }
  blocks[0] = A;
  blocks[1] = B;
  blocks[2] = C;
  blocks[3] = D;
}

REALMATRIX square_block_compose(REALMATRIX A, REALMATRIX B, REALMATRIX C, REALMATRIX D)
{
  int n = A.maxrow + C.maxrow;
  REALMATRIX composed = REALMATRIX(n,n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      if(i < n/2)
	if(j < n/2)
	  composed(i,j) = A(i,j);
	else
	  composed(i,j) = B(i,j-n/2);
      else
	if(j<n/2)
	  composed(i,j) = C(i-n/2,j);
	else
	  composed(i,j) = D(i-n/2,j-n/2);
  return composed;
}

REALMATRIX add_padding(const REALMATRIX& x)
{
  REALMATRIX padded = zeroes(x.maxrow+1, x.maxcolumn+1);
  for(int i=0; i<x.maxrow; i++)
    for(int j=0; j<x.maxcolumn; j++)
      padded(i,j) = x(i,j);
  return padded;
}

REALMATRIX remove_padding(const REALMATRIX& x)
{
  REALMATRIX unpadded = REALMATRIX(x.maxrow-1, x.maxcolumn-1);
  for(int i=0; i<x.maxrow-1; i++)
    for(int j=0; j<x.maxcolumn-1; j++)
      unpadded(i,j) = x(i,j);
  return unpadded;
}

// Strassen algorithm for where x and y are square matrices
// Currently, phase transform happens for dim = 101
// Later, it should be reduced... (n = 20 sounds good)
REALMATRIX mult (const REALMATRIX& x, const REALMATRIX& y)
{
  int n = x.maxcolumn;
  if(n < 101)
    return x*y;
  bool is_odd = false;
  REALMATRIX A, B;
  if (n % 2 == 1)
  {
    is_odd = true;
    A = add_padding(x);
    B = add_padding(y);
  }
  else
  {
    A = x;
    B = y;
  }
  
  REALMATRIX As[4];
  REALMATRIX Bs[4];
  square_block_decompose(A, As);
  square_block_decompose(B, Bs);

  REALMATRIX M1, M2, M3, M4, M5, M6, M7;
  M1 = mult(As[0] + As[3], Bs[0] + Bs[3]);
  M2 = mult(As[2] + As[3], Bs[0]);
  M3 = mult(As[0], Bs[1] - Bs[3]);
  M4 = mult(As[3], Bs[2] - Bs[0]);
  M5 = mult(As[0] + As[1], Bs[3]);
  M6 = mult(As[2] - As[0], Bs[0] + Bs[1]);
  M7 = mult(As[1] - As[3], Bs[2] + Bs[3]);

  REALMATRIX C = square_block_compose(M1+M4-M5+M7, M3+M5, M2+M4, M1-M2+M3+M6);

  if(is_odd)
  {
    C = remove_padding(C);
  }
  return C;
}

