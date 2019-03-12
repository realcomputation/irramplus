
#include <random-real.h>

using namespace iRRAM;
using namespace random;

REAL random::uniform_real()
{
  return REALRAND().asREAL();
}

REAL random::uniform_real(REAL a, REAL b)
{
  if(b > a)
    return uniform_real()*(b-a) + a ;
  else
    return 1;
}

REAL random::gaussian_real()
{
  REAL V = uniform_real();
  REAL U = uniform_real();
  return sqrt(-REAL(2)*log(U)/log(exp(1)))*cos(2*pi()*V);
}

REAL random::gaussian_real(REAL exp, REAL std)
{
  if(std > 0)
    return gaussian_real()*std + exp;
  else
    return 1;
}


// from probabilistic density 2x
REAL random::linear_real()
{
  return maximum(uniform_real(), uniform_real());
}

