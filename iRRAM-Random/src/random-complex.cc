
#include <utility>

#include <random-complex.h>

using namespace iRRAM;
using namespace random;

COMPLEX random::uniform_complex()
{
  REAL magnitude = linear_real();
  REAL angle = uniform_real(0, 2*pi());
  return COMPLEX(magnitude*cos(angle), magnitude*sin(angle));
}

COMPLEX random::uniform_complex(COMPLEX c, REAL r)
{
  return uniform_complex()*r + c;
}

COMPLEX random::gaussian_complex()
{
  return COMPLEX(gaussian_real(), gaussian_real());
}
