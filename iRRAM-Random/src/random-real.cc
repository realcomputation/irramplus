
#include <random-real.h>

using namespace iRRAM;
using namespace random;

typename random_device::result_type random_device::operator()()
{
	result_type r;
	state_t &st = *state;
	if (get_cached(r, st))
		return r;
	r = std::random_device::operator()();
	put_cached(r, st);
	return r;
}

REAL random::uniform_real(unsigned int seed)
{
  std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
  // random int generator (uniformly distributed [0,2^16-1])
  std::uniform_int_distribution<> dis(0, 65535);

  std::string r;
  int s;
  if (!(get_cached(r) && get_cached(s)))
  {
//  in continuous section, return (inconsistent) random number with no using cache
    r = std::to_string(dis(gen));
    s = 16;
    put_cached(r);
    put_cached(s);
  }

  int prec = state->ACTUAL_STACK.actual_prec;

  INTEGER result_integer = r;
  int needblock = -prec/16 + 1;
  int neededblock = needblock - s/16;
  for(int i=0; i<neededblock; i++)
    result_integer = (result_integer << 16) + dis(gen);

  std::string bits = swrite(result_integer);
  int bitlength = needblock * 16;

  modify_cached(bits);
  modify_cached(bitlength);

  sizetype err;
  sizetype_set(err,1,-bitlength);
  REAL randreal = scale(RATIONAL(bits), -bitlength);
  randreal.adderror(err);

  return randreal;
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
