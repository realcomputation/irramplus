#include "iRRAM.h"

#include <stdlib.h>     
#include <random>
#include <vector>
#include <string> 
using namespace iRRAM;


INTEGER pow(int a, int b)
{
  INTEGER result = 1;
  for(int i=0; i< b; i++)
    result = result * a;
  return result;
}


class REALRAND
{
  private:
    int bitlength;
    std::string bits;
    REAL randreal;

  public:
    REALRAND()
    {
      int prec = state.ACTUAL_STACK.actual_prec;

      sizetype err;
      std::string r;
      int s;
      INTEGER result_integer;
      int needblock, neededblock;
      std::string result_bits;
      int result_size ;



      // random int generator (uniformly distributed [0,2^16-1])
      std::random_device rd;  //Will be used to obtain a seed for the random number engine
      std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
      std::uniform_int_distribution<> dis(0, 65535);


      if (state.ACTUAL_STACK.inlimit == 0)
      {
        if(state.thread_data_address->cache_s.get(result_bits)
          && state.thread_data_address->cache_i.get(result_size))
        {
          r = result_bits;
          s = result_size;
        }
        else
        {
          result_bits = std::to_string(dis(gen));
          result_size = 16;
          state.thread_data_address->cache_s.put(result_bits);
          state.thread_data_address->cache_i.put(result_size);
          r = result_bits;
          s = result_size;
        }

      INTEGER result_integer = r;
      int needblock = (-prec)/16 +1;
      int neededblock = needblock - s/16;


      for(int i=0; i<neededblock; i++)
        result_integer = result_integer*pow(2,16) + dis(gen);


      bits = swrite(result_integer);
      bitlength = needblock * 16;

      state.thread_data_address->cache_s.modify(bits);
      state.thread_data_address->cache_i.modify(bitlength);


      sizetype_set(err,1,-bitlength);
      randreal = RATIONAL(INTEGER(bits), pow(2,bitlength));
      randreal.adderror(err);
    }


//  in continuous section, return (inconsistent) random number with no using chache
    else
    {
      result_bits = std::to_string(dis(gen));
      result_size = 16;
      r = result_bits;
      s = result_size;
        
      INTEGER result_integer = r;
      needblock = (-prec)/16 +1;
      neededblock = needblock - s/16;


      for(int i=0; i<neededblock; i++)
        result_integer = result_integer * pow(2,16) + dis(gen);


      bits = swrite(result_integer);
      bitlength = needblock * 16;

      sizetype_set(err,1,-bitlength);
      randreal = RATIONAL(INTEGER(bits), pow(2,bitlength));
      randreal.adderror(err);
    }
  }


    ~REALRAND()
    {
    }

    REAL asREAL()
    {
      return randreal;
    }
};

REAL gaussian_real()
{
  REAL V = REALRAND().asREAL();
  REAL U = REALRAND().asREAL();
  return sqrt(-REAL(2)*log(U)/log(exp(1)))*cos(2*pi()*V);
}

REAL uniform_real()
{
  return REALRAND().asREAL();
}

REALMATRIX gaussian_orthogonal_ensemble(unsigned int n)
{
  REALMATRIX M = REALMATRIX(n, n);
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < n; j++)
      if(i>j) M(i,j) = M(j,i); else M(i, j) = gaussian_real();
  return M;
}

REALMATRIX gaussian_assymetric(unsigned int n)
{
  REALMATRIX M = REALMATRIX(n, n);
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < n; j++)
      if(i<j) M(i,j) = gaussian_real(); else M(i,j) = - M(j,i); 
  return M;
}
