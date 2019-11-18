# pragma once

# include "iRRAM/lib.h"
# include <functional>
# include <algorithm>


namespace iRRAM {

using COEF = std::function<COMPLEX(INTEGER)>;
class POWERSERIES
{
public:
POWERSERIES();
POWERSERIES(COEF coef, INTEGER K, INTEGER k, COMPLEX w = COMPLEX(0));
POWERSERIES(COMPLEX(*coef)(INTEGER), INTEGER K, INTEGER k, COMPLEX w = COMPLEX(0));
POWERSERIES(const POWERSERIES& other);

friend POWERSERIES operator +(const POWERSERIES& f1, const POWERSERIES& f2);
friend POWERSERIES operator -(const POWERSERIES& f1, const POWERSERIES& f2);
friend POWERSERIES operator *(const POWERSERIES& f1, const POWERSERIES& f2);

POWERSERIES& operator+=(const POWERSERIES& f);
POWERSERIES& operator-=(const POWERSERIES& f);
POWERSERIES& operator*=(const POWERSERIES& f);

POWERSERIES& operator=(const POWERSERIES& other);

COMPLEX eval(const COMPLEX& z);//only for nonnegative d.
POWERSERIES differentiate(INTEGER d);//also for negative d. and return function has integral constant all 0.
private:

POWERSERIES differentiateHelper(INTEGER d);//for only differentiation (d>0)
POWERSERIES integralHelper(INTEGER d);// for only integral (d>0)
COMPLEX evalHelper(int p, const COMPLEX& z);

COMPLEX w;
COEF coef;
INTEGER K, k;

};

}//namespace iRRAM

