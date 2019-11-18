# include "ANALYTIC.h"

namespace iRRAM{
ANALYTIC::ANALYTIC()
{
	coef = ([] (INTEGER i) -> COMPLEX { return(0);});//constant zero function
	l = L = 1;
	x = REAL(0);
}

ANALYTIC::ANALYTIC(COEF coef, INTEGER L, INTEGER l, REAL x)
	 : coef(coef), L(L), l(l), x(x) {}
ANALYTIC::ANALYTIC(COMPLEX(* coef)(INTEGER), INTEGER L, INTEGER l, REAL x)
	 : coef(coef), L(L), l(l), x(x) {}

ANALYTIC::ANALYTIC(const ANALYTIC& other) {
	coef = other.coef;
	L = other.L;
	l = other.l;
	x = other.x;
}

COMPLEX ANALYTIC::evalHelper(int p, const COMPLEX& z) {
	COMPLEX result(0);
	COMPLEX pow(1);
	COMPLEX cur;
	COMPLEX dz = z - x;
	INTEGER t = size(L)-p;//size(L) is smallest i s.t. |L| < 2^i

	for(INTEGER i = 0; i<=t;i+=1) {
		cur = (coef(i)*pow);//i-th term of the series 
		result = result + cur;

		pow =pow * dz;
	}
	return result;
}

COMPLEX ANALYTIC::eval(const COMPLEX& z) {
	static ANALYTIC * Tfp; //using capture is not able for iRRAM::limit. so, use static copies
	Tfp = this;
	COMPLEX (*lambda)(int, const COMPLEX&)  = ([] (int p, const COMPLEX& z) -> COMPLEX{
			return Tfp->evalHelper(p,z);
	}); // to use iRRAM::limit

	return limit(lambda,z);// calculates lim p->inf evalHelper(p,z);
}

ANALYTIC operator +(const ANALYTIC& f1, const ANALYTIC& f2) {
	//ASSERT(f1.x == f2.x)
	//suppose f1 and f2 has same x
	COEF lambda = ([=] (INTEGER j) -> COMPLEX {
			return f1.coef(j)+f2.coef(j);
			}); //result coefficient sequence
	return ANALYTIC(lambda, f1.L+f2.L,std::max(f1.l,f2.l),f1.x);
}

ANALYTIC operator -(const ANALYTIC& f1, const ANALYTIC& f2) {
	//ASSERT(f1.x == f2.x)
	//suppose f1 and f2 has same x
	COEF lambda = ([=] (INTEGER j) -> COMPLEX {
			return f1.coef(j)-f2.coef(j);
			});
	return ANALYTIC(lambda, f1.L+f2.L,std::max(f1.l,f2.l),f1.x);
}
ANALYTIC operator *(const ANALYTIC& f1, const ANALYTIC& f2) {
	//ASSERT(f1.x == f2.x)
	//suppose f1 and f2 has same x
COEF lambda = ([=] (INTEGER j) -> COMPLEX {
					COMPLEX result(0);
					for(INTEGER i = 0;i<=j;i+=1)
						result= result + (f1.coef(i)*f2.coef(j-i));
					return result;
		}	);
	return ANALYTIC(lambda, f1.L*f2.L,std::max(f1.l,f2.l),f1.x);
}
ANALYTIC& ANALYTIC::operator +=(const ANALYTIC& f) {
	*this = (*this)+f;
	return *this;
}
ANALYTIC& ANALYTIC::operator -=(const ANALYTIC& f) {
	*this = (*this)-f;
	return *this;
}
ANALYTIC& ANALYTIC::operator *=(const ANALYTIC& f) {
	*this = (*this)*f;
	return *this;
}
ANALYTIC& ANALYTIC::operator=(const ANALYTIC& other) {
	if(this == (&other))
		return (*this);
	coef = other.coef;
	L = other.L;
	l = other.l;
	x = other.x;
	return (*this);
}

ANALYTIC ANALYTIC::differentiateHelper(INTEGER d) {
	COEF lambda = ([=] (INTEGER j) -> COMPLEX {
			INTEGER diffTerm(1);
			for(INTEGER i = 1;i<=d;i+=1)
			    diffTerm*=(i+j);
			return (coef(j+d))*diffTerm;
			});
	INTEGER factorial(1);
	INTEGER newl = l<<1;
	for(INTEGER i = 2;i<=d;i+=1){
		factorial*=i;
	}
	INTEGER newL = L* factorial * (newl^d);
	

	return ANALYTIC(lambda, newL, newl, x);//(l<<1 == 2l)
}

//TODO
ANALYTIC ANALYTIC::integralHelper(INTEGER d) {
	COEF lambda = ([=] (INTEGER j) -> COMPLEX {
			if(j<d)
				return COMPLEX(0);

			INTEGER intTerm(1);
			for(INTEGER i = 0;i<d;i+=1)
					intTerm*=(j-i);
			return (coef(j-d))/intTerm;
			});
	return ANALYTIC(lambda, L, l, x);
}

ANALYTIC ANALYTIC::differentiate(INTEGER d) {
	if(d==0)
		return *this;
	if(d>0)
		return this->differentiateHelper(d);
	else
		return this->integralHelper(-d);
}
ANALYTIC ANALYTIC::continuation(REAL newX) {
	COEF lambda = ([=] (INTEGER j) -> COMPLEX {
		INTEGER factorial(1);
		for(INTEGER i = 2;i<=j;i+=1)
			factorial*=j;
		return this->differentiate(j).eval(newX)/factorial;

	});
	return ANALYTIC(lambda, L, l, newX);
}

} //namespace iRRAM
