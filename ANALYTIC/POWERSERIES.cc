# include "POWERSERIES.h"

namespace iRRAM{

POWERSERIES::POWERSERIES() {
	w = COMPLEX(0);
	coef = ([] (int j) -> COMPLEX { return COMPLEX(0);});
	K = k= INTEGER(1);
}
POWERSERIES::POWERSERIES(COEF coef, INTEGER K, INTEGER k, COMPLEX w)
	 : coef(coef), K(K), k(k), w(w) {}
POWERSERIES::POWERSERIES(COMPLEX(*coef)(INTEGER), INTEGER K, INTEGER k, COMPLEX w)
	 : coef(coef), K(K), k(k), w(w) {}

POWERSERIES::POWERSERIES(const POWERSERIES& other) {
	coef = other.coef;
	K = other.K;
	k = other.k;
	w = other.w;
}

COMPLEX POWERSERIES::evalHelper(int p, const COMPLEX& z) {
	COMPLEX result(0);
	COMPLEX pow(1);
	COMPLEX cur;
	COMPLEX dz = z - w;
	INTEGER t = size(K)-p;//size(K) is smallest i s.t. |K| < 2^i

	for(INTEGER i = 0; i<=t;i+=1) {
		cur = (coef(i)*pow);//i-th term of the series 
		result = result + cur;

		pow =pow * dz;
	}
	return result;
}

COMPLEX POWERSERIES::eval(const COMPLEX& z) {
	static POWERSERIES * Tfp; //using capture is not able for iRRAM::limit. so, use static copies
	Tfp = this;
	COMPLEX (*lambda)(int, const COMPLEX& z)  = ([] (int p, const COMPLEX& z) -> COMPLEX{
			return Tfp->evalHelper(p,z);
	}); // to use iRRAM::limit

	return limit(lambda,z);// calculates lim p->inf evalHelper(p,z);
}

POWERSERIES operator +(const POWERSERIES& f1, const POWERSERIES& f2) {
	//ASSERT(f1.w == f2.w)
	//suppose f1 and f2 has same w
	COEF lambda = ([=] (INTEGER j) -> COMPLEX {
			return f1.coef(j)+f2.coef(j);
			}); //result coefficient sequence
	return POWERSERIES(lambda, f1.K+f2.K,std::max(f1.k,f2.k),f1.w);
}

POWERSERIES operator -(const POWERSERIES& f1, const POWERSERIES& f2) {
	//ASSERT(f1.w == f2.w)
	//suppose f1 and f2 has same w
	COEF lambda = ([=] (INTEGER j) -> COMPLEX {
			return f1.coef(j)-f2.coef(j);
			});
	return POWERSERIES(lambda, f1.K+f2.K,std::max(f1.k,f2.k),f1.w);
}
/**
 * @warning given two POWERSERIES should have same z0.
 */
POWERSERIES operator *(const POWERSERIES& f1, const POWERSERIES& f2) {
	//ASSERT(f1.w == f2.w)
	//suppose f1 and f2 has same w
COEF lambda = ([=] (INTEGER j) -> COMPLEX {
					COMPLEX result(0);
					for(INTEGER i = 0;i<=j;i+=1)
						result= result + (f1.coef(i)*f2.coef(j-i));
					return result;
		}	);
	return POWERSERIES(lambda, f1.K*f2.K,f1.k+f2.k,f1.w);
}
POWERSERIES& POWERSERIES::operator +=(const POWERSERIES& f) {
	*this = (*this)+f;
	return *this;
}
POWERSERIES& POWERSERIES::operator -=(const POWERSERIES& f) {
	*this = (*this)-f;
	return *this;
}
POWERSERIES& POWERSERIES::operator *=(const POWERSERIES& f) {
	*this = (*this)*f;
	return *this;
}
POWERSERIES& POWERSERIES::operator=(const POWERSERIES& other) {
	if(this == (&other))
		return (*this);
	coef = other.coef;
	K = other.K;
	k = other.k;
	w = other.w;
	return (*this);
}

POWERSERIES POWERSERIES::differentiateHelper(INTEGER d) {
	COEF lambda = ([=] (INTEGER j) -> COMPLEX {
			INTEGER diffTerm(1);
			for(INTEGER i = 1;i<=d;i+=1)
			    diffTerm*=(i+j);
			return (coef(j+d))*diffTerm;
			});

	INTEGER temp=1;
	for(INTEGER i = 1;i<=d;i+=1)
			temp*=(d+i);
	temp=temp<<d;
	INTEGER newK = K*(k^d)*temp;

	return POWERSERIES(lambda, newK, k<<1, w);//(k<<1 == 2k)
}

//TODO
POWERSERIES POWERSERIES::integralHelper(INTEGER d) {
	COEF lambda = ([=] (INTEGER j) -> COMPLEX {
			if(j<d)
				return COMPLEX(0);

			INTEGER intTerm(1);
			for(INTEGER i = 0;i<d;i+=1)
					intTerm*=(j-i);
			return (coef(j-d))/intTerm;
			});
	return POWERSERIES(lambda, K, k, w);
}

POWERSERIES POWERSERIES::differentiate(INTEGER d) {
	if(d==0)
		return *this;
	if(d>0)
		return this->differentiateHelper(d);
	else
		return this->integralHelper(-d);
}

} //namespace iRRAM
