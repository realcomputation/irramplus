#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "iRRAM.h"
#include "complement.h"


using namespace iRRAM;


class Grassmannian{
	private:
		unsigned dimension;
	public:
		REALMATRIX B;
};


void write(const REALMATRIX& M){
	for(unsigned i = 0; i < rows(M); ++i){
		for(unsigned j = 0; j < columns(M); ++j){
			cout << M(i, j) << " ";
		}
		cout << "\n";
	}
}

void gaussElim(REALMATRIX& A, const unsigned iters){
	unsigned m = rows(A), n = columns(A);
	unsigned pj;
	REAL absMax;
	for(unsigned i = 0; i < iters; ++i){
		pj = 0;
		absMax = 0;
		for(unsigned j = i; j < n; ++j){
			absMax = maximum(absMax, abs(A(i, j)));
		}
		for(unsigned j = i; j < n; ++j){
			if(choose(abs(A(i, j)) > absMax/2, abs(A(i, j)) < absMax) == 1){
				pj = j; // select the column with a nonzero entry
			}
		}
		REAL tmp;
		for(unsigned k = i; k < m; ++k){
			tmp = A(k, i);
			A(k, i) = A(k, pj);
			A(k, pj) = tmp; // swap columns i and pj
		}
		for(unsigned j = i+1; j < n; ++j){
			REAL multiplier = A(i, j) / A(i, i);
			for(unsigned k = i; k < m; ++k){
				A(k, j) -= multiplier * A(k, i);
			}
		}
	}
}


REALMATRIX complement(const REALMATRIX& A){
	unsigned d = rows(A), m = columns(A);
	REALMATRIX M(m+d, d);
	for(unsigned i = 0; i < m; ++i){
		for(unsigned j = 0; j < d; ++j){
			M(i, j) = A(j, i); // copy over entries of A transposed
		}
	}
	for(unsigned i = m; i < m+d; ++i){
		for(unsigned j = 0; j < d; ++j){
			if(i == j+m){
				M(i, j) = 1;
			}else{
				M(i, j) = 0; // copy over entries of I_d
			}
		}
	}
	gaussElim(M, m);
	REALMATRIX result(d, d-m);
	for(unsigned i = 0; i < d; ++i){
		for(unsigned j = 0; j < d-m; ++j){
			result(i, j) = M(i+m, j+m);
		}
	}
	return result;
}

REALMATRIX join(const REALMATRIX &A, const REALMATRIX& B, const unsigned l){
	unsigned d = rows(A), m = columns(A), n = columns(B);
	REALMATRIX M(2*d, m+n);
	for(unsigned i = 0; i < d; ++i){
		for(unsigned j = 0; j < m; ++j){
			M(i, j) = A(i, j);
		}
		for(unsigned j = m; j < m+n; ++j){
			M(i, j) = B(i, j-m);
		}
	}
	for(unsigned i = d; i < 2*d; ++i){
		for(unsigned j = 0; j < m; ++j){
			M(i, j) = A(i-d, j);
		}
		for(unsigned j = m; j < m+n; ++j){
			M(i, j) = 0;
		}
	}
	gaussElim(M, m+n);
	REALMATRIX result(d, l);
	for(unsigned i = 0; i < d; ++i){
		for(unsigned j = 0; j < l; ++j){
			result(i, j) = M(i, j);
		}
	}
	return result;
}


REALMATRIX meet(const REALMATRIX &A, const REALMATRIX& B, const unsigned k){
	unsigned d = rows(A), m = columns(A), n = columns(B);
	REALMATRIX M(2*d, m+n);
	for(unsigned i = 0; i < d; ++i){
		for(unsigned j = 0; j < m; ++j){
			M(i, j) = A(i, j);
		}
		for(unsigned j = m; j < m+n; ++j){
			M(i, j) = B(i, j-m);
		}
	}
	for(unsigned i = d; i < 2*d; ++i){
		for(unsigned j = 0; j < m; ++j){
			M(i, j) = A(i-d, j);
		}
		for(unsigned j = m; j < m+n; ++j){
			M(i, j) = 0;
		}
	}
	gaussElim(M, m+n);
	REALMATRIX result(d, k);
	for(unsigned i = 0; i < d; ++i){
		for(unsigned j = 0; j < k; ++j){
			result(i, j) = M(i+m, j+m+n-k);
		}
	}
	return result;
}


REALMATRIX project(const REALMATRIX& A, const REALMATRIX& B, const unsigned l){
	unsigned d = rows(A), m = columns(A), n = columns(B);
	REALMATRIX T(m, d);
	for(unsigned i = 0; i < m; ++i){
		for(unsigned j = 0; j < d; ++j){
			T(i, j) = A(j, i);
		}
	}
	REALMATRIX P(d, d);
	REALMATRIX X(m, m);
	X = eye(m)/(T*A);
	P = A*X*T;
	REALMATRIX Q(d, n);
	Q = P*B;
	gaussElim(Q, l);
	REALMATRIX result(d, l);
	for(unsigned i = 0; i < d; ++i){
		for(unsigned j = 0; j < l; ++j){
			result(i, j) = Q(i, j);
		}
	}
	return result;
}


void compute(){
	REALMATRIX p(4, 2);
	p(0, 0) = 1;
	p(1, 0) = 3;
	p(2, 0) = -1;
	p(3, 0) = 1;
	p(0, 1) = 0;
	p(1, 1) = 1;
	p(2, 1) = 1;
	p(3, 1) = 6;
	REALMATRIX q(2, 2);
	q = eye(2);
	// write(complement(p));
	REALMATRIX A(3, 2);
	A(0, 0) = 1;
	A(0, 1) = 0;
	A(1, 0) = 0;
	A(1, 1) = 1;
	A(2, 0) = 0;
	A(2, 1) = 0;
	REALMATRIX B(3, 1);
	B(0, 0) = 3;
	B(1, 0) = 2;
	B(2, 0) = 1;
	write(project(A, B, 1));
}












