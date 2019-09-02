#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "irram-random.h"
#include <vector>
#include <iostream>
#include <ctime>
using namespace iRRAM;


FASTMATRIX::FASTMATRIX(unsigned int rows,unsigned int columns)
{
  maxrow=rows;maxcolumn=columns;
};

FASTMATRIX& FASTMATRIX::operator = (const FASTMATRIX& y) {

  core = y.core;
  maxrow = y.maxrow;
  maxcolumn = y.maxcolumn;
  return(*this);
}

REALMATRIX operator + (const REALMATRIX& x, const REALMATRIX& y) {
  REALMATRIX sum(x.maxrow,x.maxcolumn);
  
  if ( x.maxrow != y.maxrow || x.maxcolumn != y.maxcolumn) { 
     fprintf(stderr,"Error in adding real matrices of different sizes \n");
     exit(1);
     }
  unsigned int size=x.maxrow*x.maxcolumn;
  for (unsigned int i=0; i<size;i++) sum.values[i]=x.values[i]+y.values[i];
  return sum;
}

REALMATRIX operator - (const REALMATRIX& x, const REALMATRIX& y) {
  REALMATRIX sum(x.maxrow,x.maxcolumn);
  
  if ( x.maxrow != y.maxrow || x.maxcolumn != y.maxcolumn) { 
     fprintf(stderr,"Error in adding real matrices of different sizes \n");
     exit(1);
     }
  unsigned int size=x.maxrow*x.maxcolumn;
  for (unsigned int i=0; i<size;i++) sum.values[i]=x.values[i]-y.values[i];
  return sum;
}

REALMATRIX operator * (const REALMATRIX& x ,const REALMATRIX& y) {
  REALMATRIX prod(x.maxrow,y.maxcolumn);
  
  if ( x.maxcolumn != y.maxrow) { 
     fprintf(stderr,"Error in multiplying real matrices of different sizes \n");
     exit(1);
     }
  for (unsigned int i=0; i<x.maxrow;i++) {
  for (unsigned int j=0; j<y.maxcolumn;j++) {
    ELEMENT(prod,i,j)=ELEMENT(x,i,0)*ELEMENT(y,0,j);
    for (unsigned int k=1; k<x.maxcolumn;k++) {
      ELEMENT(prod,i,j)=ELEMENT(prod,i,j)+ELEMENT(x,i,k)*ELEMENT(y,k,j);
    } } }
  return prod;
}

REALMATRIX operator / (const REALMATRIX& x,const REALMATRIX& y) {
  REALMATRIX rside=x,lside=y,result;
  return solve(lside,rside,1);
}

REALMATRIX solve
       (REALMATRIX& lside,
        REALMATRIX& rside,
        int               use_pivot) {
continous_begin ();
  REALMATRIX result(rside.maxrow,rside.maxcolumn);
  REAL faktor;
  int testsize,pivotsize,pivot,pivotsave;
  std::vector<int> perm (lside.maxrow);
  
  if ( lside.maxrow != rside.maxrow) { 
     fprintf(stderr,"Error, incompatible real matrices\n");
     exit(1);
     }
  for (unsigned int i=0; i<lside.maxrow;i++) perm[i]=i;
  for (unsigned int i=0; i < lside.maxrow; i++) {
    pivot=i;
    if (use_pivot) {
      pivotsize=upperbound(ELEMENT(lside,perm[i],i));
      for (unsigned int j=i+1; j < lside.maxcolumn; j++) {
        testsize=upperbound(ELEMENT(lside,perm[j],i));
        if (testsize>pivotsize) {
          pivot=j;
          pivotsize=testsize;
        }
        else if ( testsize < pivotsize && use_pivot == 1) break ;
      }
      pivotsave=perm[pivot];
      perm[pivot]=perm[i];
      perm[i]=pivotsave;  
    }
    for (unsigned int j=i+1; j < lside.maxcolumn; j++) {
      faktor=ELEMENT(lside,perm[j],i)/ELEMENT(lside,perm[i],i);
      for (unsigned int k=0; k< rside.maxcolumn; k++) 
         ELEMENT(rside,perm[j],k)=ELEMENT(rside,perm[j],k)
                   -faktor*ELEMENT(rside,perm[i],k);
      for (unsigned int l=i; l< lside.maxcolumn; l++) 
         ELEMENT(lside,perm[j],l)=ELEMENT(lside,perm[j],l)
                   -faktor*ELEMENT(lside,perm[i],l);
    }
  }

  for (int i=rside.maxrow-1; i >= 0; i--) {
     faktor=1/ELEMENT(lside,perm[i],i);
     for (unsigned int k=0; k< rside.maxcolumn; k++) {
   ELEMENT(result,i,k)=ELEMENT(rside,perm[i],k)*faktor;
   for (int j=0; j < i; j++) 
     ELEMENT(rside,perm[j],k)=ELEMENT(rside,perm[j],k)
       -ELEMENT(lside,perm[j],i)* ELEMENT(result,i,k);
  } }
continous_end ();
  return result;
}

REALMATRIX::~REALMATRIX() {
//fprintf(stderr," +- Destructor 1m: %x\n",this);
  delete []values;
};

REALMATRIX::REALMATRIX(){
  maxrow=0;maxcolumn=0;values=NULL;
};

REALMATRIX::REALMATRIX(const REALMATRIX& y){
  maxrow=y.maxrow;maxcolumn=y.maxcolumn;
  unsigned int size=maxrow*maxcolumn;
  values = new REAL[size];
  for (unsigned int i=0; i < size; i++) values[i]=y.values[i];
};

REAL&  REALMATRIX::element (unsigned int i, unsigned int j) const {
  if ((i>= maxrow) || (j>= maxcolumn) ) {
     fprintf(stderr,"Illegal indices [%d,%d] for real matrix [%d,%d]\n",
           i,j,maxrow,maxcolumn);
     exit(1);
     }
  return values[i*maxcolumn+j];
};

REAL&  REALMATRIX::operator () (unsigned int i, unsigned int j) const {
  if ((i>= maxrow) || (j>= maxcolumn) ) {
     fprintf(stderr,"Illegal indices [%d,%d] for real matrix [%d,%d]\n",
           i,j,maxrow,maxcolumn);
     exit(1);
     }
  return values[i*maxcolumn+j];
};


REALMATRIX operator * (const REALMATRIX& x, const REAL& y) {
  REALMATRIX prod(x.maxrow,x.maxcolumn); 
  for (unsigned int i=0; i<x.maxcolumn;i++) {
  for (unsigned int j=0; j<x.maxrow;j++) {
      ELEMENT(prod,i,j)=ELEMENT(x,i,j)*y;
    } }
  return prod;
}

REALMATRIX operator / (const REALMATRIX& x, const REAL& y) {
  REALMATRIX prod(x.maxrow,x.maxcolumn); 
  for (unsigned int i=0; i<x.maxcolumn;i++) {
  for (unsigned int j=0; j<x.maxrow;j++) {
      ELEMENT(prod,i,j)=ELEMENT(x,i,j)/y;
    } }
  return prod;
}

int bound (const REALMATRIX& x, const int k){
  for (unsigned int i=0;i<x.maxcolumn;i++) 
  for (unsigned int j=0;j<x.maxrow;j++) 
  if ( ! bound(ELEMENT(x,i,j),k) ) return 0;
  return 1;
}

void REALMATRIX::adderror (sizetype error)
{ 
  for (unsigned int i=0;i<(*this).maxcolumn;i++) 
  for (unsigned int j=0;j<(*this).maxrow;j++) 
      ELEMENT((*this),i,j).adderror(error);
}

void REALMATRIX::seterror (sizetype error)
{ 
  for (unsigned int i=0;i<(*this).maxcolumn;i++) 
  for (unsigned int j=0;j<(*this).maxrow;j++) 
      ELEMENT((*this),i,j).seterror(error);
}

void REALMATRIX::geterror (sizetype& error) const
{
  unsigned int i,j;
  sizetype lerror; 
  ELEMENT((*this),0,0).geterror(error); 
  for (i=0;i<(*this).maxcolumn;i++) 
  for (j=0;j<(*this).maxrow;j++) {
      ELEMENT((*this),i,j).geterror(lerror);
      sizetype_max(error,error,lerror);
  }
}


REAL       maxnorm (const REALMATRIX& x){
  REAL m=abs(x.values[0]);
  unsigned int size=x.maxrow*x.maxcolumn;
  for (unsigned int i=1; i<size;i++) m=maximum(m,abs(x.values[i]));
  return m;
}

REAL       rowsumnorm (const REALMATRIX& x){
  REAL m=0;
  REAL sum;
  for (unsigned int i=0;i<x.maxrow;i++) {
    sum=0;
    for (unsigned int j=0;j<x.maxcolumn;j++) sum=sum+abs(x(i,j));
    m=maximum(m,sum);
  }
  return m;
}

REAL       colsumnorm (const REALMATRIX& x){
  REAL m=0;
  REAL sum;
  for (unsigned int j=0;j<x.maxcolumn;j++) {
    sum=0;
    for (unsigned int i=0;i<x.maxrow;i++) sum=sum+abs(x(i,j));
    m=maximum(m,sum);
  }
  return m;
}
