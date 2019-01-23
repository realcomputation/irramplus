REAL conditional(LAZY_BOOLEAN b, const REAL& x, const REAL& y) 
{ // return:
// x if b==TRUE,
// y if b==FALSE
// x (or y) if  b==UNKNOWN  and x==y
// UNDEFINED  if b==UNKNOWN and  x<>y
    switch (b.value) {
       case 1:
           return x;
       case 0:
           return y;
       default:
           if ( iRRAM_unlikely ( x.value || y.value ) ) {
               REAL z = (x+y)/2;
               REAL delta = abs(x-y)/2;
               if ( !delta.value ) delta.mp_make_mp();
               z.adderror(delta.error);
               z.adderror(delta.vsize);
               // this is not optimal:  
               // in this solution, the error contains (x.error+y.error)
               // while  max(x.error,y.error) would suffice
               return z;
           }       
           return REAL(
                   iRRAM_double_pair(
                       fmin(x.dp.lower_pos,y.dp.lower_pos),
                       fmin(x.dp.upper_neg,y.dp.upper_neg)
                                    )
                   );
   }
}

