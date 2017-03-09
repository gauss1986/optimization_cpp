#include <Rmath.h>
#include <R.h>

void genRandN (int size, double *x, double *dens){
    GetRNGstate();
    for (int i=0;i<size;i++){
        x[i]=rnorm(1.0,3.0);
        dens[i]=dnorm(x[i],1.0,3.0,true);
    }
    PutRNGstate();
}
