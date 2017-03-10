#include <iostream>
#include <Rmath.h>
#include <R.h>

int main(){
    int N = 10;
    double x;
    double dens;   
 
    //double* x = new double[N];
    //double* dens = new double[N];

    //genRandN(N,x,dens);
    GetRNGstate();
    for (int i=0;i<N;i++){
        x=rnorm(1.0,3.0);
        //dens=dnorm(x,1.0,3.0,true);
        std::cout << x << "," << dens << std::endl; 
    }
    PutRNGstate();

    //delete[] x;
    //delete[] dens;
}
