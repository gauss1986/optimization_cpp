#include <testR.h>

int main(){
    int N = 10;
    
    double* x = new double[N];
    double* dens = new double[N];

    genRandN(N,x,dens);

    delete[] x;
    delete[] dens;
}
