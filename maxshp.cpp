#include <mkl_lapacke.h>
#include <stdlib.h>
#include <stdio.h>
#include "armadillo"
#include <maxshp.h>

using namespace arma;

void maxshp(const int N, const int n, const int m, std::vector<std::vector<double> >& x0, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y){
    // covert x0,x,y to arma format
    mat mx(N,n+1);
    mat my(N,n+1);
    vec vx0(N);
    for (int i=0;i<N;i++){
        for (int j=0;j<n+1;j++){
            mx(i,j) = x[i][j];
            my(i,j) = y[i][j];
        }
        vx0(i) = x0[i][0];
    }

    // construct xy
    mat xy(N,n+1);
    xy.col(0) = sum(my,1);
    xy.col(1) = vx0%xy.col(0);
    xy.col(2) = sum(mx%my,1);

    //std::cout << "Cov(xy):" << std::endl;
    //cov(xy).print();
    //std::cout << "Cov(xy)^(-1):" << std::endl;
    //(cov(xy).i()).print();
    //std::cout << "mean(xy):" << std::endl;
    //mean(xy).print();

    mat A = cov(xy).i()*mean(xy).t();

    std::cout << "A:" << std::endl;
    A.print();
}
