#include <mkl_lapacke.h>
#include <stdlib.h>
#include <stdio.h>
#include "armadillo"
#include <maxshp.h>

using namespace arma;

vec maxshp(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0){
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

    return conv_to<vec>::from(A);
}
