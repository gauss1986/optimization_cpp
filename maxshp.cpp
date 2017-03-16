#include <mkl_lapacke.h>
#include <stdlib.h>
#include <stdio.h>
#include "armadillo"
#include <maxshp.h>

using namespace arma;
using namespace std;

vec maxshp(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0){
    // construct xy
    mat xy(N,n+1);
    xy.col(0) = sum(my,1);
    xy.col(1) = vx0%xy.col(0);
    xy.col(2) = sum(mx%my,1);
    
    // print out some details if wanted
    //cout << "Cov(xy):" << endl;
    //cov(xy).print();
    //cout << "Cov(xy)^(-1):" << endl;
    //(cov(xy).i()).print();
    //cout << "mean(xy):" << endl;
    //mean(xy).print();

    // compute max sharpe setting
    mat A = cov(xy).i()*mean(xy).t();

    return conv_to<vec>::from(A);
}
