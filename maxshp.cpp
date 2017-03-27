#include <mkl_lapacke.h>
#include <stdlib.h>
#include <stdio.h>
#include "armadillo"
#include <maxshp.h>
#include <ticktock.h>

using namespace arma;
using namespace std;

vec maxshp(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0){
    // construct xy
    mat xy(N,n+1);
    xy.col(0) = sum(my,1);
    xy.col(1) = vx0%xy.col(0);
    xy.col(2) = sum(mx%my,1);
    
    // compute max sharpe setting
    //TickTock t;
    //t.tick();
    mat A = cov(xy).i()*mean(xy).t();
    //t.tock("Core of maxshp costs:");

    // compute stats

    return conv_to<vec>::from(A);
}
