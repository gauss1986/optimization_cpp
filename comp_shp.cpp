#include <stdlib.h>
#include <stdio.h>
#include "armadillo"

using namespace std;
using namespace arma;

vec comp_shp(double& shp, const int N, const int m, const int n, vec A, mat& mx, mat& my, vec& vx0){
    // compute the sharpe ratio
    vec pnl = zeros(N);
    vec shp_contract(m); 
    for (int i=0;i<m;i++){
        mat completex(N,n+1);
        completex.col(0).ones();
        completex.col(1) = vx0;
        completex.col(2) = mx.col(i); 
        vec yfit = conv_to<vec>::from(completex*A);
        vec pnl_contract = my.col(i)%yfit;
        pnl = pnl + pnl_contract;
        shp_contract(i) = mean(pnl_contract)/stddev(pnl_contract);
    }
    shp = mean(pnl)/stddev(pnl);
    std::cout << "Shp:" << shp << std::endl; 

    return shp_contract;
}
