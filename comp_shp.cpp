#include <stdlib.h>
#include <stdio.h>
#include "armadillo"

using namespace std;
using namespace arma;

vec comp_shp(double& shp, const int N, const int m, const int n, vec A, mat& mx, mat& my, vec& vx0){
    // compute the sharpe ratio
    vec pnl(N);
    vec shp_contract(m); 
    for (int i=1;i<m;i++){
        mat completex(N,n+1);
        completex.col(0).ones();
        completex.col(1) = vx0;
        completex.col(2) = mx.col(i); 
        mat yfit = completex*A;
        vec pnl_contract = conv_to<vec>::from(my.col(i)%yfit);
        pnl = pnl + pnl_contract;
        shp_contract(i) = mean(pnl_contract)/stddev(pnl_contract);
    }
    shp = mean(pnl)/stddev(pnl);

    return shp_contract;
}
