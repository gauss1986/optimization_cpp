#include <numeric>
#include <math.h>
#include <sstream>
#include <vector>
#include <iterator>
#include <fstream>
#include <iostream>
#include <mkl_lapacke.h>
#include <txtIO.h>
#include <stdlib.h>
#include <stdio.h>
#include <OLS.h>
#include <ticktock.h>
#include <misc.h>
#include "armadillo"
#include <boost/math/distributions/students_t.hpp>

using namespace arma;
using namespace boost::math;
using namespace std;

/* OLS by day */
vec OLS_day(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0){
	// generate newx, newy
	vec newy = sum(my,1);
	mat newx(N,n+1);
	newx.col(0).fill(m);
	newx.col(1) = vx0*m;
	newx.col(2) = sum(mx,1);

	vec vc = solve(newx,newy);

    // Compute and print the stats on OLS
    //comp_stat(newy,vc,newx);

    return vc;
}

/* OLS by record */
vec OLS_record(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0){
	// construct newx, newy
    mat newx(3*N,n+1,fill::zeros);
    vec newy(3*N,fill::zeros);
    for (int i=0;i<m;i++){
    	mat newx_block(N,n+1,fill::ones);
    	newx_block.col(1) = vx0;
    	newx_block.col(2) = mx.col(i);
        newx.rows(i*N,(i+1)*N-1) = newx_block;
		newy.rows(i*N,(i+1)*N-1) = my.col(i);
    }

    vec vc = solve(newx,newy);

    // Compute and print the stats on OLS
    //comp_stat(newy,vc,newx);

    return vc;
}

/* compute and print stats */
void comp_stat(vec vy, vec vc, mat mx){
    int N = vy.n_elem;
    // sigmasq
    vec vy_copy(vy);
    vy = vy-mx*vc;
    vy = vy%vy;
    double sigmasq = sum(vy)/(N-mx.n_cols); 
    // standard error
    vec se(mx.n_cols);
    se = sqrt(diagvec(sigmasq*(mx.t()*mx).i())); 
    // t stats
    vec t(mx.n_cols);
    t = vc/se; 
    students_t dist(N-1);
    vec q(mx.n_cols);
    for (int i=0;i<mx.n_cols;i++){
        q(i) = 2*cdf(complement(dist, fabs(t(i))));
    }
    // R2 stats
    double R2 = 1-sum(vy)/var(vy_copy)/(N-1);
    // adjusted R2 stats
    double R2_adj = R2-(1-R2)*(mx.n_cols-1)/(N-mx.n_cols);
    // print out the results
    cout << "Estimat   " << "Std. Error    " << "t value   " << "Pr(>|t|)" << endl;
    mat results(mx.n_cols,4);
    results.col(0) = vc;
    results.col(1) = se;
    results.col(2) = t;
    results.col(3) = q;
    results.print();
    cout << "Residual standard error: " << sqrt(sigmasq) << " on " << N-mx.n_cols << " degrees of freedom"<< endl;
    cout << "Multiple R-squared: " << R2 << endl;
    cout << "Adjusted R-squared: " << R2_adj << endl;
}
