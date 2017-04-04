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

void comp_stat(vec vy, vec vc, mat mx){
    // compute and print stats
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
    std::cout << "Estimat   " << "Std. Error    " << "t value   " << "Pr(>|t|)" << std::endl;
    mat results(mx.n_cols,4);
    results.col(0) = vc;
    results.col(1) = se;
    results.col(2) = t;
    results.col(3) = q;
    results.print();
    std::cout << "Residual standard error: " << sqrt(sigmasq) << " on " << N-mx.n_cols << " degrees of freedom"<< std::endl;
    std::cout << "Multiple R-squared: " << R2 << std::endl;
    std::cout << "Adjusted R-squared: " << R2_adj << std::endl;
}

vec OLS_day(const int N, const int n, const int m, std::vector<std::vector<double> >& x0, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y){
    // construct newx, newy
    double **newx = matrix(N,n+1);
    double *newy = new double[N];
    std::vector<double> newy_copy;
    std::vector<std::vector<double> > newx_copy;
    for (int i=0;i<N;i++){
        newx[i][0] = m;
        newx[i][1] = x0[i][0]*m;
        newx[i][2] = std::accumulate(x[i].begin(),x[i].end(),0.0);
        newy[i] = std::accumulate(y[i].begin(),y[i].end(),0.0);
    }
    for (int i=0;i<N;i++){
        newy_copy.push_back(newy[i]);
        std::vector<double> newx_copy_1D;
        newx_copy_1D.push_back(newx[i][0]);
        newx_copy_1D.push_back(newx[i][1]);
        newx_copy_1D.push_back(newx[i][2]);
        newx_copy.push_back(newx_copy_1D);
    }

    // DGELS is general purpose and most efficient.
    // Use DGELSY to handle rank-deficient problems more realiably than DGELS.
    // Refer to http://www.netlib.org/lapack/lug/node71.html for details.
    //TickTock T;
    //T.tick();
    LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', N, n+1, 1, &(newx[0][0]), n+1, &(newy[0]), 1);
    //T.tock("DGELS costs ");

    // Compute and print the stats on OLS
    mat mx;
    vec vc = vec(newy,n+1);
    vec vy = conv_to<vec>::from(newy_copy);
    vec2D_to_arma(newx_copy, mx); 
    //comp_stat(vy, vc, mx);
    //OLS_stat1.print();

    /* Print residual sum of squares for the solution */
    //print_vector_norm( "Residual sum of squares for the solution", N-n-1, 1, &newy[n+1], 1 );
    /* Print details of QR factorization */
    //print_matrix( "Details of QR factorization", N, n+1, &(newx[0][0]), n+1 );

    // free matrix and array
    free_matrix(newx);
    vec vcoeff = vec(newy,n+1);
    return vcoeff;
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

    vec A = solve(newx,newy);
    cout << "The coeffs are " << endl;
    A.t().print(); 

    // Compute and print the stats on OLS
    comp_stat(newy,A,newx);

    return A;
}

