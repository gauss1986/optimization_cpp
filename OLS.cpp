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
#include "armadillo"
#include <boost/math/distributions/students_t.hpp>

using namespace arma;
using namespace boost::math;

//class OLS_stat{
    // Compute OLS stats
    // Reference 1: https://en.wikipedia.org/wiki/Simple_linear_regression#Normality_assumption 
    // Reference 2: http://stats.stackexchange.com/questions/44838/how-are-the-standard-errors-of-coefficients-calculated-in-a-regression/44841#44841
    // Reference 3: http://www.netlib.org/lapack/lawnspdf/lawn193.pdf
    // Reference 4: http://qed.econ.queensu.ca/pub/faculty/abbott/econ351/351note04.pdf
//    double sigmasq; // MSE
//    double R2; // R-square
//    double R2_adj; // adjusted R-square
//    int N; // No. of sample points in linear fit
//    mat mx_;// arma format of x
//    vec vy_;// arma format of y
//    vec vc_;// arma format of coefficient
//    vec se;// arma format of standard error
//    vec t; // arma format of t stats
//    vec q; // arma format of probability that value with t stats is due to chance
 //   public:
        // convert x,y and coeff to arma format to prepare for stats calculation
        //void conv_form(const int n, const std::vector<std::vector<double> >& x, const std::vector<double>& y, const double* coeff);
        // compute stats of the coeffs
//        void comp_stat(vec vy, vec vc, mat mx);
        // print out the stats
//        void print();
//};

void convert_to_arma(const int n, const std::vector<std::vector<double> >& x, const std::vector<double>& y, const double* coeff, int N, vec& vy, vec& vc, mat& mx){
    N = x.size();
    // convert vectors and arrays to arma mat and vec 
    vy = conv_to<vec>::from(y);
    mx.set_size(N,n+1);
    for (int i=0;i<N;i++){
        for (int j=0;j<n+1;j++){
            mx(i,j) = x[i][j];
        }
    }
    vc = vec(coeff,n+1);
}

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
    vec vy;
    vec vc;
    mat mx;
    convert_to_arma(n, newx_copy, newy_copy, newy, N, vy, vc, mx); 
    comp_stat(vy, vc, mx);
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
    //OLS_stat OLS_stat1;
    //OLS_stat1.conv_form(n, newx_copy, newy_copy, newy); 
    comp_stat(newy,A,newx);
    //OLS_stat1.print();

    /* Print residual sum of squares for the solution */
    //print_vector_norm( "Residual sum of squares for the solution", N-n-1, 1, &newy[n+1], 1 );
    /* Print details of QR factorization */
    //print_matrix( "Details of QR factorization", N, n+1, &(newx[0][0]), n+1 );

    // free matrix and array
    return A;
}

double **matrix(int n,int m) {
    double **a = new double * [n];
    a[0] = new double [n*m];
    
    for (int i=1; i<n; i++)
        a[i] = &a[0][i*m];
    
    return a;
}

void free_matrix(double **a) {
    delete[] a[0];
    delete[] a;
}

/* Auxiliary routine from Intel MKL example: printing norms of matrix columns */
void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
    MKL_INT i, j;
    double norm;
    printf( "\n %s\n", desc );
    for( j = 0; j < n; j++ ) {
        norm = 0.0;
        for( i = 0; i < m; i++ ) norm += a[i*lda+j] * a[i*lda+j];
        printf( " %6.2f", norm );
    }
    printf( "\n" );
}
