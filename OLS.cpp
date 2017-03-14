#include <numeric>
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
#include "armadillo"

class OLS_stat{
    // Compute OLS stats
    // Reference 1: https://en.wikipedia.org/wiki/Simple_linear_regression#Normality_assumption 
    // Reference 2: http://stats.stackexchange.com/questions/44838/how-are-the-standard-errors-of-coefficients-calculated-in-a-regression/44841#44841
    // Reference 3: http://www.netlib.org/lapack/lawnspdf/lawn193.pdf
    // Reference 4: http://qed.econ.queensu.ca/pub/faculty/abbott/econ351/351note04.pdf
    double sigmasq;
    public:
        void compute_sigmasq(const int n, const std::vector<std::vector<double> >& x, const std::vector<double>& y, const double* coeff);
        void print();
};

using namespace arma;

void OLS_stat::compute_sigmasq(const int n, const std::vector<std::vector<double> >& x, const std::vector<double>& y, const double* coeff){
    //vec(y)-mat(x)*vec(coeff) 
    vec vy = conv_to<vec>::from(y);
    std::cout << "vector y is:" << std::endl;
    vy.print();
    mat mx(x.size(),n+1);
    std::cout << "matrix x is of size " << x.size() << "x" << n+1 << std::endl;
    for (int i=0;i<x.size();i++){
        for (int j=0;j<n+1;j++){
            mx(i,j) = x[i][j];
            std::cout << x[i][j] << " ";
        }
        std::cout << std::endl;
    }
    mx.print();
    std::cout << "vector coeff is:" << std::endl;
    vec vc(x[0].size());
    for (int i=0;i<n+1;i++){
        vc(i) = coeff[i];
    }
    vc.print();
    vy = vy-mx*vc;
    vy = vy%vy;
    sigmasq = sum(vy)/(x.size()-n); 
}

void OLS_stat::print(){
    std::cout << "Sigma sq=" << sigmasq << std::endl;
}

double *OLS(const int N, const int n, const int m, std::vector<std::vector<double> >& x0, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y){
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
        newy_copy.push_back(newy[i]);
        std::vector<double> newx_copy_1D;
        newx_copy_1D.push_back(newx[i][0]);
        newx_copy_1D.push_back(newx[i][1]);
        newx_copy_1D.push_back(newx[i][2]);
        newx_copy.push_back(newx_copy_1D);
        //std::cout << "newx=" << newx[i][0] << "," << newx[i][1] << "," << newx[i][2] << ". newy=" << newy[i] << "." << std::endl;
    }
    printdata_2D<double **>(newx,"First 10 rows in newx are:",10,3);
    printdata_1D<double *>(newy,"First 10 entries in newy are:",10);

    // DGELS is general purpose and most efficient.
    // Use DGELSY to handle rank-deficient problems more realiably than DGELS.
    // Refer to http://www.netlib.org/lapack/lug/node71.html for details.
    LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', N, n+1, 1, &(newx[0][0]), n+1, &(newy[0]), 1);

    /* Print least squares solution */
    printdata_1D<double *>(newy,"Least squares solution",n+1);

    OLS_stat OLS_stat1;
    OLS_stat1.compute_sigmasq(n, newx_copy, newy_copy, newy); 
    OLS_stat1.print();

    /* Print residual sum of squares for the solution */
    //print_vector_norm( "Residual sum of squares for the solution", N-n-1, 1, &newy[n+1], 1 );
    /* Print details of QR factorization */
    //print_matrix( "Details of QR factorization", N, n+1, &(newx[0][0]), n+1 );

    // free matrix and array
    free_matrix(newx);
    return newy;
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
