#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <numeric>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <txtIO.h>
#include <mkl_lapacke.h>
#include <gsl/gsl_cblas.h>

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
double **matrix(int n,int m);
void free_matrix(double **a);
template <class T> void printdata_tmp(T& data, int N_row, int N_col);

int main(int argc, char *argv[])
{
    // initialize variables
    int N = 0; // No. of records
    int n = 2; // No. of vars. in linear model
    int m = 0; // No. of contracts    
    double **newx;
    double *newy; 
    int *jpvt;

    // load and check files
    int N_x0_row,N_x0_col;
    int N_x_row,N_x_col;
    int N_y_row,N_y_col;
    std::string f_x0("x0.csv");
    std::string f_x("x.csv");
    std::string f_y("y.csv");
    std::vector<std::vector<double> >  x0 = readtxt(f_x0,N_x0_row,N_x0_col);
    std::vector<std::vector<double> >  x = readtxt(f_x,N_x_row,N_x_col);
    std::vector<std::vector<double> >  y = readtxt(f_y,N_y_row,N_y_col);
    printdata_tmp<std::vector<std::vector<double> > >(x0,10,1);
    // output size of data
    //std:: cout<< "Size of x0 is "<<N_x0_row<<" rows, "<<N_x0_col<< " cols.\n";
    //std:: cout<< "Size of x is "<<N_x_row<<" rows, "<<N_x_col<< " cols.\n";
    //std:: cout<< "Size of y is "<<N_y_row<<" rows, "<<N_y_col<< " cols.\n";
    // re-organize data (swap col & row)
    //std::vector<std::vector<double> > x0_reog = reorgdata(x0, N_x0_row, N_x0_col);
    //std::vector<std::vector<double> > x_reog = reorgdata(x, N_x_row, N_x_col);
    //std::vector<std::vector<double> > y_reog = reorgdata(y, N_y_row, N_y_col);
    // sanity check
    if ((N_x0_row != N_x_row) | (N_x_row != N_y_row)){
        std::cout << "The No. or records in x0, x and y are not consistent!\n";
    } 
    else{
        N = N_x0_row; 
    }
    if (N_x_col != N_y_col){
        std::cout << "The No. of columns in x and y are not consistent!\n";
    }
    else{
        m = N_x_col;
    }

    // construct newx, newy
    newx = matrix(N,n+1);
    newy = new double[N];
    for (int i=0;i<N;i++){
        newx[i][0] = m;
        newx[i][1] = x0[i][0]*m;
        newx[i][2] = std::accumulate(x[i].begin(),x[i].end(),0);
        newy[i] = std::accumulate(y[i].begin(),y[i].end(),0); 
        //std::cout << "newx=" << newx[i][0] << "," << newx[i][1] << "," << newx[i][2] << ". newy=" << newy[i] << "." << std::endl;
    }
    //printdata(newx,10,3);
    //printdata(newy,10,1);

    // uses DGELSY to handle rank-deficient problems more realiably than DGELS.
    // Refer to http://www.netlib.org/lapack/lug/node71.html for details.
    jpvt = new int[n+1];
    for (int i=0;i<n+1;i++){
        jpvt[i] = 0;
    }
    int rank;
    //LAPACKE_dgelsy(LAPACK_ROW_MAJOR, N, n+1, 1, &(newx[0][0]), n+1, &(newy[0]), 1, &(jpvt[0]), 1e-11, &rank);
    LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', N, n+1, 1, &(newx[0][0]), n+1, &(newy[0]), 1);
    /* Print least squares solution */
    print_matrix( "Least squares solution", n+1, 1, newy, 1 );
    /* Print residual sum of squares for the solution */
    print_vector_norm( "Residual sum of squares for the solution", N-n-1, 1, &newy[n+1], 1 );
    /* Print details of QR factorization */
    //print_matrix( "Details of QR factorization", N, n+1, &(newx[0][0]), n+1 );

    // free matrix and array
    free_matrix(newx);
    delete[] newy;
    delete[] jpvt;
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

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing norms of matrix columns */
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
//void printdata(const std::vector<std::vector<double> >& data, int N_row, int N_col){
template <class T> void printdata_tmp(T& data, int N_row, int N_col){
    // output content 
    std::cout << "Printing the first " << N_row << " rows of the data" << std::endl;
    for (int i=0;i<N_row;i++){    
        for (int j=0;j<N_col;j++){
            std::cout << data[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

