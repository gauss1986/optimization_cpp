#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <stdio.h>
#include <math.h>
#include <txtIO.h>
#include <stat.h>
#include <algorithm>
#include <gsl/gsl_cblas.h>
#include <OLS.h>

// Author Lin Gao, lingao.gao@mail.utoronto.ca
// Created March 5, 2017 for max sharpe problem at Bluewater Technologies Inc.

// The OLS routine is calling mkl_lapack package
// The statistics is calling BOOST package

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
double **matrix(int n,int m);
void free_matrix(double **a);

int main(int argc, char *argv[])
{
    // initialize variables
    int N = 0; // No. of records
    int n = 2; // No. of vars. in linear model
    int m = 0; // No. of contracts    
    double **newx;
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

    // compute and output the statistics of inputs
    //std::cout << "x0" << std::endl;
    std::vector<std::vector<double> > x0_col = reorgdata(x0,N_x0_row,N_x0_col);
    std::vector<std::vector<double> > x0_stat = simplestat(x0_col, 1);
    //std::cout << "x" << std::endl;
    std::vector<std::vector<double> > x_col = reorgdata(x,N_x_row,N_x_col);
    std::vector<std::vector<double> > x_stat = simplestat(x_col, n+1);
    //std::cout << "y" << std::endl;
    std::vector<std::vector<double> > y_col = reorgdata(y,N_y_row,N_y_col);
    std::vector<std::vector<double> > y_stat = simplestat(y_col, n+1);
    std::cout << std::endl;

    // OLS
    std::cout << "OLS" << std::endl;
    double *newy = OLS(N, n, m, x0, x, y);

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


