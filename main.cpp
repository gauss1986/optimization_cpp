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
#include <maxshp.h>
#include <comp_shp.h>
#include <ticktock.h>
#include "armadillo"

// Author Lin Gao, lingao.gao@mail.utoronto.ca
// Created March 5, 2017 for max sharpe problem at Bluewater Technologies Inc.

// The OLS routine is calling mkl_lapack package
// The statistics is calling BOOST package

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
double **matrix(int n,int m);
void free_matrix(double **a);

using namespace arma;
using namespace std;

int main(int argc, char *argv[])
{
    // initialize variables
    int N = 0; // No. of records
    int n = 2; // No. of vars. in linear model
    int m = 0; // No. of contracts    
    int N_bs = 500; // Bootstrap No. 
    double **newx;
    int *jpvt;

    // load and check files
    int N_x0_row,N_x0_col;
    int N_x_row,N_x_col;
    int N_y_row,N_y_col;
    string f_x0("x0.csv");
    string f_x("x.csv");
    string f_y("y.csv");
    vector< vector<double> >  x0 = readtxt(f_x0,N_x0_row,N_x0_col);
    vector< vector<double> >  x = readtxt(f_x,N_x_row,N_x_col);
    vector< vector<double> >  y = readtxt(f_y,N_y_row,N_y_col);
    // sanity check
    if ((N_x0_row != N_x_row) | (N_x_row != N_y_row)){
        cout << "The No. or records in x0, x and y are not consistent!\n";
    } 
    else{
        N = N_x0_row; 
    }
    if (N_x_col != N_y_col){
        cout << "The No. of columns in x and y are not consistent!\n";
    }
    else{
        m = N_x_col;
    }

    // compute and output the statistics of inputs
    //cout << "x0" <<  endl;
    vector< vector<double> > x0_col = reorgdata(x0,N_x0_row,N_x0_col);
    vector< vector<double> > x0_stat = simplestat(x0_col, 1);
    //cout << "x" <<  endl;
    vector< vector<double> > x_col = reorgdata(x,N_x_row,N_x_col);
    vector< vector<double> > x_stat = simplestat(x_col, n+1);
    //cout << "y" <<  endl;
    vector< vector<double> > y_col = reorgdata(y,N_y_row,N_y_col);
    vector< vector<double> > y_stat = simplestat(y_col, n+1);
    cout <<  endl;

    // covert x0,x,y to arma format
    mat mx(N,n+1);
    mat my(N,n+1);
    vec vx0(N);
    for (int i=0;i<N;i++){
        for (int j=0;j<n+1;j++){
            mx(i,j) = x[i][j];
            my(i,j) = y[i][j];
        }
        vx0(i) = x0[i][0];
    }

    // OLS
    cout << "OLS" << endl;    
    TickTock T;
    mat mA_OLS(N_bs,n+1);
    mat mshp_contract_OLS(N_bs,m);
    vec vshp(N_bs);
    // bootstrapping
    T.tick();
    for (int i=0;i<N_bs;i++){ 
        cout << double(i)/N_bs*100 << "%" << endl;
        vector<vector<double> > x0_sample;
        vector<vector<double> > x_sample;
        vector<vector<double> > y_sample;
        resample(x0_sample,x_sample,y_sample,x0,x,y);
        vec A_OLS = OLS(N, n, m, x0_sample, x_sample, y_sample);
        mA_OLS.row(i) = A_OLS.t();
        double shp_OLS;
        vec shp_contract_OLS = comp_shp(shp_OLS, N, m, n, A_OLS, mx, my, vx0);
        mshp_contract_OLS.row(i) = shp_contract_OLS.t();
        vshp(i) = shp_OLS;
    	cout << "Shp:" << shp_OLS <<  endl; 
	
    }
    //cout << "Shp:" << shp_OLS <<  endl; 
    //cout << "Shp per contract is:" <<  endl;
    //shp_contract_OLS.print();
    //cout << endl;
    T.tock("Bootstrapping costs ");
    vec tstat_A(n+1);
    for (int i=0;i<n+1;i++){
        tstat_A(i) = mean(mA_OLS.col(i))/stddev(mA_OLS.col(i));
    }
    cout << "Mean of the OLS parameters are:" << endl;
    mean(mA_OLS).print();
    cout << "T stats for the OLS parameters are:" << endl;
    tstat_A.print();
    double shp_OLS_mean;
    vec shp_contract_OLS_mean = comp_shp(shp_OLS_mean, N, m, n, mean(mA_OLS).t(), mx, my, vx0);
    cout << "Sharpe ratio is:" << shp_OLS_mean << endl;
    cout << "Shp per contract is:" <<  endl;
    shp_contract_OLS_mean.print();

    // maxsharpe
    cout << "Max sharpe:" << endl;    
    T.tick();
    vec A_MAX = maxshp(N, n, m, mx, my, vx0); 
    T.tock();
    cout << "coeff:" <<  endl;
    A_MAX.print();
    double shp_MAX;
    vec shp_contract_MAX = comp_shp(shp_MAX, N, m, n, A_MAX, mx, my, vx0);
    cout << "Shp:" << shp_MAX <<  endl; 
    cout << "Shp per contract is:" <<  endl;
    shp_contract_MAX.print();

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


