#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <numeric>
#include <math.h>
#include <stdlib.h>
#include <readtxt.h>
#include <mkl_lapacke.h>
#include <gsl/gsl_cblas.h>

double **matrix(int n,int m);
void free_matrix(double **a);

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
    // re-organize data (swap col & row)
    //std::vector<std::vector<double> > x0_reog = reorgdata(x0, N_x0_row, N_x0_col);
    //std::vector<std::vector<double> > x_reog = reorgdata(x, N_x_row, N_x_col);
    //std::vector<std::vector<double> > y_reog = reorgdata(y, N_y_row, N_y_col);
    // sanity check
    if ((N_x0_col != N_x_col) | (N_x_col != N_y_col)){
        std::cout << "The No. or records in x0, x and y are not consistent!\n";
    } 
    else{
        N = N_x0_col; 
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
        newx[i][1] = x0[i][1]*m;
        newx[i][2] = std::accumulate(x[i].begin(),x[i].end(),0);
        newy[i] = std::accumulate(y[i].begin(),y[i].end(),0); 
    }

    // uses DGELSY to handle rank-deficient problems more realiably than DGELS.
    // Refer to http://www.netlib.org/lapack/lug/node71.html for details.
    jpvt = new int[n+1];
    for (int i=0;i<n+1;i++){
        jpvt[i] = 0;
    }
    int rank;
    LAPACKE_dgelsy(LAPACK_ROW_MAJOR, N, n+1, 1, &(newx[0][0]), n+1, newy, 1, jpvt, 1e-11, &rank);
    std::cout << "Results of OLS:" << newy[0] << "," << newy[1] << "," << newy[2] << "\n";     
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
