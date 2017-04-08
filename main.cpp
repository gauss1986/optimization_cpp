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
#include <OLS.h>
#include <maxshp.h>
#include <misc.h>
#include <ticktock.h>
#include "armadillo"

// Author Lin Gao, lingao.gao@mail.utoronto.ca
// Created March 5, 2017 for max sharpe problem at Bluewater Technologies Inc.

// The OLS routine is calling mkl_lapack package
// The statistics is calling BOOST package


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
    int N_x0_row,N_x0_col,N_x_row,N_x_col,N_y_row,N_y_col;
    string f_x0("x0.csv");
    string f_x("x.csv");
    string f_y("y.csv");

	char sep = ',';
    cout << "x0" <<  endl;
    mat mx0 = readtxt(f_x0,N_x0_row,N_x0_col,sep);
	vec vx0 = conv_to<vec>::from(mx0);
    cout << "x" <<  endl;
    mat mx = readtxt(f_x,N_x_row,N_x_col,sep);
    cout << "y" <<  endl;
    mat my = readtxt(f_y,N_y_row,N_y_col,sep);
    cout <<  endl;

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

    TickTock T;
    mat mA_OLS(N_bs,n+1);
    mat mshp_contract_OLS(N_bs,m);
    vec vshp_OLS(N_bs);
    mat mA_MAX(N_bs,n+1);
    mat mshp_contract_MAX(N_bs,m);
    vec vshp_MAX(N_bs);
    mat mA_MAX_norm(N_bs,n+1);
    mat mshp_contract_MAX_norm(N_bs,m);
    vec vshp_MAX_norm(N_bs);
    cout << "Bootstrapping size " << N_bs << endl;
    // bootstrapping
    T.tick();
    for (int i=0;i<N_bs;i++){ 
	if ((i%(N_bs/100)==0)) cout << double(i)/N_bs*100 << "%" << endl;
        // resampling
        vec vx0_sample(N);
        mat mx_sample(N,m);
        mat my_sample(N,m);
        resample_arma(vx0_sample,mx_sample,my_sample,vx0,mx,my);
        // OLS
        vec A_OLS = OLS_day(N, n, m, mx_sample, my_sample, vx0_sample);
        mA_OLS.row(i) = A_OLS.t();
        double shp_OLS;
        vec shp_contract_OLS = comp_shp(shp_OLS, n, A_OLS, mx, my, vx0);
        mshp_contract_OLS.row(i) = shp_contract_OLS.t();
        vshp_OLS(i) = shp_OLS;
        // maxsharpe
        vec A_MAX = maxshp(N, n, m, mx_sample, my_sample, vx0_sample);
        mA_MAX.row(i) = A_MAX.t(); 
        double shp_MAX;
        vec shp_contract_MAX = comp_shp(shp_MAX, n, A_MAX, mx, my, vx0);
        mshp_contract_MAX.row(i) = shp_contract_MAX.t();
        vshp_MAX(i) = shp_MAX;
		// normalize maxsharpe coeffs w.r.t. OLS coeffs.
		mat A(n+1,2,fill::ones);
		A.col(1) = A_MAX;
		vec c = solve(A,A_OLS);
		c(0) = 0;
		vec A_MAX_norm = A*c; 
        mA_MAX_norm.row(i) = A_MAX_norm.t(); 
        double shp_MAX_norm;
        vec shp_contract_MAX_norm = comp_shp(shp_MAX_norm, n, A_MAX, mx, my, vx0);
        mshp_contract_MAX_norm.row(i) = shp_contract_MAX_norm.t();
        vshp_MAX_norm(i) = shp_MAX_norm;
    }
    T.tock("Bootstrapping costs ");
    cout << "Statistics: "<< endl;

    cout << "OLS" << endl;
    bsstat(mA_OLS, vx0, mx, my, vshp_OLS, mshp_contract_OLS, N, n, m);
    cout << "MAX Sharpe" << endl;
    bsstat(mA_MAX, vx0, mx, my, vshp_MAX, mshp_contract_MAX, N, n, m);
    cout << "MAX Shapre normalized" << endl;
    bsstat(mA_MAX_norm, vx0, mx, my, vshp_MAX_norm, mshp_contract_MAX_norm, N, n, m);

    // OLS_by_day just once
    cout << "OLS just once" << endl;
    vec A_OLS_1 = OLS_day(N, n, m, mx, my, vx0);
    double shp_OLS_1;
    vec shp_contract_OLS_1 = comp_shp(shp_OLS_1, n, A_OLS_1, mx, my, vx0);
    cout << "Coeffs are " << endl;
    A_OLS_1.t().print();
    cout << "Portfolio Shp " << shp_OLS_1 << endl;
    cout << "Contract Shp " << endl;
    shp_contract_OLS_1.t().print();
    cout << endl;

    // OLS_by_record just once
    cout << "Naive OLS" << endl;
    vec A_OLS_record = OLS_record(N,n,m,mx,my,vx0); 
    double shp_OLS_2;
    vec shp_contract_OLS_2 = comp_shp(shp_OLS_2, n, A_OLS_record, mx, my, vx0);
    cout << "Coeffs are " << endl;
    A_OLS_record.t().print();
    cout << "Portfolio Shp " << shp_OLS_2 << endl;
    cout << "Contract Shp " << endl;
    shp_contract_OLS_2.t().print();
    cout << endl;

    // maxsharpe just once
    cout << "Max Sharpe just once" << endl;
    vec A_MAX_1 = maxshp(N, n, m, mx, my, vx0);
    mat A1(m,2,fill::ones); // normalize
    A1.col(1) = A_MAX_1;
    vec c1 = solve(A1,A_OLS_1);
    cout << "Regression coeffs are " << endl;
    //c1.print();
    c1(0) = 0;
    vec A_MAX_1norm = A1*c1; 
    //cout << "Raw coeffs are " << endl;
    A_MAX_1.t().print();
    double shp_MAX_1;
    vec shp_contract_MAX_1 = comp_shp(shp_MAX_1, n, A_MAX_1, mx, my, vx0);
    //cout << "Portfolio Shp " << shp_MAX_1 << endl;
    //cout << "Contract Shp " << endl;
    //shp_contract_MAX_1.t().print();
    cout << "Normalized coeffs are " << endl;
    A_MAX_1norm.t().print();
    shp_contract_MAX_1 = comp_shp(shp_MAX_1, n, A_MAX_1norm, mx, my, vx0);
    cout << "Portfolio Shp " << shp_MAX_1 << endl;
    cout << "Contract Shp " << endl;
    shp_contract_MAX_1.t().print();
}

