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
#include <boost/math/distributions/students_t.hpp>

// Author Lin Gao, lingao.gao@mail.utoronto.ca
// Created March 5, 2017 for max sharpe problem at Bluewater Technologies Inc.

// The OLS routine is calling mkl_lapack package
// The statistics is calling BOOST package

using namespace arma;
using namespace std;
using namespace boost::math;

int main(int argc, char *argv[])
{
	int n = 21;
	int n0 = 13;

    // read data
    string f_x0("data_example_boot.txt");
	char sep = ' ';
	int N_row,N_col;
	vector<string> col_names; // column names
	vector<string> contracts; // contract names on each row
	vector<int> D;
    mat raw = readtxt_real(f_x0,N_row,N_col,sep,D,col_names,contracts);
	//cout << "N_row=" << N_row << ",N_col=" << N_col << endl;
	cout << "There are " << N_row << " records in total." << endl;

	// parsing data
	ivec date = conv_to<ivec>::from(D);
	vec y = raw.col(0);
	mat x = raw.cols(1,n);
	mat x0 = raw.cols(n+1,n+n0); 
	cout << "Size of x=21," << "x0=13" << endl;

	// find unique dates
	uvec ind = find_unique(date);
	int N = ind.n_elem;
	cout << "There are "<< N << " unique days starting " << min(date(ind)) << " ending " << max(date(ind)) << endl;;

	// parsing data w.r.t. dates
	uvec c_D(N,fill::zeros); // contracts on record per day
	mat x_OLS(N,1+n+n0);
	vec y_OLS(N);
	mat xy_MS(N,1+n+n0);
	for (int i=0;i<N;i++){
		// index for records on the same day
		uvec m_i = find(date==date(ind(i)));	
		c_D(i) = m_i.n_elem;
		// day x, x0 and y
		mat x_D = x.rows(m_i);
		mat x0_D = x0.rows(m_i);
		vec y_D = y(m_i);
		// sanity check x0 values on the same day, should be identical for different contracts
		if (as_scalar(max(max(x0_D)-min(x0_D),1)) > 1e-5) 
			cout << "x0 values on day " << date(ind(i)) << " are not identical for different contracts!" << endl;	
		// set y_OLS
		y_OLS(i) = sum(y_D);
		// set x_OLS
		x_OLS(i,0) = c_D(i);
		x_OLS(i,span(1,n)) = sum(x_D);
		x_OLS(i,span(n+1,n+n0)) = sum(x0.rows(m_i));
		// set xy_MS
		for (int j=0;j<n;j++){
			xy_MS(i,j+1) = sum(x_D.col(j)%y_D);
		}
		for (int j=0;j<n0;j++){
			xy_MS(i,1+n+j) = x0_D(0,j)*y_OLS(i);
		}
	}
	xy_MS.col(0) = y_OLS;
	cout << "Min contracts per day " << min(c_D) << ", max " << max(c_D) << ", mean " << mean(c_D) << endl;

	vec c_OLS = solve(x_OLS,y_OLS);
	
	vec c_MS = cov(xy_MS).i()*mean(xy_MS).t();

	cout << "c_OLS=" << endl;
	c_OLS.t().print();
	cout << "c_MS=" << endl;
	c_MS.t().print();

	double shp_D_OLS;
	double shp_D_MS;
	comp_shp_real(shp_D_OLS, c_OLS, x, x0, y, date, ind, N);
	comp_shp_real(shp_D_MS, c_MS, x, x0, y, date, ind, N);
	cout << "OLS shp=" << shp_D_OLS << endl;
	cout << "MS shp=" << shp_D_MS << endl;

	return 1;
}
