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
	mat newx(N,n+n0);
	vec newy(N);
	for (int i=0;i<N;i++){
		uvec m_i = find(date==date(ind(i)));	
		c_D(i) = m_i.n_elem;
		// set value for newx & newy
		newy(i) = c_D(i);
	    vec temp(n,fill::zeros);
		for (int j=0;j<n;j++) 	
	}
	cout << "Min contracts per day " << min(c_D) << ", max " << max(c_D) << ", mean " << mean(c_D) << endl;;

	return 1;
}
