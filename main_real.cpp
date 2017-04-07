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
	mat x = raw.cols(1,21);
	mat x0 = raw.cols(22,34); 
	cout << "Size of x=21," << "x0=13" << endl;

	uvec ind = find_unique(date);
	cout << "There are "<< ind.n_elem << " unique days starting " << min(date(ind)) << " ending " << max(date(ind)) << endl;;

	return 1;
}
