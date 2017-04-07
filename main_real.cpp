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

/* Auxiliary routines prototypes */
void linearstat(mat mA, vec vx0, mat mx, mat my, vec vshp, mat mshp_contract, const int N, const int n, const int m);

using namespace arma;
using namespace std;
using namespace boost::math;

int main(int argc, char *argv[])
{
    // initialize variables

    string f_x0("data_example_boot.txt");

	char sep = ' ';
	int N_row,N_col;
    mat raw = readtxt(f_x0,N_row,N_col,sep);

	

	return 1;
}
