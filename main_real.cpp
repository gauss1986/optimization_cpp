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
#include <iterator>
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

/* find end of a vector */
template<typename T, size_t N>
T * end(T (&ra)[N]){
	return ra+N;
}

/* select features from x and x0 */
mat selectx(vector<string>& x_select, vector<string>& x_pool, mat& raw, const int offset){
	vector<string>::iterator it1;
	vector<string>::iterator it;
	uvec x_ind(x_select.size(),fill::zeros);
	int i=0;
	for (it1=x_select.begin();it1!=x_select.end();++it1){
		for (it=x_pool.begin();it!=x_pool.end();++it){
			if (strcmp((*it1).c_str(),(*it).c_str())==0){
				x_ind(i) = distance(x_pool.begin(),it); 	
				cout << "selected dof " << x_ind(i) << endl;
				i++;
			}
		}
	}
	mat x = raw.cols(x_ind+offset);
	cout << "First record of selection is" << endl;
	x.row(0).print();

	return x;
}

int main(int argc, char *argv[])
{
	int Nx = 21;
	int Nx0 = 13;
	int N_bs = 500; // bootstrapping number
	const char *x_names[] = {"O.1","H.1","C.1"};
	vector<string> x_select(x_names,end(x_names));	
	const char *x0_names[] = {"SP.CC.1","TY.CC.1"};
	vector<string> x0_select(x0_names,end(x0_names));	

    // read data
    string f_x0("data_example_boot.txt");
	char sep = ' ';
	int N_row,N_col;
	vector<string> col_names; // column names
	vector<string> contracts_all; // contract names on each row
	vector<int> D;
    mat raw = readtxt_real(f_x0,N_row,N_col,sep,D,col_names,contracts_all);
	// print out some statistics of the data
	cout << "There are " << col_names.size()-3 << " number of variables." << endl;
	cout << "The variables are ";
	vector<string>::iterator it;
	for (it=col_names.begin()+3;it!=col_names.end();++it)
		cout << " " << *it;
	cout << endl;
	cout << "There are " << N_row << " records in total." << endl;
	vector<string> contracts = contracts_all;
	//contracts = contracts_all;
	it = unique(contracts.begin(),contracts.end());
	contracts.erase(it, contracts.end());
	cout << "There are in total " << contracts.size() << " contracts" << endl;
	cout << "The contracts are ";
	for (it=contracts.begin();it!=contracts.end();++it)
		cout << " " << *it;
	cout << endl;

	// parsing data
	ivec date = conv_to<ivec>::from(D);
	vec y = raw.col(0);
	mat x_all = raw.cols(1,Nx);
	mat x0_all = raw.cols(Nx+1,Nx+Nx0); 

	// select x from list of x_select 
	cout << "Selected x features: ";
	for (it=x_select.begin();it!=x_select.end();++it)
		cout << " " << *it;
	cout << endl;
	mat x = selectx(x_select, col_names, raw, -2); // offset -2 is needed to take care of the contract name and date
	int n = x.n_cols;
	cout << "Selected x0 features: ";
	for (it=x0_select.begin();it!=x0_select.end();++it)
		cout << " " << *it;
	cout << endl;
	mat x0 = selectx(x0_select, col_names, raw, -2);
	int n0 = x0.n_cols;
	cout << "n=" << n << ", n0=" << n0 << endl;
	
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

	// bootstrapping by day
	cout << "bootstrapping by day" << endl;
    TickTock T;
    T.tick();
	for (int i=0;i<N_bs;i++){
		if ((i%(N_bs/20)==0)) cout << double(i)/N_bs*100 << "%" << endl;
        // resampling
		uvec samplepoints = resample(N);
        vec c_OLS = solve(x_OLS.rows(samplepoints),y_OLS(samplepoints));
		vec c_MS = cov(xy_MS.rows(samplepoints)).i()*mean(xy_MS.rows(samplepoints)).t();
	}
    T.tock("Bootstrapping costs ");

	// just once
	vec c_OLS = solve(x_OLS,y_OLS);
	vec c_MS = cov(xy_MS).i()*mean(xy_MS).t();

	mat A(1+n+n0,2,fill::ones);
	A.col(1) = c_MS;
	vec c_norm = solve(A,c_OLS);
	vec c_MS_norm = c_MS*c_norm(1); 
	vec c_MS_norm2 = A*c_norm; 

	cout << endl;
	cout << "Results:" << endl;
	cout << "c_OLS	" << "c_MS	"<< "c_MS_norm	" << "c_MS_norm2"<< endl;
	cout << "intersection" << endl;
	cout << c_OLS(0) << "	" << c_MS(0) << "	" << c_MS_norm(0) << "	" <<  c_MS_norm2(0) << endl; 

	cout << "c_x" << endl;
	mat c_x = join_horiz(join_horiz(join_horiz(c_OLS(span(1,n)),c_MS(span(1,n))),c_MS_norm(span(1,n))),c_MS_norm2(span(1,n)));
	c_x.print();

	cout << "c_x0" << endl;
	mat c_x0 = join_horiz(join_horiz(join_horiz(c_OLS(span(n+1,n+n0)),c_MS(span(n+1,n+n0))),c_MS_norm(span(n+1,n+n0))),c_MS_norm2(span(n+1,n+n0)));
	c_x0.print();

	double shp_D_OLS;
	double shp_D_MS;
	vec shp_C_OLS = comp_shp_real(shp_D_OLS, c_OLS, x, x0, y, date, contracts_all, contracts, N);
	vec shp_C_MS = comp_shp_real(shp_D_MS, c_MS, x, x0, y, date, contracts_all, contracts, N);
	cout << "OLS shp=" << shp_D_OLS << endl;
	cout << "Contract shp: " << endl;
	shp_C_OLS.t().print();
	cout << "MS shp=" << shp_D_MS << endl;
	cout << "Contract shp: " << endl;
	shp_C_MS.t().print();
	

	return 1;
}
