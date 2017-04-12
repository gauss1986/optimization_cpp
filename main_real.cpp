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

/* print out stats of bootstrapping */
void printbs(const vec& m_OLS, const vec& s_OLS, const vec& t_OLS, const vector<string>& x_select, const vector<string>& x0_select, int n, int n0){
	cout << "		" << "Estimate	" << "Std. Error	" << "t value	" << endl;  
	cout << "(Intercept)	";
	vec output;
	output <<  m_OLS(0) << s_OLS(0) << t_OLS(0); 
	output.t().print();
	for (int i=0;i<n;i++){
		cout << x_select[i] << "		";
		vec output;
		output <<  m_OLS(i+1) << s_OLS(i+1) << t_OLS(i+1); 
		output.t().print();
	}
	for (int i=0;i<n0;i++){
		//cout << x0_select[i] << "		"<< m_OLS(i+1+n)<< "	" << s_OLS(i+1+n) << "		" << t_OLS(i+1+n) << endl; 
		cout << x0_select[i] << "		";
		vec output;
		output <<  m_OLS(i+n+1) << s_OLS(i+n+1) << t_OLS(i+n+1); 
		output.t().print();
	}
}

int main(int argc, char *argv[])
{
	int Nx = 21;
	int Nx0 = 13;
	int N_bs = 5; // bootstrapping number
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
		//if (as_scalar(max(max(x0_D)-min(x0_D),1)) > 1e-5) 
		//	cout << "x0 values on day " << date(ind(i)) << " are not identical for different contracts!" << endl;	
		// set y_OLS
		xy_MS(i,0) = sum(y_D);
		// set x_OLS
		//x_OLS(i,0) = c_D(i);
		//x_OLS(i,span(1,n)) = sum(x_D);
		//x_OLS(i,span(n+1,n+n0)) = sum(x0.rows(m_i));
		// set xy_MS
		for (int j=0;j<n;j++){
			xy_MS(i,j+1) = sum(x_D.col(j)%y_D);
		}
		for (int j=0;j<n0;j++){
			xy_MS(i,1+n+j) = x0_D(0,j)*sum(y_D);
		}
	}
	//xy_MS.col(0) = y_OLS;
	cout << "Min contracts per day " << min(c_D) << ", max " << max(c_D) << ", mean " << mean(c_D) << endl;
	cout << "sum(c_D)=" << sum(c_D) << endl;

	// bootstrapping by day
	cout << "bootstrapping by day" << endl;
    TickTock T;
	mat mc_OLS(N_bs,1+n+n0);
	mat mc_MS(N_bs,1+n+n0);
	for (int i=0;i<N_bs;i++){
		//if ((i%(N_bs/100)==0)) cout << double(i)/N_bs*100 << "%" << endl;
        cout << "N_bs=" << i << endl;
		// resampling
		uvec samplepoints = resample(N);
		int N_record = sum(c_D(samplepoints));
		//cout << "N_record=" << N_record << endl;
		mat x_OLS(N_record,1+n+n0,fill::ones);
		vec y_OLS(N_record);
		int j = 0;
		for (int k=0;k<N;k++){
			//cout << "Day " << k << "/" << N ;
			// index for records on the same day
			uvec m_i = find(date==date(ind(samplepoints(k))));	
			//cout << "C_D(k)=" << c_D(samplepoints(k)) << "vs m_i.n_elem=" << m_i.n_elem << endl;
			x_OLS(span(j,j+m_i.n_elem-1),span(1,n)) = x.rows(m_i);
			x_OLS(span(j,j+m_i.n_elem-1),span(n+1,n+n0)) = x0.rows(m_i);
			y_OLS(span(j,j+m_i.n_elem-1)) = y(m_i);
			j = j+m_i.n_elem;
			//cout << ", Record " << j << "/" << N_record << endl;
		}
		//if (j!=N_row) cout << "Bootstrapping has some issues!" << endl;
		// OLS
		T.tick();
        vec c_OLS = solve(x_OLS,y_OLS);
    	T.tock("OLS costs ");
		mc_OLS.row(i) = c_OLS.t();
		// MS
		T.tick();
		vec c_MS = cov(xy_MS.rows(samplepoints)).i()*mean(xy_MS.rows(samplepoints)).t();
    	T.tock("MS costs ");
		mc_MS.row(i) = c_MS.t();
	}
	vec m_OLS = arma::mean(mc_OLS).t();
	vec m_MS = arma::mean(mc_MS).t();
	vec s_OLS = stddev(mc_OLS).t();
	vec s_MS = stddev(mc_MS).t();
    vec t_OLS = m_OLS/s_OLS;
    vec t_MS = m_MS/s_MS;
	cout << "Stats of OLS:" << endl;
	printbs(m_OLS, s_OLS, t_OLS, x_select, x0_select, n, n0);
	cout << "Stats of MS:" << endl;
	printbs(m_MS, s_MS, t_MS, x_select, x0_select, n, n0);

	double shp_D_OLS;
	double shp_D_MS;
	vec shp_C_OLS = comp_shp_real(shp_D_OLS, m_OLS, x, x0, y, date, contracts_all, contracts, N);
	vec shp_C_MS = comp_shp_real(shp_D_MS, m_MS, x, x0, y, date, contracts_all, contracts, N);
	cout << "OLS shp=" << shp_D_OLS << endl;
	cout << "Contract shp: " << endl;
	//shp_C_OLS.t().print();
	cout << "MS shp=" << shp_D_MS << endl;
	cout << "Contract shp: " << endl;
	//shp_C_MS.t().print();

	return 1;
}
