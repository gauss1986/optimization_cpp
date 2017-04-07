#include <iostream>
#include <numeric>
#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/tuple/tuple.hpp>
#include <boost/bind/bind.hpp>
#include <boost/ref.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <stat.h>
#include <txtIO.h>
#include "armadillo"

using namespace arma;
using namespace std;
using namespace boost::accumulators;
using namespace boost::math;

class stat{
    double min, max, mean, median, std, skew, kurt;
    public:
        vector<double>  compute(const vector<double> &data); // compute the statistics
        void report(); // report the statistics
};

vector<double> stat::compute(const vector<double> &data){
    // use boost.minmax to reduce cost associated with min and max
    min = *min_element(data.begin(),data.end());
    max = *max_element(data.begin(),data.end());

    // using boost accumulators library
    accumulator_set<double, features<tag::mean, tag::variance(lazy), tag::median, tag::skewness, tag::kurtosis> > acc;
    for_each(data.begin(), data.end(), boost::bind<void>(boost::ref(acc),boost::lambda::_1)); // put data on each dim into acc
    mean = boost::accumulators::mean(acc); 
    median =  boost::accumulators::median(acc); 
    std = sqrt(boost::accumulators::variance(acc));
    skew = boost::accumulators::skewness(acc);
    kurt = boost::accumulators::kurtosis(acc);
        
    vector<double> stat_1D_val;
    stat_1D_val.push_back(min);
    stat_1D_val.push_back(max);
    stat_1D_val.push_back(mean);
    stat_1D_val.push_back(median);
    stat_1D_val.push_back(std);
    stat_1D_val.push_back(skew);
    stat_1D_val.push_back(kurt);

    return stat_1D_val;
}

void stat::report (){
    // report statistics
    cout << "Statistics:" << endl;

    // modify report details here
    cout << " Min      ";
    cout << min << endl;
    cout << " Max      ";
    cout << max << endl;
    cout << " Mean     ";
    cout << mean << endl;
    cout << " Median   ";
    cout << median << endl;
    cout << " Std      ";
    cout << std << endl;
    cout << " Skewness ";
    cout << skew << endl;
    cout << " Kurtosis ";
    cout << kurt << endl;
    cout << endl;
}

vector<vector<double> > simplestat(const vector<vector<double> >& data, int n_row, int n_col){
    // compute and output mean,std,skew and kurt on each dim
    vector<vector<double> > stat_2D;
    vector<vector<double> > data_col = reorgdata(data,n_row,n_col);
    for (int i=0;i< n_col;i++){
        stat stat_1D;
        vector<double> stat_1D_val = stat_1D.compute(data_col[i]);

        cout << " DOF " << i << endl;
        stat_1D.report();

        stat_2D.push_back(stat_1D_val); // store stat_1D into stat
    }

    return stat_2D;
}

void resample_arma(vec& vx0_sample, mat& mx_sample, mat& my_sample, vec& vx0, mat& mx, mat& my){
    // sample arma form of matrix and vectors
    // sample x0,x and y
    arma_rng::set_seed_random(); 
    int N = vx0.n_elem;
    int n = mx.n_cols;
    vec sample(N); 
    sample.randu();
    sample = sample * N;
    vec sample_processed = arma::trunc(sample);

    for (int i=0;i<vx0.n_rows;i++){
        vx0_sample(i) = vx0(sample_processed(i));
        mx_sample.row(i) = mx.row(sample_processed(i));
        my_sample.row(i) = my.row(sample_processed(i));
    }
}

vec comp_shp(double& shp, const int n, vec& A, mat& mx, mat& my, vec& vx0){
    // compute the sharpe ratio
	int m = my.n_cols; // No. of contracts
	int N = vx0.n_elem; // No. of days

    vec pnl = zeros(N);
    vec shp_contract(m); 
    for (int i=0;i<m;i++){
        mat completex(N,n+1);
        completex.col(0).ones();
        completex.col(1) = vx0;
        completex.col(2) = mx.col(i); 
        vec yfit = conv_to<vec>::from(completex*A);
        vec pnl_contract = my.col(i)%yfit;
        pnl = pnl + pnl_contract;
        shp_contract(i) = arma::mean(pnl_contract)/stddev(pnl_contract);
    }
    shp = arma::mean(pnl)/stddev(pnl);

    return shp_contract;
}

void bsstat(mat mA, vec vx0, mat mx, mat my, vec vshp, mat mshp_contract, const int N, const int n, const int m){
    // compute statistics
    vec tstat(n+1);
    vec q(n+1);
    for (int i=0;i<n+1;i++){
        tstat(i) = arma::mean(mA.col(i))/stddev(mA.col(i));
    }
    students_t dist(N-1); // double check if it is N_bs or N
    for (int i=0;i<n+1;i++){
        q(i) = 2*cdf(complement(dist, fabs(tstat(i))));
    }
    vec sigmasq(m);
    vec R2(m);
    vec R2_adj(m);

    for (int i=0;i<m;i++){
    	// sigmasq
    	vec vy = my.col(i);
    	vec vy_copy(vy);
	mat mx_sym(N,n+1,fill::ones);
	mx_sym.col(1) = vx0;
	mx_sym.col(2) = mx.col(i);
    	vy = vy-mx_sym*arma::mean(mA).t();
    	vy = vy%vy;
    	sigmasq(i) = arma::sum(vy)/(N-n-1); 
    	// R2 stats
    	R2(i) = 1-arma::sum(vy)/var(vy_copy)/(N-1);
    	// adjusted R2 stats
    	R2_adj(i) = R2(i)-(1-R2(i))*n/(N-n-1);
    }

    // report statistics
    cout << " Estimate (mean) " << " t value " << " Pr(>|t|)"<<endl;
    mat report(n+1,3);
    report.col(0) = arma::mean(mA).t();
    report.col(1) = tstat;
    report.col(2) = q;
    report.print();
    cout << "Portfolio Shp " << arma::mean(vshp) << endl;
    cout << "Contract Shp " << endl;
    arma::mean(mshp_contract).print();
    for (int i=0;i<m;i++){
        cout << "Sym " << i << endl;
    	cout << "Residual standard error=" << sigmasq(i) << " on " << N-n-1 << " dof." << endl;
    	cout << "R2=" << R2 (i)<< ", Adj R2=" << R2_adj(i) << endl;
    }
    cout << endl;
}
