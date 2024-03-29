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

uvec resample(int N){
    arma_rng::set_seed_random(); 
    vec sample(N); 
    sample.randu();
    sample = sample * N;
    vec sample_processed = arma::trunc(sample);
	return conv_to<uvec>::from(sample_processed);
}

//void resample_mat(vec& samplepoints, mat& m_sample, mat& m){
//    for (int i=0;i<vx0.n_rows;i++){
//        vx0_sample(i) = vx0(sample_processed(i));
//        mx_sample.row(i) = mx.row(sample_processed(i));
//        my_sample.row(i) = my.row(sample_processed(i));
//    }
//}

vec comp_shp(double& shp, const int n, const vec& A, const mat& mx, const mat& my, const vec& vx0){
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

vec comp_shp_real(double& shp_D, const vec& A, const mat& x, const mat& x0, const vec& y, const ivec& date, vector<string>& contracts_all, vector<string>& contracts, const int N){
	vec p = y%(join_horiz(ones(y.n_elem),join_horiz(x,x0))*A);	
	// compute day shp
	vec p_D(N); 
	uvec ind = find_unique(date);
	for (int i=0;i<N;i++){
		// index for records on the same day
		uvec m_i = find(date==date(ind(i)));	
		// aggregrate profit on days
		p_D(i) = arma::sum(p(m_i));
	}
	shp_D = arma::mean(p_D)/stddev(p_D);	
	// compute shp by contract
	vec shp_C(contracts.size(),fill::zeros);
	vector<string>::iterator i;
	vector<string>::iterator j;
	for (i=contracts.begin();i!=contracts.end();++i){
		int indi = distance(contracts.begin(),i);
		vector<double> temp;
		for (j=contracts_all.begin();j!=contracts_all.end();++j){
			if (strcmp((*i).c_str(),(*j).c_str())==0){
				temp.push_back(p(distance(contracts_all.begin(),j)));	
			}
		}
		vec p_C = conv_to<vec>::from(temp);
		shp_C(indi) = arma::mean(p_C)/stddev(p_C);
	}
	return shp_C;	
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

void bsstat_real(mat mA, mat x, mat x0, vec y, vec vshp, mat mshp_contract, const int N, const int n, const int m){
}
