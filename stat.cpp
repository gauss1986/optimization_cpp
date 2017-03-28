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
#include <stat.h>
#include "armadillo"

using namespace arma;
using namespace std;
using namespace boost::accumulators;

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

vector<vector<double> > simplestat(const vector<vector<double> >& data, int N){
    // compute and output mean,std,skew and kurt on each dim
    vector<vector<double> > stat_2D;
    for (int i=0;i<N;i++){
        stat stat_1D;
        vector<double> stat_1D_val = stat_1D.compute(data[i]);

        //cout << " DOF " << i << endl;
        stat_1D.report();

        stat_2D.push_back(stat_1D_val); // store stat_1D into stat
    }

    return stat_2D;
}

vec resample_vec(vector<vector<double> >& x0_sample, vector<vector<double> >& x_sample, vector<vector<double> >& y_sample, const vector<vector<double> >& x0, const vector<vector<double> >& x, const vector<vector<double> >& y){
    // sample x0,x and y
    arma_rng::set_seed_random(); 
    int N = x0.size();
    int n = x[0].size();
    vec sample(N); 
    sample.randu();
    sample = sample * N;
    vec sample_processed = arma::trunc(sample);
    //sample_processed.print();
    
    for (int i=0;i<N;i++){
        vector<double> x0_temp;
        vector<double> x_temp;
        vector<double> y_temp;
        for (int j=0;j<n;j++){
            x_temp.push_back(x[sample_processed(i)][j]);
            y_temp.push_back(y[sample_processed(i)][j]);
        }
        x0_temp.push_back(x0[sample_processed(i)][0]);
        x_sample.push_back(x_temp);
        y_sample.push_back(y_temp);
        x0_sample.push_back(x0_temp);
    } 

    return sample_processed;
}

void resample_arma(vec& sample, vec& vx0_sample, mat& mx_sample, mat& my_sample, vec& vx0, mat& mx, mat& my){
    // sample arma form of matrix and vectors
    for (int i=0;i<vx0.n_rows;i++){
        vx0_sample(i) = vx0(sample(i));
        mx_sample.row(i) = mx.row(sample(i));
        my_sample.row(i) = my.row(sample(i));
    }
}
