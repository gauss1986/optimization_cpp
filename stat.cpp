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

class stat{
    double min, max, mean, median, std, skew, kurt;
    public:
        std::vector<double>  compute(const std::vector<double> &data);
        void report();
};

std::vector<double> stat::compute(const std::vector<double> &data){
    // use boost.minmax to reduce cost associated with min and max
    min = *std::min_element(data.begin(),data.end());
    max = *std::max_element(data.begin(),data.end());

    // using boost accumulators library
    using namespace boost::accumulators;
    accumulator_set<double, features<tag::mean, tag::variance(lazy), tag::median, tag::skewness, tag::kurtosis> > acc;
    std::for_each(data.begin(), data.end(), boost::bind<void>(boost::ref(acc),boost::lambda::_1)); // put data on each dim into acc
    mean = boost::accumulators::mean(acc); 
    std = sqrt(boost::accumulators::variance(acc));
    skew = boost::accumulators::skewness(acc);
    kurt = boost::accumulators::kurtosis(acc);
        
    std::vector<double> stat_1D_val;
    stat_1D_val.push_back(min);
    stat_1D_val.push_back(max);
    stat_1D_val.push_back(mean);
    stat_1D_val.push_back(median);
    stat_1D_val.push_back(std);
    stat_1D_val.push_back(skew);
    stat_1D_val.push_back(kurt);

    return stat_1D_val;
}

void stat:report (){
    // report statistics
    std::cout << "Statistics:" << std::endl;

    // modify report details here
    std::cout << " Mean ";
    for (int i=0;i<N;i++){
        std::cout << mean << ",";
    }
    std::cout << std::endl;
    std::cout << " Std ";
    for (int i=0;i<N;i++){
        std::cout << stat_2D[i][1] << ",";
    }
    std::cout << std::endl;
    std::cout << " Skewness ";
    for (int i=0;i<N;i++){
        std::cout << stat_2D[i][2] << ",";
    }
    std::cout << std::endl;
    std::cout << " Kurtosis ";
    for (int i=0;i<N;i++){
        std::cout << stat_2D[i][3] << ",";
    }
    std::cout << std::endl;
}

std::vector<std::vector<double> > simplestat(const std::vector<std::vector<double> >& data, int N){
    // compute and output mean,std,skew and kurt on each dim
    std::vector<std::vector<double> > stat_2D;
    for (int i=0;i<N;i++){
        stat stat_1D;
        std::vector<double> stat_1D_val = stat_1D.compute(data[i]);
        stat_2D.push_back(stat_1D_val); // store stat_1D into stat
    }

    
    return stat_2D;
}
