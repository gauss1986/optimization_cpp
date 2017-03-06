#include <numeric>
#include <boost/bind/bind.hpp>
#include <boost/ref.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <stat.h>

std::vector<std::vector<double> > simplestat(const std::vector<std::vector<double> >& data, int N){
    // compute and output mean,std,skew and kurt on each dim
    std::vector<std::vector<double> > stat;

    // using boost accumulators library
    using namespace boost::accumulators;
    accumulator_set<double, features<tag::mean, tag::variance(lazy), tag::skewness, tag::kurtosis> > acc;
    for (int i=0;i<N;i++){
        std::vector<double> stat_1D;  // initialize stat on each dim

        std::for_each(data[i].begin(), data[i].end(), boost::bind<void>(boost::ref(acc),boost::lambda::_1)); // put data on each dim into acc

        stat_1D.push_back(boost::accumulators::mean(acc)); // store stats of acc in stat_1D
        stat_1D.push_back(sqrt(boost::accumulators::variance(acc)));
        stat_1D.push_back(boost::accumulators::skewness(acc));
        stat_1D.push_back(boost::accumulators::kurtosis(acc));

        stat.push_back(stat_1D); // store stat_1D into stat
    }

    // print out the stats
    std::cout << "Statistics:" << std::endl;
    std::cout << " Mean ";
    for (int i=0;i<N;i++){
        std::cout << stat[i][0] << ",";
    }
    std::cout << std::endl;
    std::cout << " Std ";
    for (int i=0;i<N;i++){
        std::cout << stat[i][1] << ",";
    }
    std::cout << std::endl;
    std::cout << " Skewness ";
    for (int i=0;i<N;i++){
        std::cout << stat[i][2] << ",";
    }
    std::cout << std::endl;
    std::cout << " Kurtosis ";
    for (int i=0;i<N;i++){
        std::cout << stat[i][3] << ",";
    }
    
    return stat;
}
