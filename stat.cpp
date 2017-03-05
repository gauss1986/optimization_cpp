#include <numeric>
#include <boost/bind/bind.hpp>
#include <boost/ref.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <stat.h>

std::vector<std::vector<double> > simplestat(const std::vector<std::vector<double> >& data, int N){
    // get stat data
    // mean and std on each dim (for now)
    std::vector<std::vector<double> > stat;
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean,boost::accumulators::tag::variance(boost::accumulators::lazy)> > acc;
    for (int i=0;i<N;i++){
        std::for_each(data[i].begin(), data[i].end(), boost::bind<void>(boost::ref(acc),boost::lambda::_1));
        std::vector<double> stat_1D;
        stat_1D.push_back(boost::accumulators::mean(acc));
        stat_1D.push_back(sqrt(boost::accumulators::variance(acc)));
        stat.push_back(stat_1D);
    }
    
    return stat;
}
