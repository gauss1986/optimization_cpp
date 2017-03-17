#include "armadillo"

using namespace arma;

std::vector<std::vector<double> > simplestat(const std::vector<std::vector<double> >& data, int N);
vec resample_vec(std::vector<std::vector<double> >& x0_sample, std::vector<std::vector<double> >& x_sample, std::vector<std::vector<double> >& y_sample, const std::vector<std::vector<double> >& x0, const std::vector<std::vector<double> >& x, const std::vector<std::vector<double> >& y);
void resample_arma(vec& sample, vec& vx0_sample, mat& mx_sample, mat& my_sample, vec& vx0, mat& mx, mat& my);
