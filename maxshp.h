#include <stdlib.h>
#include <stdio.h>
#include "armadillo"

using namespace arma;

//arma::vec maxshp(const int N, const int n, const int m, std::vector<std::vector<double> >& x0, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y);
vec maxshp(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0);
