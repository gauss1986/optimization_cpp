#include <stdlib.h>
#include <stdio.h>
#include "armadillo"

using namespace arma;

void comp_stat(vec vy, vec vc, mat mx);
vec OLS_day(const int N, const int n, const int m, std::vector<std::vector<double> >& x0, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y);
vec OLS_record(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0);
