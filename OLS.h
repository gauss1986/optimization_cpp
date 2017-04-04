#include <stdlib.h>
#include <stdio.h>
#include "armadillo"

using namespace arma;

vec OLS_day(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0);
vec OLS_record(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0);
void comp_stat(vec vy, vec vc, mat mx);
