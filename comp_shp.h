#include <stdlib.h>
#include <stdio.h>
#include "armadillo"

using namespace arma;

vec comp_shp(double& shp, const int N, const int m, const int n, vec A, mat& mx, mat& my, vec& vx0);
