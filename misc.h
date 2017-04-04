#include <vector>
#include <stdlib.h>
#include "armadillo"

using namespace arma;
using namespace std;

void vec2D_to_arma(const vector<vector<double> >& x, mat& mx);
void free_matrix(double **a);
double **matrix(int n,int m);
