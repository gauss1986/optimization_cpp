#include <vector>
#include "armadillo"

using namespace arma;
using namespace std;

vector<vector<double> > simplestat(const vector<vector<double> >& data, int n_row, int n_col);
//vec resample_vec(vector<vector<double> >& x0_sample, vector<vector<double> >& x_sample, vector<vector<double> >& y_sample, const vector<vector<double> >& x0, const vector<vector<double> >& x, const vector<vector<double> >& y);
void resample_arma(vec& vx0_sample, mat& mx_sample, mat& my_sample, vec& vx0, mat& mx, mat& my);
vec comp_shp(double& shp, const int n, vec& A, mat& mx, mat& my, vec& vx0);
void bsstat(mat mA, vec vx0, mat mx, mat my, vec vshp, mat mshp_contract, const int N, const int n, const int m);
