#include <vector>
#include "armadillo"

using namespace arma;
using namespace std;

vector<vector<double> > simplestat(const vector<vector<double> >& data, int n_row, int n_col);
//vec resample_vec(vector<vector<double> >& x0_sample, vector<vector<double> >& x_sample, vector<vector<double> >& y_sample, const vector<vector<double> >& x0, const vector<vector<double> >& x, const vector<vector<double> >& y);
//void resample_arma(vec& vx0_sample, mat& mx_sample, mat& my_sample, vec& vx0, mat& mx, mat& my);
vec comp_shp(double& shp, const int n, const vec& A, const mat& mx, const mat& my, const vec& vx0);
vec comp_shp_real(double& shp_D, const vec& A, const mat& x, const mat& x0, const vec& y, const ivec& date, vector<string>& contracts_all, vector<string>& contracts, const int N);
void bsstat(mat mA, vec vx0, mat mx, mat my, vec vshp, mat mshp_contract, const int N, const int n, const int m);
void bsstat_real(mat mA, mat x, mat x0, vec y, vec vshp, mat mshp_contract, const int N, const int n, const int m);
uvec resample(int N);
