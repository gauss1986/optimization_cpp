#include <mkl_lapacke.h>
#include <stdlib.h>
#include <stdio.h>
#include "armadillo"

using namespace arma;

void comp_stat(vec vy, vec vc, mat mx);
void convert_to_arma(const int n, const std::vector<std::vector<double> >& x, const std::vector<double>& y, const double* coeff, int N, vec& vy, vec& vc, mat& mx);
void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
void free_matrix(double **a);
double **matrix(int n,int m);
vec OLS_day(const int N, const int n, const int m, std::vector<std::vector<double> >& x0, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y);
vec OLS_record(const int N, const int n, const int m, mat& mx, mat& my, vec& vx0);
