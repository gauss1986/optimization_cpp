#include <mkl_lapacke.h>
#include <stdlib.h>
#include <stdio.h>

void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
void free_matrix(double **a);
double **matrix(int n,int m);
double *OLS(const int N, const int n, const int m, std::vector<std::vector<double> >& x0, std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y);
