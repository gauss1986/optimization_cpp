#include <vector>
#include <stdlib.h>
#include "armadillo"

using namespace std;
using namespace arma;

void vec2D_to_arma(const vector<vector<double> >& x, mat& mx){
    // conver array to arma format
    int N1 = x.size();
    int N2 = x[0].size();
    mx.set_size(N1,N2);
    for (int i=0;i<N1;i++){
        for (int j=0;j<N2;j++){
            mx(i,j) = x[i][j];
        }
    }
}

double **matrix(int n,int m) {
    // 2D matrix
    double **a = new double * [n];
    a[0] = new double [n*m];
    
    for (int i=1; i<n; i++)
        a[i] = &a[0][i*m];
    
    return a;
}

void free_matrix(double **a) {
    // free 2D matrix
    delete[] a[0];
    delete[] a;
}

