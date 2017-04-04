#include <vector>
#include <stdlib.h>
#include "armadillo"

using namespace std;
using namespace arma;

// conver array to arma format
void vec2D_to_arma(const vector<vector<double> >& x, mat& mx){
    int N1 = x.size();
    int N2 = x[0].size();
    mx.set_size(N1,N2);
    for (int i=0;i<N1;i++){
        for (int j=0;j<N2;j++){
            mx(i,j) = x[i][j];
        }
    }
}

// 2D matrix
double **matrix(int n,int m) {
    double **a = new double * [n];
    a[0] = new double [n*m];
    
    for (int i=1; i<n; i++)
        a[i] = &a[0][i*m];
    
    return a;
}

// free 2D matrix
void free_matrix(double **a) {
    delete[] a[0];
    delete[] a;
}

