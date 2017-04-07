#include <stdlib.h>
#include <stdio.h>
#include "armadillo"

using namespace arma;
using namespace std;

mat readtxt(const string& file_name,int& N_row, int& N_col, const char sep);
vector<vector<double> > reorgdata(const vector<vector<double> >& data, const int N_row, const int N_col);
void printdata(const vector<vector<double> >& data, int N_row, int N_col);
template <class T> void printdata_2D(T& data, char* desc, int N_row, int N_col){
    // output content 
    cout << "\n" << desc << endl;
    for (int i=0;i<N_row;i++){    
        for (int j=0;j<N_col;j++){
            cout << " " << data[i][j];
        }
        cout << endl;
    }
}
template <class T> void printdata_1D(T& data, char* desc, int N){
    // output content 
    cout << "\n" << desc << endl;
    for (int i=0;i<N;i++){    
        cout << " " << data[i] << endl;
    }
}
