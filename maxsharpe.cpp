#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <readtxt.h>

int main(int argc, char *argv[])
{
    // load and check files
    std::ifstream f("x0.dat");
    if(!f) {
        std::cout << "Cannot open input file "<< filename << ".\n";
        return 1;
    }

    // read data 
    vector<vector<double>> data(f);

    // output data to verify
    cout << "The data is:\n";
    for (int i=0;i<N_row;i++){
        for (int j=0;j<N_col;j++){
            cout << data[i][j] << " ";
        }
        cout << "\n";
    } 
}
