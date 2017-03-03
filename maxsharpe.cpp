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
    std::ifstream f("x1.csv");
    if(!f) {
        std::cout << "Cannot open input file "<< f << ".\n";
        return 1;
    }

    // read data 
    int N_row;
    int N_col;
    std::vector<std::vector<double> > data = readtxt("x1.csv",N_row,N_col);
    std:: cout<< "Size of data is "<<N_row<<" rows, "<<N_col<< " cols.\n";

    // output data to verify
    std::cout << "The data is:\n";
    for (int i=0;i<N_row;i++){
        for (int j=0;j<N_col;j++){
            std::cout << data[i][j] << " ";
        }
        std::cout << "\n";
    } 
}
