#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <readtxt.h>

int main(int argc, char *argv[])
{
    // load and check files
    std::string x0("x0.csv");
    std::string x("x.csv");
    std::string y("y.csv");
    std::ifstream f_x0(x0.c_str());
    if(!f_x0) {
        std::cout << "Cannot open input file "<< x0.c_str() << ".\n";
        return 1;
    }
    std::ifstream f_x(x.c_str());
    if(!f_x) {
        std::cout << "Cannot open input file "<< x.c_str() << ".\n";
        return 1;
    }
    std::ifstream f_y(y.c_str());
    if(!f_y) {
        std::cout << "Cannot open input file "<< y.c_str() << ".\n";
        return 1;
    }

    // read data 
    int N_x0_row;
    int N_x0_col;
    std::vector<std::vector<double> > data = readtxt(x0,N_x0_row,N_x0_col);
    std:: cout<< "Size of x0 is "<<N_x0_row<<" rows, "<<N_x0_col<< " cols.\n";

}
