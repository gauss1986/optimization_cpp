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
    std::string f_x0("x0.csv");
    std::string f_x("x.csv");
    std::string f_y("y.csv");

    // read data file
    int N_x0_row,N_x0_col;
    int N_x_row,N_x_col;
    int N_y_row,N_y_col;
    std::vector<std::vector<double> >  x0 = readtxt(f_x0,N_x0_row,N_x0_col);
    std::vector<std::vector<double> >  x = readtxt(f_x,N_x_row,N_x_col);
    std::vector<std::vector<double> >  y = readtxt(f_y,N_y_row,N_y_col);
}
