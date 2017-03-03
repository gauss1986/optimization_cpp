#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <readtxt.h>

std::vector<std::vector<double> > readtxt(const std::string& filename,int& N_row, int& N_col){
    // read 2D vector from file 
    std::vector<std::vector<double> > rows; 
    std::vector<int> col_per_row;
    std::string temp;
    N_row = 0;
    std::ifstream f;
    f.open(filename.c_str());
    std::cout << "The file name is "<< filename << "\n";
    while(getline(f,temp)){
        std::istringstream buffer(temp);
        std::vector<double> row;
        std::string temp2;
        int col=0;
        while(getline(buffer,temp2,',')){
            row.push_back(atof(temp2.c_str()));
            col++;
        }
        rows.push_back(row);
        col_per_row.push_back(col);
        N_row++; 
    }

    // verify if the No. of col. on each row is the same
    N_col = col_per_row[0];
    int i = 0;
    while(i<N_row-1 && col_per_row[i+1]==N_col){
        N_col = col_per_row[i+1];
        i++;
    }
    if (i==N_row-1){
        std::cout<<"The No. of Col. on each row is the same.\n";
    }
    else{
        std::cout<<"The No. of Col. on each row is NOT the same.\n";
    }

    return rows;        
}
