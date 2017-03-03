#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <readtxt.h>

vector<vector<double>> data(const ifstream& f){
    // read 2D vector from file 
    vector<vector<double>> rows; 
    vector<int> col_per_row;
    string temp;
    int N_row=0;
    while(getline(f,temp)){
        istringstrem buffer(temp);
        vector<double> row;
        string temp2;
        int col=0;
        while(getline(buffer,temp2,' ')){
            row.push_back(atof(temp2.c_str()))
            col++;
        }
        rows.push_back(row);
        col_per_row.push_back(col);
        N_row++; 
    }

    // verify if the No. of col. on each row is the same
    int N_col = col_per_row[0];
    int i = 0;
    while(i<N_row-1 && col_per_row[i+1]==N_col){
        N_col = col_per_row[i+1];
        i++;
    }
    if (i==N_row-1){
        cout<<"The No. of Col. on each row is the same.\n";
    }
    else{
        cout<<"The No. of Col. on each row is NOT the same.\n";
    }
        
}
