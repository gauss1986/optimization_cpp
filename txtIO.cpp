#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <txtIO.h>

using namespace std;

/* read text file into 2D vector and feedback N_row/N_co */
vector<vector<double> > readtxt(const string& filename, int& N_row, int& N_col){
    // open and check file
    ifstream f;
    f.open(filename.c_str());
    if(!f) {
        cout << "Cannot open input file "<< filename.c_str() << ".\n";
    }

    // read 2D vector from file 
    vector<vector<double> > rows; 
    vector<int> col_per_row;
    string temp;
    N_row = 0;
    while(getline(f,temp)){
        istringstream buffer(temp);
        vector<double> row;
        string temp2;
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
    if (i!=N_row-1){
        cout<<"The No. of Col. on each row is NOT the same.\n";
    }

    // output size of data
     cout<< "Size of "<<filename.c_str()<<" is "<<N_row<<" rows, "<<N_col<< " cols.\n";

    //printdata(rows,10,N_col);

    // return data
    return rows;        
}

vector<vector<double> > reorgdata(const vector<vector<double> >& data, const int N_row, const int N_col){
    vector<vector<double> > newdata;
    for (int i=0;i<N_col;i++){
        vector<double> coldata;
        for (int j=0;j<N_row;j++){
            coldata.push_back(data[j][i]);
        }
        newdata.push_back(coldata);
    }
    
    return newdata;
}
    
void printdata(const vector<vector<double> >& data, int N_row, int N_col){
    // output content 
    cout << "Printing the first " << N_row << " rows of the data" << endl;
    for (int i=0;i<N_row;i++){    
        for (int j=0;j<N_col;j++){
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

