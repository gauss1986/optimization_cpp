#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stat.h>
#include <txtIO.h>
#include <misc.h>

using namespace std;

/* read text file into 2D vector and feedback N_row/N_co */
mat readtxt(const string& filename, int& N_row, int& N_col, const char sep){
    // open and check file
    ifstream f;
    f.open(filename.c_str());
    if(!f) {
        cout << "Cannot open input file "<< filename.c_str() << ".\n";
    }

	vector<vector<double> > rows = read2Dvec(f, sep, N_row, N_col);

    // output size of data
     cout<< "Size of "<<filename.c_str()<<" is "<<N_row<<" rows, "<<N_col<< " cols.\n";

	// print stats
    // vector< vector<double> > stat = simplestat(rows, N_row, N_col);

    // covert x0,x,y to arma format
	mat m;
    vec2D_to_arma(rows,m);

    cout << "Size of m is " << m.n_rows	<< "x" << m.n_cols << endl;

	return m;
}

mat readtxt_real(const string& filename, int& N_row, int& N_col, const char sep, vector<int>& date, vector<string>& col_names, vector<string>& contracts){
    // open and check file
    ifstream f;
    f.open(filename.c_str());
    if(!f) {
        cout << "Cannot open input file "<< filename.c_str() << ".\n";
    }

	vector<vector<double> > rows = read2Dvec_real(f, sep, N_row, N_col, date, col_names, contracts);

    // output size of data
     cout<< "Size of "<<filename.c_str()<<" is "<<N_row<<" rows, "<<N_col<< " cols.\n";

	// print stats
    // vector< vector<double> > stat = simplestat(rows, N_row, N_col);

    // covert x0,x,y to arma format
	mat m;
    vec2D_to_arma(rows,m);

    cout << "Size of m is " << m.n_rows	<< "x" << m.n_cols << endl;

	return m;
}

/* read 2D vector from file */ 
vector<vector<double> > read2Dvec(ifstream& f, const char sep, int& N_row, int& N_col){
    vector<vector<double> > rows; 
    vector<int> col_per_row;
    string temp;
    N_row = 0;
    while(getline(f,temp)){
        istringstream buffer(temp);
        vector<double> row;
        string temp2;
        int col=0;
        while(getline(buffer,temp2,sep)){
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

	return rows;
}

/* read 2D vector from file real*/ 
vector<vector<double> > read2Dvec_real(ifstream& f, const char sep, int& N_row, int& N_col, vector<int>& date, vector<string>& col_names, vector<string>& contracts){
    vector<vector<double> > rows; 
    vector<int> col_per_row;

    string temp;
    N_row = 0;
    while(getline(f,temp)){
        istringstream buffer(temp);
        vector<double> row;
        string temp2;
        int col=0;
        while(getline(buffer,temp2,sep)){
            if (N_row > 0 && col > 1){
				row.push_back(atof(temp2.c_str()));
			}
			else{
				if(N_row==0){
					col_names.push_back(temp2);
				}
				else{
					if (col==0){
						contracts.push_back(temp2); 
					}
					else{
						date.push_back(atoi(temp2.c_str()));
					}
				}
			}
            col++;
        }
		if (N_row > 0){
        	rows.push_back(row);
		}
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

	N_row = N_row-1;
	N_col = N_col-2;

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

