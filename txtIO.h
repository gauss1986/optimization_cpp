std::vector<std::vector<double> > readtxt(const std::string& file_name,int& N_row, int& N_col);
std::vector<std::vector<double> > reorgdata(const std::vector<std::vector<double> >& data, const int N_row, const int N_col);
void printdata(const std::vector<std::vector<double> >& data, int N_row, int N_col);
template <class T> void printdata_2D(T& data, char* desc, int N_row, int N_col){
    // output content 
    std::cout << "\n" << desc << std::endl;
    for (int i=0;i<N_row;i++){    
        for (int j=0;j<N_col;j++){
            std::cout << " " << data[i][j];
        }
        std::cout << std::endl;
    }
}
template <class T> void printdata_1D(T& data, char* desc, int N){
    // output content 
    std::cout << "\n" << desc << std::endl;
    for (int i=0;i<N;i++){    
        std::cout << " " << data[i] << std::endl;
    }
}
