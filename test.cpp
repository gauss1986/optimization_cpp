#include "armadillo"
#include <vector>
#include <ctime>
 
using namespace arma;
using namespace std;
 
 int main()
 {
   srand((unsigned)time(0));
 
   mat X = randu<mat>(8, 9);
 
   X.print("source");
 
   mat U;
   vec s;
   mat V;
 
   svd(U,s,V,X);        // use standard algorithm by default
 
   U.print("U");
   s.print("s");
   V.print("V");
   return 0;
 }
