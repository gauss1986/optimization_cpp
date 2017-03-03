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
    std::ifstream f(filename);
    if(!f) {
        std::cout << "Cannot open input file "<< filename << ".\n";
        return 1;
    }

}
