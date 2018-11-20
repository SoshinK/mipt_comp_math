#include <iostream>
#include <fstream>
#include "vector.h"
#include "vector.cpp"
#include "matrix.h"
#include "matrix.cpp"

int main(int argc, char** argv)
    {
    try
        {
        if(argc < 2)return -1;
        std::ifstream fstream;
        fstream.open(argv[1]);
    
        size_t Arows, Acols;
        fstream >> Arows;
        fstream >> Acols;
    
        matrix<double> A (Arows, Acols);
        matrix<double> b (Arows, 1);
    
        fstream >> A;
        fstream >> b;
        
        std::cout << A << b;

        matrix<double> x = solve_gauss(A, b);

        std::cout << x;

        return 0;
        }
    catch(...)
        {
        std::cout << "Fatal error\n";
        return -2;
        }

    

    }
