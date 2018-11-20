#include <iostream>
#include <exception>
#include "linalg.h"
#include "linalg.cpp"


matrix<double> SNLE(vector<double> x);

matrix<double> jacobian(vector<double> x);
   

int main()
    {
    try
        {
        vector<double> start(4);
        start[0] = start[1] = start[2] = start[3] = 1;
        double eps = 1e-7;

        vector<double> answer = newton((*SNLE), (*jacobian), start, eps, 4);

        std::cout <<"Answer [accuracy: " << eps << "] :\n" <<  answer << '\n';
         
        return 0;
        
        }
    catch(...)
        {
        std::cerr << "Something unexpected happened\n";
        }
    }



matrix<double> SNLE(vector<double> x)
    {
    matrix<double> res(4, 1);
    res[0][0] = 4 * x[0] - x[1] + x[2] - x[0] * x[3];
    res[1][0] = x[0] - 2 * x[1] + 3 * x[2] - x[2] * x[3];
    res[2][0] = -x[0] + 3 * x[1] - 2 * x[2] - x[1] * x[3];
    res[3][0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1;

    return res;
    }


matrix<double> jacobian(vector<double> x)
    {
    matrix<double> jacobian (4, 4);
    jacobian[0][0] = 4 - x[3];
    jacobian[0][1] = -1;
    jacobian[0][2] = 1;
    jacobian[0][3] = -x[0];
    jacobian[1][0] = 1;
    jacobian[1][1] = -2;
    jacobian[1][2] = 3 - x[3];
    jacobian[1][3] = -x[2];
    jacobian[2][0] = -1;
    jacobian[2][1] = 3 - x[3];
    jacobian[2][2] = -2;
    jacobian[2][3] = -x[1];
    jacobian[3][0] = 2 * x[0];
    jacobian[3][1] = 2 * x[1];
    jacobian[3][2] = 2 * x[2];
    jacobian[3][3] = 0;
    return jacobian;
    }


