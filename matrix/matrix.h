#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cstdio>
#include <exception>
#include <new>

#include "vector.h"

struct mshape
    {
    int rows;
    int cols;
    };

std::ostream & operator << (std::ostream & o, const mshape s);

template <typename T>
class matrix
    {
    public:
        matrix();
        matrix(int rows, int cols);
        matrix(const matrix & src);
        //matrix(matrix && src);
        ~matrix() {}

        vector<T> & operator [](const int index);

        const matrix<T> & operator = (const matrix & that);
                
        matrix<T> transpose();
        
        mshape shape() const;
        
    private:
        int Rows_;
        int Cols_;
        vector < vector<T> > Mtrx_;

    };


//! 
template <typename T>
std::ostream & operator<< (std::ostream & out, matrix<T> & m); 

template <typename T>
std::istream & operator>>(std::istream& in, matrix<T> & m);

template <typename T>
matrix<T> operator* (matrix<T> & m1, matrix<T> & m2);


//! Solve SLE with Gaussian elimination
template <typename T>
matrix<T> solve_gauss(matrix<T>  A, matrix<T> b);


//! Solve SLE with successive over-relaxation method           // TODO
template <typename T>
matrix<T> solve_SOR(matrix<T> A, matrix<T> b, double residual, double relax_factor);




#endif
