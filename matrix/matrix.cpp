#include <iostream>
#include <cstdio>
#include <exception>
#include <new>
#include <cmath>

#include "vector.h"
#include "matrix.h"

template<typename T> 
matrix<T>::matrix():
    Rows_(VECTOR_CAPACITY),
    Cols_(VECTOR_CAPACITY),
    Mtrx_(vector <vector<T>> (Rows_))
        {
        for(int i = 0; i < Rows_; i++)
            {
            Mtrx_[i] = vector<T>(Cols_);
            }
        }


template<typename T>
matrix<T>::matrix(int rows, int cols):
    Rows_(rows),
    Cols_(cols),
    Mtrx_(vector<vector<T>>(Rows_))
        {
        for(int i = 0; i < Rows_; i++)
            {
            Mtrx_[i] = vector<T>(Cols_);
            }
        }

template<typename T>
matrix<T>::matrix(const matrix & src):
    Rows_(src.Rows_),
    Cols_(src.Cols_),
    Mtrx_(src.Mtrx_)
    {}



template<typename T>
mshape matrix<T>::shape() const
    {
    mshape s;
    s.rows = Rows_;
    s.cols = Cols_;
    return s;
    }

template<typename T>
const matrix<T> & matrix<T>::operator=(const matrix & src)
    {
    Rows_ = src.Rows_;
    Cols_ = src.Cols_;
    Mtrx_ = src.Mtrx_;
    return *this;
    }


template <typename T>
vector<T> & matrix<T>::operator [](const int index)
    {
    if(index >= Rows_)
        {
        vectorErrBadIndex error;
        throw error;
        }
    return Mtrx_[index];
    }


template <typename T>
std::ostream & operator << (std::ostream & out, matrix<T> & m)
    {
    out << "[\n";
    for(int i = 0; i < m.shape().rows; i++)
        {
        for (int j = 0; j < m.shape().cols; j++)
            out << m[i][j] << ' ';
        out << '\n';
        }
    out << "]\n";
    return out;
    }


template <typename T>
std::istream & operator>>(std::istream& in, matrix<T> & m)
    {
    for(int i = 0; i < m.shape().rows; i++)
        for(int j = 0; j < m.shape().cols; j++)
            in >> m[i][j];
    return in;
    }




template <typename T>
matrix<T>  matrix<T>::transpose() 
    {
    matrix<T> res (Cols_, Rows_);
    for(int i = 0; i < Rows_; i++)
        for(int j = 0; j < Cols_; j++)
            res[j][i] = Mtrx_[i][j];
    return res;
    }

std::ostream & operator << (std::ostream & o, const mshape s)
    {
    o << '(' << s.rows << ' ' << s.cols << ')' << '\n';
    return o;
    }

template <typename T>
matrix<T> operator*(matrix<T>& m1, matrix<T>& m2)
    {
    if(m1.shape().cols != m2.shape().rows)
        {
        vectorErrBadShape badsh;
        throw badsh;
        }
    double a;
    matrix<T> res (m1.shape().rows, m2.shape().cols);
    matrix<T> m2_T (m2.transpose());
    //std::cout << res;
    //std::cout << m2_T;
    for(int i = 0; i < m1.shape().rows; i++)
        for(int j = 0; j < m2_T.shape().rows; j++)
            {
            a = res[i][j] = (m1[i]).dot(m2_T[j]);
            //printf("%d\n", a);
            }
    return res;
    }


template <typename T>
matrix<T> solve_gauss(matrix<T> A, matrix<T> b)
    {
    if((A.shape().rows != A.shape().cols) || 
                    (A.shape().rows != b.shape().rows))
        {
        vectorErrBadShape badsh;
        throw badsh;
        }
    int n = A.shape().rows;
    matrix<T> ans(n, 1);
    int i = 0, j = 0, k = 0;
    vector<T> row;
    T val;
    for(i = 0; i < n; i++)
        {
        for(j = i; j < n; j++)
            if(A[j][i] != 0)break;
        if(j == n)return ans;
        row = A[j];
        A[j] = A[i];
        A[i] = row;
        row = b[j];
        b[j] = b[i];
        b[i] = row;
        

        //std::cout <<"!\n" << A << b;

        for(j = i + 1; j < n; j++)
            {
            if(A[j][i])
                {
                val = -A[j][i] / A[i][i];
                for(k = i; k < n; k++) A[j][k] += A[i][k] * val;
                b[j][0] += b[i][0] * val;
                }
            }
        }
    for(i = n - 1; i > -1; i--)
        {
        for(j = i + 1; j < n; j++) b[i][0] -= A[i][j] * ans[j][0];
        ans[i][0] = b[i][0] / A[i][i];
        }

    return ans;
    }


//=============================

//===========
/*
int main()
    { 
    matrix <int> m(2, 2);
    m[0][0] = 1;
    m[1][0] = 6;
    m[0][1] = 0;
    m[1][1] = 1;
    vector<int> s1 = m[0] * 3;
    std::cout<< m[0].eu_norm() ;
    }
*/

