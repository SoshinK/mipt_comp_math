/////////////////////////////////////////////////////////////////////
//! Konstantin Soshin, 2018, MIPT
//! <soshinkonstantinv@gmail.com>
//! <soshin.kv@phystech.edu>
//! <https://github.com/SoshinK>
/////////////////////////////////////////////////////////////////////



#include <iostream>
#include <cstdio>
#include <exception>
#include <new>
#include <algorithm>
#include <cmath>

#include "linalg.h"

template <typename T>
vector<T>::vector():
    Capacity_(VECTOR_CAPACITY),
    Size_(0)
        {
        Data_ = new T[Capacity_];
        if(!Data_)
            {
            vectorErrNoMem enomem;
            throw enomem;
            }
        }


template <typename T>
vector<T>::vector(const int size):
    Capacity_(size),
    Size_(size)
        {
        Data_ = new T[Capacity_];
        if(!Data_)
            {
            vectorErrNoMem enomem;
            throw enomem;
            }
        for(int i = 0; i < Capacity_; ++i)
            {
            Data_[i] = 0;
            }
        }

template <typename T>
vector<T>::~vector()
    {
    Capacity_ = POISON_DATA;
    Size_ = POISON_DATA;
    delete [] Data_;
    }

template <typename T>
vector<T>::vector(const vector & that):
    Capacity_(that.Capacity_),
    Size_(that.Size_),
    Data_(new T [Capacity_])
    {
    std::copy(that.Data_, that.Data_ + Capacity_, Data_);
    }

template <typename T>
const vector<T> & vector<T>::operator = (const vector & that)
    {
    if(this == &that)
        {
        return *this;
        }
    this->~vector();
    new(this)vector(that);
    return *this;
    }


template <typename T>
T & vector<T>::operator [](const int index)
    {
    if(index >= Size_)
        {
        vectorErrBadIndex error;
        throw error;
        }
    return Data_[index];
    }

template <typename T>
void vector<T>::push_back(const T val)
    {
    if(Size_ == Capacity_)
        {
        Capacity_ *= 2;
        Data_ = (T*)realloc(Data_, sizeof(T) * Capacity_);
        for(int i = Size_; i < Capacity_; ++i)Data_[i] = POISON_DATA;
        }
    Data_[Size_++] = val;
    } 

template <typename T>
void vector<T>::dump() const
    {
    printf("\n======================================\n");
    printf("Vector:dump\nCapacity_ = %zu\nSize_ = %zu\nData_ = %p\n", Capacity_, Size_, Data_);
    for(int i = 0; i < Capacity_; i++)
        {
        std::cout << '[' << i << ']' << Data_ + i << " - " << Data_[i];
        if(Data_[i] == POISON_DATA)std::cout << " // poisoned data. Cell is empty or broken";
        std::cout << '\n';
        }
    printf("\n");
    }
    

template <typename T>
int vector<T>::capacity() const{return Capacity_;}

template <typename T>
int vector<T>::size() const {return Size_;}

template <typename T>
vector<T>  operator * (vector<T> & v, T val)
    {
    vector<T> res (v);
    for(int i = 0; i < v.size(); i++) res[i] *= val;
    return res;
    }

template <typename T>
vector<T>  operator + ( vector<T> & v1, vector<T> & v2) 
    {
    if(v1.size() != v2.size())
        {
        vectorErrBadShape badsh;
        throw badsh;
        }
    vector<T> res(v1);
    for(int i = 0; i < v1.size(); i++) res[i] += v2[i];
    return res;
    }

template <typename T>
vector<T>  operator - (vector<T> & v1, vector<T> & v2)
    {
    if(v1.size() != v2.size())
        {
        vectorErrBadShape badsh;
        throw badsh;
        }
    vector<T> res(v1);
    for(int i = 0; i < v1.size(); i++) res[i] -= v2[i];
    return res;
    }
   
template <typename T>
std::ostream & operator << (std::ostream & out, vector<T> & v)
    {
    out << "[ ";
    for(int i = 0; i < v.size(); i++) out << v[i] << ' ';
    out << "]";
    return out;
    }

template <typename T> 
const T vector<T>::dot(vector & v)
    {
    if(Size_ != v.size())
        {
        vectorErrBadShape badsh;
        throw badsh;
        }
    T res = T(0);
    for(int i = 0; i < Size_; i++)res += T(Data_[i] * v[i]);
    return res;
    }

template <typename T>
double vector<T>::eu_norm() const
    {
    double res = 0;
    for(int i = 0; i < Size_; i++)
        {
        res += Data_[i] * Data_[i];
        }
     return sqrt(res);
    }


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
    for(int i = 0; i < m1.shape().rows; i++)
        for(int j = 0; j < m2_T.shape().rows; j++)
            a = res[i][j] = (m1[i]).dot(m2_T[j]);
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

vector<double> newton(matrix<double> (*snle)(vector<double> x),
                    matrix<double> (*jacob)(vector<double> x),
                    vector<double> x0,
                    double eps, int size)
    {
    vector<double> x(size);
    x = x0;            
    vector<double> delta_x(size); 
    matrix<double> f(size, 1);
    double residual = eps + 1000;
    matrix<double> W(size, size); //jacobian
    while(residual > eps)
        {
        W = (*jacob)(x);
        f = (*snle)(x);
        delta_x = ((solve_gauss(W, f)).transpose())[0];
        x = x - delta_x;
        residual = delta_x.eu_norm();
        }
    return x;
    }



