#include <iostream>
#include <cstdio>
#include <exception>
#include <new>
#include <algorithm>
#include <cmath>

#include "vector.h"

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

/*        
template <typename T>
vector<T>::vector(vector && that):
    Capacity_(that.Capacity_),
    Size_(that.Size_),
    Data_(new T [Capacity_])
    {
    std::copy(that.Data_, that.Data_ + Capacity_, Data_);
    }
*/

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

/*==========

template <typename T>
const vector<T> & vector<T>::operator = (const vector && that)
    {
    if(this == &that)
        {
        return *this;
        }
    this->~vector();
    new(this)vector(that);
    return *this;
    }


//==========*/
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





