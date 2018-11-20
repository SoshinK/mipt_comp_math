#ifndef VECTOR_H
#define VECTOR_H


#include <cstdio>
#include <exception>
#include <iostream>




enum CONSTANTS
    {
    VECTOR_CAPACITY = 16,
    POISON_DATA = -666,
    };

/*************************************/

//! Class vector<T>. T - typename(float, double, int, etc)

template <typename T>
class vector
    {
    public:
	//! Default constructor. Fills empty cells with poison data (POISON_DATA)
        vector();

	//! Copy constructor
        vector(const vector & source);

    //! Constructor with given size
        vector(const int size);


        //vector(vector && src);


        //! Destructor
        ~vector();
	
	//! Operator [] for indexing     
        T & operator [](const int index);

	//! Assigning operator
        const vector & operator = (const vector & that);    
        //const vector & operator = (const vector && that); // 
        //! Dot product of two vectors
        const T dot(vector & v);
        
        double eu_norm() const;

        //! Service function which prints all meta-information
        void dump() const;

        //! Return size and capacity of the vector respectively
        int capacity() const;
        int size() const;
        
        //! Add element at the end of the vector
        void push_back(const T val);

    private:
        int Capacity_;
        int Size_;
        T* Data_;
  };

//! Scalar mutiplication
template <typename T>
vector<T> operator * (const vector<T> & v, T val);

//! Sum of two vectors
template <typename T>
vector<T>  operator + (vector<T> & v1, vector<T> & v2);

//! Minus
template <typename T>
vector<T>  operator - ( vector<T> & v1, vector<T> & v2);

//! 
template <typename T>
std::ostream & operator << (std::ostream & out, vector<T> &  v);    

//! Exception is thrown if given index is greater than size
class vectorErrBadIndex: public std::exception      
    {                                                 
    virtual const char* what() const throw()
        {
        return "Index is out of range";
        }
    };

//! Is thrown if memory can't be allocated
class vectorErrNoMem: public std::exception                 
    {                                                         
    virtual const char* what() const throw()
        {
        return "Can't allocate memory";
        }
    };

//! Is thrown if two vectors have inappropriate sizes
class vectorErrBadShape: public std::exception
    {
    virtual const char* what() const throw()
        {
        return "Vector sizes don't fit";
        }
    };



#endif
