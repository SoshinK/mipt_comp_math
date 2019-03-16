/////////////////////////////////////////////////////////////////////
//! Konstantin Soshin, 2018, MIPT
//! <soshinkonstantinv@gmail.com>
//! <soshin.kv@phystech.edu>
//! <https://github.com/SoshinK>
/////////////////////////////////////////////////////////////////////



#ifndef LINALG_H
#define LINALG_H

#include <iostream>
#include <cstdio>
#include <exception>
#include <new>



enum LINALG_CONSTANTS
    {
    VECTOR_CAPACITY = 16,
    POISON_DATA = -666,
    };


//!==================================================================
//!
//! Class vector<T>. T - typename(float, double, int, etc)
//!
//!==================================================================

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

//! Destructor
        ~vector();
	
//! Operator [] for indexing     
        T & operator [](const int index);

//! Assigning operator
        const vector & operator = (const vector & that);    
 
//! Dot product of two vectors
        const T dot(vector & v);

//! Euclidean norm
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
vector<T> operator * (const vector<T> v, T val);

//! Sum of two vectors
template <typename T>
vector<T>  operator + (vector<T>  v1, vector<T>  v2);

//! Minus
template <typename T>
vector<T>  operator - ( vector<T> & v1, vector<T> & v2);

//! Insertion operator. Applied to an output stream
template <typename T>
std::ostream & operator << (std::ostream & out, vector<T> &  v);    


//! Exception is thrown if given index is greater than size
class vectorErrBadIndex: public std::exception      
    {
    public:
    vectorErrBadIndex(const char* what):
	what_(what)
	{}	
    virtual const char* what() const throw()
        {
        return what_;
        }
    const char* what_;
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



//!==================================================================
//!
//! Struct mshape defines the size of the matrix
//! 
//!==================================================================

struct mshape
    {
    int rows;
    int cols;
    };



//!==================================================================
//!
//! Class matrix<T> represents matrix. It is implemented as
//! vector<vector<T>>, where T is typename(int, double etc.)
//!
//!==================================================================

template <typename T>
class matrix
    {
    public:

//! Default constructor
        matrix();

//! Creates matrix with (rows x cols) shape, filled with zeroes
        matrix(int rows, int cols);

//! Copying constructor
        matrix(const matrix & src);

//! Destructor
        ~matrix() {}

//! Operator [] for indexing  
        vector<T> & operator [](const int index);

//! Assigning operator
        const matrix<T> & operator = (const matrix & that);
                
//! Returns transposed matrix                
        matrix<T> transpose();
        
//! Returns shape of the matrix
        mshape shape() const;
        

    private:
        int Rows_;
        int Cols_;
        vector < vector<T> > Mtrx_;

    };


//! Insertion operator for mshape. Allows to print shape in 
//! convenient form. Applied to an output stream.
std::ostream & operator<< (std::ostream & o, const mshape s);


//! Insertion operator for matrix.
template <typename T>
std::ostream & operator<< (std::ostream & out, matrix<T> & m); 

//! Extraction operator for matrix. Matrix should have an 
//! appropriate shape. Applied to an input stream.
template <typename T>
std::istream & operator>> (std::istream& in, matrix<T> & m);

//! Matrix multiplication
template <typename T>
matrix<T> operator* (matrix<T> & m1, matrix<T> & m2);


//! Solve SLE with Gaussian elimination
template <typename T>
matrix<T> solve_gauss(matrix<T>  A, matrix<T> b);


//! Solve SLE with successive over-relaxation method           // TODO
template <typename T>
matrix<T> solve_SOR(matrix<T> A, matrix<T> b, double residual, double relax_factor);


//! Solve system of non-linear equation with Newton's method.
//!
//! *snle - points to function, which correpsonds the original system.
//! *jacob = points to functions, which corresponds the jacobian
//! x0 - start value
//! eps - desired accuracy
//! size - size of the system
vector<double> newton(matrix<double> (*snle)(vector<double> x),
                    matrix<double> (*jacob)(vector<double> x),
                    vector<double> x0,
                    double eps, int size);




template<typename T>
T vmax(vector<T> v);


#endif



