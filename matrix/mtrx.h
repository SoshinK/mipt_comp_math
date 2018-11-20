#ifndef MTRX_H
#define MTRX_H


typedef struct matrix
	{
	int rows;
	int columns;
	double** mtrx;
	}
matrix;

//! Constructor
matrix* mk_mtrx(const int rows, const int cols);

//! Destructor
void del_mtrx(matrix* m);

//! Copy(to solve ownership problem during assignment) <
matrix* cpy_mtrx(matrix* original);

//! Addition
matrix* add(matrix a, matrix b);

//! Multiplicate with scalar
matrix* mulscalar(matrix m, double scalar);

//! Solve SLE with Gaussian elimination 
matrix* solve_gauss(matrix A, matrix b);


//! Solve SLE with successive over-relaxation method
matrix* solve_SOR(matrix A, matrix b, double residual, double relax_factor);


//! Print matrix in stdin
void printm(const matrix m);

//! Multiplicate two matrices
matrix* mul(const matrix a, const matrix b);

//! Constructor + load from file
matrix* load_mtrx(FILE* in);

#endif
