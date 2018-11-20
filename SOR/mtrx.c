#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mtrx.h"


#define SQ(val) \
	(val)*(val)



matrix* mk_mtrx(const int rows, const int cols)
	{
	if(rows < 0 || cols < 0)return NULL;

	matrix* m = (matrix *)calloc(1, sizeof(matrix));
	if(!m)return NULL;

	m->mtrx = (double**)calloc(rows, sizeof(double*));
	if(!(m->mtrx))
		{
		free(m);
		return NULL;
		}
	int i;
	for(i = 0; i < rows; i++)
		{
		m->mtrx[i] = (double*)calloc(cols, sizeof(double));
		if(!m->mtrx[i])
			{
			for(i = i - 1; i > -1; i--)free(m->mtrx[i]);
			free(m->mtrx);
			free(m);
			return NULL;
			}
		}
	m->rows = rows;
	m->columns = cols;
	return m;
	}

matrix* solve_gauss(matrix A, matrix b)
	{
	if((A.rows != A.columns) || (A.rows != b.rows))return NULL;
	int n = A.rows;
	int i = 0, j = 0, k = 0;
	double* row;
	double val;
	for (i = 0; i < n; i++)
		{
		for(j = i; j < n; j++)
			if(A.mtrx[j][i] != 0)break;
		if(j == n)return NULL;
		row = A.mtrx[j];
		A.mtrx[j] = A.mtrx[i];
		A.mtrx[i] = row;
		row = b.mtrx[j];
		b.mtrx[j] = b.mtrx[i];
		b.mtrx[i] = row;


		for(j = i + 1; j < n; j++)
			{
			if(A.mtrx[j][i])
				{
				val = -A.mtrx[j][i] / A.mtrx[i][i];
				for(k = i; k < n; k++)A.mtrx[j][k] += A.mtrx[i][k] * val;
				b.mtrx[j][0] += b.mtrx[i][0] * val;
				}
			}
		}


	matrix* ans = mk_mtrx(A.rows, 1);
    if(!ans)return NULL;
	for(i = n - 1; i > -1; i--)
		{
		for(j = i + 1; j < n; j++)b.mtrx[i][0] -= A.mtrx[i][j] * ans->mtrx[j][0];
		ans->mtrx[i][0] = b.mtrx[i][0] / A.mtrx[i][i];
		}

	return ans;
	}

void del_mtrx(matrix* m)
	{
	int i;
	for(i = m->rows - 1; i > -1; i--)free(m->mtrx[i]);
	free(m->mtrx);
	free(m);
	}

void printm(const matrix m)
	{
	int i, j;
	for(i = 0; i < m.rows; i++)
		{
		for(j = 0; j < m.columns; j++)
			printf("%f ", m.mtrx[i][j]);
		printf("\n");
		}
	}

matrix* mul(const matrix a, const matrix b)
	{
	if(a.rows < 0 || a.columns < 0 ||
		b.rows < 0 || b.columns < 0	||
		a.columns != b.rows)return NULL;
	matrix* m = mk_mtrx(a.rows, b.columns);
	if(!m)return NULL;
	int i = 0, j = 0, k = 0;
	for(i = 0; i < m->rows; i++)
		for(j = 0; j < m->columns; j++)
			for(k = 0; k < a.columns; k++)
				m->mtrx[i][j] += a.mtrx[i][k] * b.mtrx[k][j];

	return m;
	}

matrix* load_mtrx(FILE* in)
	{
	int n;
	double curval;
	if(!fscanf(in, "%d", &n))return NULL;
	matrix* m = mk_mtrx(n, n);
	int i, j;
	for(i = 0; i < n; i++)
		{
		for(j = 0; j < n; j++)
			{
			if(fscanf(in, "%lf", &curval) < 1)
				{
				del_mtrx(m);
				return NULL;
				}
			m->mtrx[i][j] = curval;
			}
		}
	return m;
	}


matrix* add(matrix a, matrix b)
    {
    if((a.rows != b.rows) || (a.columns != b.columns))return NULL;
    matrix* sum = mk_mtrx(a.rows, a.columns);
    if(!sum)return NULL;
    int i = 0, j = 0;
    for(i = 0; i < a.rows; i++)
        for(j = 0; j < a.columns; j++)
            sum->mtrx[i][j] = a.mtrx[i][j] + b.mtrx[i][j];
    return sum;
    }

matrix* mulscalar(matrix a, double val)
    {
    matrix* res = mk_mtrx(a.rows, a.columns);
    if(!res)return NULL;
    int i = 0, j = 0;
    for(i = 0; i < a.rows; i++)
        for(j = 0; j < a.columns; j++)
            res->mtrx[i][j] = a.mtrx[i][j] * val;
    return res;
    }


double count_residual(matrix A, matrix u, matrix f)
    {
    matrix* mulAu = mul(A, u);
    if(!mul)return NAN;
    
    matrix* minusF = mulscalar(f, -1.0);
    if(!minusF)
        {
        del_mtrx(mulAu);
        return NAN;
        }
    
    matrix* R = add(*mulAu, *minusF);
    if(!R)
        {
        del_mtrx(mulAu);
        del_mtrx(minusF);
        return NAN;
        }
    
    double r = 0;
    int i = 0;
    
    for(i = 0; i < R->rows; i++)r += SQ(R->mtrx[i][0]);
    del_mtrx(mulAu);
    del_mtrx(minusF);
    del_mtrx(R);
    return sqrt(r);
    }

matrix* solve_SOR(matrix A, matrix b, double residual, double relax_factor)
    {
    if((A.rows != A.columns) || (A.rows != b.rows))return NULL;
    matrix* x = mk_mtrx(A.rows, 1);
    if(!x)return NULL;
    double val;
    int i = 0, j = 0, k = 0; 

    while(count_residual(A, *x, b) > residual)
        { 
        for(i = 0; i < A.rows; i++)
            {
            val = b.mtrx[i][0];

            for(j = 0; j < i; j++) val -= A.mtrx[i][j] * x->mtrx[j][0];
            
            for(j = i + 1; j < A.columns; j++) val -= A.mtrx[i][j] * x->mtrx[j][0];
            
            x->mtrx[i][0] = (1 - relax_factor) * x->mtrx[i][0] + 
                            relax_factor / A.mtrx[i][i] * val;
            }
        }

    return x;
    }
    
    
#undef SQ

