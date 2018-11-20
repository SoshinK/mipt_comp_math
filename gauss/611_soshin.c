//! Konstantin Soshin, 611
//! MIPT


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

//! Solve SLE with Gaussian elimination
matrix* solve_gauss(matrix A, matrix b);

//! Print matrix in stdin
void printm(const matrix m);

//! Multiplicate two matrices
matrix* mul(const matrix a, const matrix b);

//! Constructor + load from file
matrix* load_mtrx(FILE* in);

//! Norm of residual between x and y(see task)
double norm(matrix a, matrix b);



#define ERR(errmsg, errcode) \
	{ \
	printf(errmsg); \
	return errcode; \
	}

int main(int argc, char* argv[])
{
	if(argc < 2)ERR("enter input file\n", 1);

	FILE* f = fopen(argv[1], "rb");
	if(!f)ERR("can't open file\n", 2);

	matrix* A = load_mtrx(f);
	if(!A)ERR("can't allocate memory\n", 3);

	printf("A:\n");
	printm(*A);

	int n = A->rows;

	matrix* x = mk_mtrx(n, 1);
	if(!x)ERR("can't make x\n", 4);

	int i = 0;
	for(i = 0; i < n; i++)x->mtrx[i][0] = i + 1;

	printf("x:\n");
	printm(*x);

	matrix* b = mul(*A, *x);

	printf("b:\n");
	printm(*b);

	if(!b)ERR("can't multiplicate\n", 5);

	matrix* y = solve_gauss(*A, *b);

	if(!y)ERR("can't solve\n", 6);

	printf("y:\n");
	printm(*y);

	double nrm = norm(*x, *y);
	if(nrm == -1)ERR("Can't calculate norm\n", 7);

	printf("nrm:\n");
	printf("%f\n", nrm);

	printf("Successful!\n");

    fclose(f);
    del_mtrx(A);
    del_mtrx(x);
    del_mtrx(b);
    del_mtrx(y);
	return 0;
}
#undef ERR



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
		// перемещаеем строку с ненулевым элементом наверх
		for(j = i; j < n; j++)
			if(A.mtrx[j][i] != 0)break;
		if(j == n)return NULL;
		row = A.mtrx[j];
		A.mtrx[j] = A.mtrx[i];
		A.mtrx[i] = row;
		row = b.mtrx[j];
		b.mtrx[j] = b.mtrx[i];
		b.mtrx[i] = row;

		// прямой ход
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

	// обратный ход
	matrix* ans = mk_mtrx(A.rows, 1);
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

#define SQ(val) \
	(val)*(val)

double norm(matrix a, matrix b)
	{
	if((a.rows != b.rows) || a.columns != 1 || b.columns != 1)
		return -1;
	double res = 0;
	int i = 0;
	for(i = 0; i < a.rows; i++)res += SQ(a.mtrx[i][0] - b.mtrx[i][0]);
	return sqrt(res);
	}
#undef SQ
