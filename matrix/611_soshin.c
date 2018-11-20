//! Konstantin Soshin, 611
//! MIPT


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mtrx.h"

#define ERR(errmsg, errcode) \
	{ \
	printf(errmsg); \
	return errcode; \
	}

int main(int argc, char* argv[])
{
	if(argc < 2)ERR("enter input file\n", 1);

	FILE* inpfile = fopen(argv[1], "rb");
	if(!inpfile)ERR("can't open file\n", 2);

	matrix* A = load_mtrx(inpfile);
	if(!A)ERR("can't allocate memory\n", 3);

	printf("A:\n");
	printm(*A);

	int n = A->rows;

	matrix* b = mk_mtrx(n, 1);
	if(!b)ERR("can't make x\n", 4);

	int i = 0;
	for(i = 0; i < n; i++)b->mtrx[i][0] = i + 1;

	printf("b:\n");
	printm(*b);

    
    matrix* f = mul(*A, *b);
	if(!f)ERR("can't multiplicate\n", 5);

    printf("f:\n");
    printm(*f);
    

	matrix* x = solve_SOR(*A, *f, 0.0001, 1.5);

	if(!x)ERR("can't solve\n", 6);

	printf("x:\n");
	printm(*x);

	printf("Successful!\n");

    fclose(inpfile);
    del_mtrx(A);
    del_mtrx(x);
    del_mtrx(b);
    del_mtrx(f);
	return 0;
}
#undef ERR


