#include <malloc.h>
#include "array.h"
#include <math.h>
#define MAX(x,y)  (x)>(y)?(x):(y)
#define MIN(x,y)  (x)<(y)?(x):(y)

double*  Make1DArray(int m)
{
	double *a;
	a=(double*)malloc(sizeof(double)*m);
	return a;
}

double** Make2DArray(int m,int n)
{
	int i;
	double **a;
	a= (double**)malloc(sizeof(double*)*m);
	for(i=0;i<m;i++)
		a[i]= (double*)malloc(sizeof(double)*n);
	return a;
}

void free2DArray(double **a, int m)
{
	int i;
	for( i=0; i<m; i++ )
		free(a[i]);
	free(a);
}

/* --------- vector operation ---------- */
// vector manipulation. should be defined as inline function
void   vec_init  ( double a[], int n, double val )  // a[0:n-1] = val
{
	int i;
    for( i=0; i<n; i++ )
        a[i] = val;
}

void   vec_add(double *x1, double *x2, double *x3, int n) // x1= x2 + x3
{
	int i;
	for (i = 0; i<n; i++)
		x1[i] = x2[i] + x3[i];
}

void   vec_minus (double *x1, double *x2, double *x3, int n) // x1= x2 - x3
{
	int i;
    for(i=0; i<n; i++ )
        x1[i] = x2[i] - x3[i];
}

double vec_dot   (double *a, double *b, int n)  // Return = a . b
{
	int i;
    double s=0.;
    for( i=0; i<n; i++ )
        s += a[i]*b[i];
    return s;
}

double vec_len   (double *a, int n)
{
	int i;
    double s=0.;
    for( i=0; i<n; i++ )
        s += a[i]*a[i];
    return sqrt(s);
}

void   vec_cross (double a[], double b[], double c[])  // only for a[3]= b[3] x c[3];
{
    a[0]= b[1]*c[2] - b[2]*c[1];
    a[1]=-b[0]*c[2] + b[2]*c[0];
    a[2]= b[0]*c[1] - b[1]*c[0];
}

double vec_max   (double *a, int n)
{
	int i;
    double s=0.;
    for( i=0; i<n; i++ )
        s = MAX( s, a[i] );
    return s;
}

double vec_min (double *a, int n)
{
	int i;
	double vmin = 1.e10;
	for( i=0; i<n; i++ )
		vmin = MIN(vmin, a[i]);
	return vmin;
}

void vec_copy(double *a, double *b, int n)
{
	int i;
	for (i = 0; i < n; i++)
		b[i] = a[i];
}

void vec_normalize(double *a, int n)
{
	int i;
	double s = 0.;
	for (i = 0; i < n; i++)	s += a[i]*a[i];
	s = sqrt(s);
	for (i = 0; i < n; i++) a[i] /= s;
}

void vec_mul(double *a, double *b, double x, int n)
{
	int i;
	for (i = 0; i < n; i++)
		a[i] = b[i] * x;
}

void vec_mul_s(double *a, double x, int n)
{
	int i;
	for (i = 0; i < n; i++)
		a[i] *= x;
}

//以下是新写的
/* --------- matrix operation ---------- */

//矩阵各项赋值为va
void mat_init(double** a, int m, int n, double va)
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			a[i][j] = va;
}

//矩阵赋值 b=a
void mat_asm(double** a, double** b, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			b[i][j] = a[i][j];
}

//矩阵相加 a=b+c
void mat_add(double** a, double** b, double** c, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			a[i][j] = b[i][j] + c[i][j];
}

//矩阵数乘 b=ka
void mat_num_mul(double** a, double** b, int m, int n, double k)
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			b[i][j] = k * a[i][j];
}

//矩阵乘法 c=ab
//a为m*n的矩阵，b为n*l的矩阵，c为m*l的矩阵
void mat_mul(double** a, double** b, double** c, int m, int n, int l)
{
	int i, j, k;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < l; j++)
		{
			c[i][j] = 0;
			for (k = 0; k < n; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
}

//矩阵转置 b=aT
void mat_trp(double** a, double** b, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			b[j][i] = a[i][j];
}

//三维矩阵???
double*** Make3DArray(int m, int n, int p)
{
	int i,j;
	double*** a;
	a = (double***)malloc(sizeof(double**) * m);
	for (i = 0; i < m; i++)
	{
		a[i] = (double**)malloc(sizeof(double*) * n);
		for (j = 0; j < n; j++)
			a[i][j] = (double*)malloc(sizeof(double) * p);
	}
	return a;
}

void free3DArray(double*** a, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
			free(a[i][j]);
		free(a[i]);
	}
	free(a);
}

void _3Dmat_init(double*** a, int m, int n, int p, double va)
{
	int j, k, l;
	for (j = 0; j < m; j++)
	{
		for (k = 0; k < n; k++)
		{
			for (l = 0; l < p; l++)
				a[j][k][l] = va;
		}
	}
}

void _3Dmat_asm(double*** a, double*** b, int m, int n, int p)
{
	int j, k, l;
	for (j = 0; j < m; j++)
	{
		for (k = 0; k < n; k++)
		{
			for (l = 0; l < p; l++)
				b[j][k][l] = a[j][k][l];
		}
	}
}