#ifndef ARRAY_H
#define ARRAY_H

#include <malloc.h>

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

	double* Make1DArray(int m);
	double** Make2DArray(int m, int n);
	void free2DArray(double** a, int m);

	// vector manipulation. should be defined as inline function for higher efficiency
	void   vec_init(double a[], int n, double val);             // a[]=val
	void   vec_add(double c[], double a[], double b[], int n); // c= a + b
	void   vec_minus(double c[], double a[], double b[], int n); // c= a - b
	double vec_dot(double a[], double b[], int n);             // a[].b[] dot product
	double vec_len(double a[], int);                          // length of a vector
	void   vec_cross(double a[], double b[], double c[]); // only for C[3]= A[3] x B[3]; cross product
	double vec_max(double a[], int n);                  // maximum value of a vector
	double vec_min(double a[], int n);                  // minimum value of a vector
	void   vec_copy(double a[], double b[], int n);      // a[] = b[]
	void   vec_normalize(double a[], int n);               // normalize a vector
	void   vec_mul(double a[], double b[], double x, int n);  // a[] = b[] * x
	void   vec_mul_s(double a[], double x, int n);              // a[] *= x


	//��������д��
	/* --------- matrix operation ---------- */
	void mat_init(double** a, int m, int n, double va);//������ֵΪva
	void mat_asm(double** a, double** b, int m, int n);//����ֵ b=a
	void mat_add(double** a, double** b, double** c, int m, int n);//������� a=b+c
	void mat_num_mul(double** a, double** b, int m, int n, double k);//�������� b=ka
	void mat_mul(double** a, double** b, double** c, int m, int n, int l);//����˷� c=ab  aΪm*n�ľ���bΪn*l�ľ���cΪm*l�ľ���

	//��ά����???
	double*** Make3DArray(int m, int n, int p);
	void free3DArray(double*** a, int m, int n);
	void _3Dmat_init(double*** a, int m, int n, int p, double va);
	void _3Dmat_asm(double*** a, double*** b, int m, int n, int p);


#ifdef __cplusplus
}
#endif // __cplusplus

#endif

