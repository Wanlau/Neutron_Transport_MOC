#ifndef _MATERIALS_H
#define _MATERIALS_H
/*
struct material
{
	int group_num;
	double* total;
	double** scattering;
	double* fission;
	double* niu;
	double* chi;
};
*/

//test
struct materialG2
{
	int group_num = 2;
	double total[2];
	double scattering[2][2];
	double fission[2];
	double niu[2];
	double chi[2];
};

struct materialG7
{
	int group_num = 7;
	double total[7];
	double scattering[7][7];
	double fission[7];
	double niu[7];
	double chi[7];
};

union material
{
	materialG2 G2;
	materialG7 G7;
};

#endif // !_MATERIALS_H

