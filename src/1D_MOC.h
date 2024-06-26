#ifndef _1D_MOC_H
#define _1D_MOC_H
#include "materials.h"

//初始化网格结构
int count_grids(int n_regions, double* region_widths, double mesh_size, int* rgn, int* rgnu, double* rglu);
int initialize_grids(int n_regions, double* region_widths, double mesh_size, int* rgn, int* rgnu, double* rglu, material* materials, double* grid_centers, double* grid_lengths, material* grid_material);
int eigen_solver(int ng, int m, int G, double* grid_centers, double* grid_lengths, material* grid_material, int* bdr);

void _1D_transport_GN_output(int G, int n, double** phi, double* x, char* nom);
void _1D_transport_SN_GN_output(int G, int n, int m, double*** phi, double* x, double* miu, char* nom);
void _1D_transport_GN_output2(int G, int n, double** phi, double* x, char* nom);

#endif // !_1D_MOC_H
