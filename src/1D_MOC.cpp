#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "1D_MOC.h"
#include "array.h"
#include "../data/LG_QJZ.h"
#include "../data/materials_data.h"

//一维输运方程特征线求解

//计算网格数目
int count_grids(int n_regions, double* region_widths, double mesh_size, int* rgn, int* rgnu, double* rglu)
{
	/*
	传入参数：
	n_regions ：区域数目
	region_widths ：各区域宽度
	mesh_size ：期望网格宽度

	传出数组：
	rgn ：各区域完整网格数目
	rgnu ：各区域不完整网格数目
	rglu ：各区域不完整网格宽度

	返回网格总数
	*/

	int num_grids = 0;
	int i;

	double* xl;
	xl = Make1DArray(n_regions);

	for (i = 0; i < n_regions; i++)
	{
		//计算完整网格数目
		rgn[i] = (int)floor(region_widths[i] / mesh_size);
		xl[i] = region_widths[i] - mesh_size * rgn[i];
		num_grids += rgn[i];

		//判断是否存在非完整网格
		if (xl[i] > 1e-8) {
			rgnu[i] = 1;
			num_grids += 1;
		}
		else
			rgnu[i] = 0;
	}

	vec_copy(xl, rglu, n_regions);

	free(xl);

	return num_grids;
}

//初始化网格结构
int initialize_grids(int n_regions, double* region_widths, double mesh_size, int* rgn, int* rgnu, double* rglu, material* materials, double* grid_centers, double* grid_lengths, material* grid_material)
{
	/*
	传入参数：
	n_regions ：区域数目
	region_widths ：各区域宽度
	mesh_size ：完整网格宽度
	rgn ：各区域完整网格数目
	rgnu ：各区域不完整网格数目
	rglu ：各区域不完整网格宽度
	materials ：各区域材料

	传出数组：
	grid_centers ：各网格中心位置
	grid_lengths ：各网格长度
	grid_material ：各网格材料
	*/

	int i, j;
	int xg = 0;
	double ltt, lre;
	lre = 0;

	for (i = 0; i < n_regions; i++)
	{
		ltt = lre;
		//添加完整的网格
		for (j = 0; j < rgn[i]; j++)
		{
			grid_centers[xg] = ltt + mesh_size / 2;
			grid_lengths[xg] = mesh_size;
			grid_material[xg] = materials[i];
			xg += 1;
			ltt += mesh_size;
		}
		//添加不完整的网格
		if (rgnu[i] == 1)
		{
			grid_centers[xg] = ltt + rglu[i] / 2;
			grid_lengths[xg] = rglu[i];
			grid_material[xg] = materials[i];
			xg += 1;
		}
		lre += region_widths[i];
	}

	return 1;
}


//源迭代求解
int eigen_solver(int ng, int m, int G, double* grid_centers, double* grid_lengths, material* grid_material, int* bdr)
{
	/*
	传入参数：
	ng ：网格数目
	m ：角度离散数目
	G ：能群数目
	grid_centers ：各网格中心位置
	grid_lengths ：各网格长度
	grid_material ：各网格材料
	bdr ：边界条件
	*/

	int i, j, k, g, g2;
	double k_eff = 1;
	double k0 = 1;
	int mh = m / 2;

	double x00, x01, x02;
	double phiin, phiout;
	double err_k;

	double* miu, * omg;
	miu = Make1DArray(m);
	omg = Make1DArray(m);

	double** L;
	L = Make2DArray(m, ng);

	double** phi_lb, ** phi_rb, ** pphi, ** pphi0, ** Q;
	phi_lb = Make2DArray(m, G);
	phi_rb = Make2DArray(m, G);
	pphi = Make2DArray(ng,G);
	pphi0 = Make2DArray(ng, G);
	Q = Make2DArray(ng,G);

	double*** phi0, *** phin;
	phi0 = Make3DArray(m, ng, G);
	phin = Make3DArray(m, ng, G);

	mat_init(phi_lb, m, G, 0.);
	mat_init(phi_rb, m, G, 0.);

	char nom[] = "phi_1D_transport_MOC.dat";

	//获取高斯勒让德求积组
	switch (m)
	{
	case 2:
		miu[0] = LG_S2[0];
		miu[1] = -LG_S2[0];
		omg[0] = LG_S2[1];
		omg[1] = LG_S2[1];
		break;
	case 4:
		for (i = 0; i < mh; i++) {
			miu[i] = LG_S4[i];
			miu[m - i - 1] = -LG_S4[i];
			omg[i] = LG_S4[i + mh];
			omg[m - i - 1] = LG_S4[i + mh];
		}
		break;
	case 8:
		for (i = 0; i < mh; i++) {
			miu[i] = LG_S8[i];
			miu[m - i - 1] = -LG_S8[i];
			omg[i] = LG_S8[i + mh];
			omg[m - i - 1] = LG_S8[i + mh];
		}
		break;
	case 16:
		for (i = 0; i < mh; i++) {
			miu[i] = LG_S16[i];
			miu[m - i - 1] = -LG_S16[i];
			omg[i] = LG_S16[i + mh];
			omg[m - i - 1] = LG_S16[i + mh];
		}
		break;
	default:
		printf("求积组未定义\n");
		goto fail;
	}

	//根据特征线方向及网格计算特征线段长度
	for (j = 0; j < mh; j++)
	{
		for (k = 0; k < ng; k++)
		{
			L[j][k] = grid_lengths[k] / miu[m - j - 1];
			L[m - j - 1][k] = L[j][k];
		}
	}
	
	//phi初始化
	_3Dmat_init(phi0, m, ng,G, 1.);

	//计算初始中子通量密度
	for (k = 0; k < ng; k++)
	{
		for (g = 0; g < G; g++)
		{
			pphi[k][g] = 0;
			for (j = 0; j < m; j++)
				pphi[k][g] += omg[j] * phi0[j][k][g];
		}
	}
	mat_asm(pphi, pphi0, ng, G);

	//计算初始源项
	switch (G)
	{
	case 2:
		for (k = 0; k < ng; k++)
		{
			for (g = 0; g < G; g++)
			{
				//散射源
				x00 = 0;
				for (g2 = 0; g2 < G; g2++)
					x00 += pphi[k][g2] * grid_material->G2.scattering[g2][g];
				//裂变源
				x01 = 0;
				for (g2 = 0; g2 < G; g2++)
					x01 += pphi[k][g2] * grid_material->G2.niu[g2] * grid_material->G2.fission[g2];
				Q[k][g] = x00 / 2 + x01 * grid_material->G2.chi[g] / (2 * k_eff);
			}
		}
		break;
	case 7:
		for (k = 0; k < ng; k++)
		{
			for (g = 0; g < G; g++)
			{
				//散射源
				x00 = 0;
				for (g2 = 0; g2 < G; g2++)
					x00 += pphi[k][g2] * grid_material->G7.scattering[g2][g];
				//裂变源
				x01 = 0;
				for (g2 = 0; g2 < G; g2++)
					x01 += pphi[k][g2] * grid_material->G7.niu[g2] * grid_material->G7.fission[g2];
				Q[k][g] = x00 / 2 + x01 * grid_material->G7.chi[g] / (2 * k_eff);
			}
		}
		break;
	default:
		printf("能群数目与材料不符\n");
		goto fail;
	}

	//计算初始边界条件
	//brd[]为边界条件类型，0为真空边界，1为反射边界
	//左边界
	switch (bdr[0])
	{
	case 0:
		for (j = mh; j < m; j++)
		{
			for (g = 0; g < G; g++)
				phi_lb[j][g] = 0;
		}
		break;
	case 1:
		for (j = mh; j < m; j++)
		{
			for (g = 0; g < G; g++)
				phi_lb[j][g] = phi0[m - j - 1][0][g];
		}
		break;
	default:
		printf("边界条件未定义\n");
		goto fail;
	}
	//右边界
	switch (bdr[1])
	{
	case 0:
		for (j = 0; j < mh; j++)
		{
			for (g = 0; g < G; g++)
				phi_rb[j][g] = 0;
		}
		break;
	case 1:
		for (j = 0; j < mh; j++)
		{
			for (g = 0; g < G; g++)
				phi_rb[j][g] = phi0[m-j-1][ng-1][g];
		}
		break;
	default:
		printf("边界条件未定义\n");
		goto fail;
	}


	//源迭代
	for (i = 1;; i++)
	{
		for (g = 0; g < G; g++)
		{
			//miu<0，从右往左扫
			for (j = 0; j < mh; j++)
			{
				phiin = phi_rb[j][g];
				for (k = ng - 1; k > -1; k -= 1)
				{
					switch (G)
					{
					case 2:
						phiout = phiin * exp(-L[j][k] * grid_material->G2.total[g]) + (1 - exp(-L[j][k] * grid_material->G2.total[g])) * Q[k][g] / grid_material->G2.total[g];
						phin[j][k][g] = Q[k][g] / grid_material->G2.total[g] + (phiin - phiout) / (L[j][k] * grid_material->G2.total[g]);
						phiin = phiout;
						break;
					case 7:
						phiout = phiin * exp(-L[j][k] * grid_material->G7.total[g]) + (1 - exp(-L[j][k] * grid_material->G7.total[g])) * Q[k][g] / grid_material->G7.total[g];
						phin[j][k][g] = Q[k][g] / grid_material->G7.total[g] + (phiin - phiout) / (L[j][k] * grid_material->G7.total[g]);
						phiin = phiout;
						break;
					default:
						goto fail;
					}
				}
				phi_lb[j][g] = phiout;
			}
			//miu>0，从左往右扫
			for (j = mh; j < m; j++)
			{
				phiin = phi_lb[j][g];
				for (k = 0; k < ng; k += 1)
				{
					switch (G)
					{
					case 2:
						phiout = phiin * exp(-L[j][k] * grid_material->G2.total[g]) + (1 - exp(-L[j][k] * grid_material->G2.total[g])) * Q[k][g] / grid_material->G2.total[g];
						phin[j][k][g] = Q[k][g] / grid_material->G2.total[g] + (phiin - phiout) / (L[j][k] * grid_material->G2.total[g]);
						phiin = phiout;
						break;
					case 7:
						phiout = phiin * exp(-L[j][k] * grid_material->G7.total[g]) + (1 - exp(-L[j][k] * grid_material->G7.total[g])) * Q[k][g] / grid_material->G7.total[g];
						phin[j][k][g] = Q[k][g] / grid_material->G7.total[g] + (phiin - phiout) / (L[j][k] * grid_material->G7.total[g]);
						phiin = phiout;
						break;
					default:
						goto fail;
					}
				}
				phi_rb[j][g] = phiout;
			}

			//计算中子通量密度
			for (k = 0; k < ng; k++)
			{
				pphi[k][g] = 0;
				for (j = 0; j < m; j++)
					pphi[k][g] += omg[j] * phin[j][k][g];
			}
		}
		//计算新的有效增殖因子
		switch (G)
		{
		case 2:
			x01 = 0;
			x02 = 0;
			for (k = 0; k < ng; k++)
			{
				for (g = 0; g < G; g++)
				{
					x01 += pphi[k][g] * grid_material->G2.niu[g] * grid_material->G2.fission[g];
					x02 += pphi0[k][g] * grid_material->G2.niu[g] * grid_material->G2.fission[g];
				}
			}
			k_eff = k0 * x01 / x02;
			break;
		case 7:
			x01 = 0;
			x02 = 0;
			for (k = 0; k < ng; k++)
			{
				for (g = 0; g < G; g++)
				{
					x01 += pphi[k][g] * grid_material->G7.niu[g] * grid_material->G7.fission[g];
					x02 += pphi0[k][g] * grid_material->G7.niu[g] * grid_material->G7.fission[g];
				}
			}
			k_eff = k0 * x01 / x02;
			break;
		default:
			goto fail;
		}

		//检验k_eff是否收敛
		err_k = fabs(k_eff - k0) / k0;
		if (err_k < 1e-7)
		{
			printf("源迭代次数：%d\n", i);
			printf("err_k=%e\n", err_k);
			break;
		}

		if (i > 1e4)
		{
			printf("迭代次数超过指定值：%d", i);
			goto fail;
		}

		//计算新的源项
		switch (G)
		{
		case 2:
			for (k = 0; k < ng; k++)
			{
				for (g = 0; g < G; g++)
				{
					//散射源
					x00 = 0;
					for (g2 = 0; g2 < G; g2++)
						x00 += pphi[k][g2] * grid_material->G2.scattering[g2][g];
					//裂变源
					x01 = 0;
					for (g2 = 0; g2 < G; g2++)
						x01 += pphi[k][g2] * grid_material->G2.niu[g2] * grid_material->G2.fission[g2];
					Q[k][g] = x00 / 2 + x01 * grid_material->G2.chi[g] / (2 * k_eff);
				}
			}
			break;
		case 7:
			for (k = 0; k < ng; k++)
			{
				for (g = 0; g < G; g++)
				{
					//散射源
					x00 = 0;
					for (g2 = 0; g2 < G; g2++)
						x00 += pphi[k][g2] * grid_material->G7.scattering[g2][g];
					//裂变源
					x01 = 0;
					for (g2 = 0; g2 < G; g2++)
						x01 += pphi[k][g2] * grid_material->G7.niu[g2] * grid_material->G7.fission[g2];
					Q[k][g] = x00 / 2 + x01 * grid_material->G7.chi[g] / (2 * k_eff);
				}
			}
			break;
		default:
			printf("能群数目与材料不符\n");
			goto fail;
		}

		//处理边界条件
		//左边界
		switch (bdr[0])
		{
		case 0:
			for (j = mh; j < m; j++)
			{
				for (g = 0; g < G; g++)
					phi_lb[j][g] = 0;
			}
			break;
		case 1:
			for (j = mh; j < m; j++)
			{
				for (g = 0; g < G; g++)
					phi_lb[j][g] = phi_lb[m - j - 1][g];
			}
			break;
		default:
			printf("边界条件未定义\n");
			goto fail;
		}
		//右边界
		switch (bdr[1])
		{
		case 0:
			for (j = 0; j < mh; j++)
			{
				for (g = 0; g < G; g++)
					phi_rb[j][g] = 0;
			}
			break;
		case 1:
			for (j = 0; j < mh; j++)
			{
				for (g = 0; g < G; g++)
					phi_rb[j][g] = phi_rb[m - j - 1][g];
			}
			break;
		default:
			printf("边界条件未定义\n");
			goto fail;
		}

		//为下一次迭代做准备
		mat_asm(pphi, pphi0, ng, G);
		_3Dmat_asm(phin, phi0, m, ng, G);
		k0 = k_eff;
	}

	//输出
	printf("k_eff=%f\n", k_eff);
	//输出中子角通量密度至指定文件
	//_1D_transport_SN_GN_output(G, ng, m, phin, grid_centers, miu, nom);
	//输出中子通量密度至指定文件
	_1D_transport_GN_output2(G, ng, pphi, grid_centers, nom);

	free(miu);
	free(omg);

	free2DArray(L, m);
	free2DArray(phi_rb, m);
	free2DArray(phi_lb, m);
	free2DArray(Q, ng);
	free2DArray(pphi, ng);
	free2DArray(pphi0, ng);

	free3DArray(phi0, m, ng);
	free3DArray(phin, m, ng);

	return 1;

fail:
	free(miu);
	free(omg);

	free2DArray(L, m);
	free2DArray(phi_rb, m);
	free2DArray(phi_lb, m);
	free2DArray(Q, ng);
	free2DArray(pphi, ng);
	free2DArray(pphi0, ng);

	free3DArray(phi0, m, ng);
	free3DArray(phin, m, ng);

	return 0;
}

void _1D_transport_GN_output(int G, int n, double** phi, double* x, char* nom)
{
	FILE* fp;
	int k, g;

	fp = fopen(nom, "w+");
	for (g = 0; g < G; g++)
	{
		for (k = 0; k < n; k++)
			fprintf(fp, "%d\t%f\t%f\n", g,x[k], phi[k][g]);
	}
	fclose(fp);
}

void _1D_transport_SN_GN_output(int G, int n, int m, double*** phi, double* x, double* miu, char* nom)
{
	FILE* fp;
	int j, k, g;

	fp = fopen(nom, "w+");
	for (g = 0; g < G; g++)
	{
		for (k = 0; k < n; k++)
		{
			for (j = 0; j < m; j++)
				fprintf(fp, "%d\t%f\t%f\t%f\n", g, x[k], miu[j], phi[j][k][g]);
		}
	}
	fclose(fp);
}

void _1D_transport_GN_output2(int G, int n, double** phi, double* x, char* nom)
{
	FILE* fp;
	int k, g;

	fp = fopen(nom, "w+");

	fprintf(fp,"x\t");
	for (g = 0; g < G; g++)
		fprintf(fp,"G%d\t", g + 1);
	fprintf(fp,"\n");

	for (k = 0; k < n; k++)
	{
		fprintf(fp,"%f\t", x[k]);
		for (g = 0; g < G; g++)
			fprintf(fp,"%f\t", phi[k][g]);
		fprintf(fp,"\n");
	}
	fclose(fp);
}