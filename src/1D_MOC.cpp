#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "1D_MOC.h"
#include "array.h"
#include "../data/LG_QJZ.h"
#include "../data/materials_data.h"

//һά���˷������������

//����������Ŀ
int count_grids(int n_regions, double* region_widths, double mesh_size, int* rgn, int* rgnu, double* rglu)
{
	/*
	���������
	n_regions ��������Ŀ
	region_widths ����������
	mesh_size ������������

	�������飺
	rgn ������������������Ŀ
	rgnu ������������������Ŀ
	rglu ������������������

	������������
	*/

	int num_grids = 0;
	int i;

	double* xl;
	xl = Make1DArray(n_regions);

	for (i = 0; i < n_regions; i++)
	{
		//��������������Ŀ
		rgn[i] = (int)floor(region_widths[i] / mesh_size);
		xl[i] = region_widths[i] - mesh_size * rgn[i];
		num_grids += rgn[i];

		//�ж��Ƿ���ڷ���������
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

//��ʼ������ṹ
int initialize_grids(int n_regions, double* region_widths, double mesh_size, int* rgn, int* rgnu, double* rglu, material* materials, double* grid_centers, double* grid_lengths, material* grid_material)
{
	/*
	���������
	n_regions ��������Ŀ
	region_widths ����������
	mesh_size ������������
	rgn ������������������Ŀ
	rgnu ������������������Ŀ
	rglu ������������������
	materials �����������

	�������飺
	grid_centers ������������λ��
	grid_lengths �������񳤶�
	grid_material �����������
	*/

	int i, j;
	int xg = 0;
	double ltt, lre;
	lre = 0;

	for (i = 0; i < n_regions; i++)
	{
		ltt = lre;
		//�������������
		for (j = 0; j < rgn[i]; j++)
		{
			grid_centers[xg] = ltt + mesh_size / 2;
			grid_lengths[xg] = mesh_size;
			grid_material[xg] = materials[i];
			xg += 1;
			ltt += mesh_size;
		}
		//��Ӳ�����������
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


//Դ�������
int eigen_solver(int ng, int m, int G, double* grid_centers, double* grid_lengths, material* grid_material, int* bdr)
{
	/*
	���������
	ng ��������Ŀ
	m ���Ƕ���ɢ��Ŀ
	G ����Ⱥ��Ŀ
	grid_centers ������������λ��
	grid_lengths �������񳤶�
	grid_material �����������
	bdr ���߽�����
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

	//��ȡ��˹���õ������
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
		printf("�����δ����\n");
		goto fail;
	}

	//���������߷���������������߶γ���
	for (j = 0; j < mh; j++)
	{
		for (k = 0; k < ng; k++)
		{
			L[j][k] = grid_lengths[k] / miu[m - j - 1];
			L[m - j - 1][k] = L[j][k];
		}
	}
	
	//phi��ʼ��
	_3Dmat_init(phi0, m, ng,G, 1.);

	//�����ʼ����ͨ���ܶ�
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

	//�����ʼԴ��
	switch (G)
	{
	case 2:
		for (k = 0; k < ng; k++)
		{
			for (g = 0; g < G; g++)
			{
				//ɢ��Դ
				x00 = 0;
				for (g2 = 0; g2 < G; g2++)
					x00 += pphi[k][g2] * grid_material->G2.scattering[g2][g];
				//�ѱ�Դ
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
				//ɢ��Դ
				x00 = 0;
				for (g2 = 0; g2 < G; g2++)
					x00 += pphi[k][g2] * grid_material->G7.scattering[g2][g];
				//�ѱ�Դ
				x01 = 0;
				for (g2 = 0; g2 < G; g2++)
					x01 += pphi[k][g2] * grid_material->G7.niu[g2] * grid_material->G7.fission[g2];
				Q[k][g] = x00 / 2 + x01 * grid_material->G7.chi[g] / (2 * k_eff);
			}
		}
		break;
	default:
		printf("��Ⱥ��Ŀ����ϲ���\n");
		goto fail;
	}

	//�����ʼ�߽�����
	//brd[]Ϊ�߽��������ͣ�0Ϊ��ձ߽磬1Ϊ����߽�
	//��߽�
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
		printf("�߽�����δ����\n");
		goto fail;
	}
	//�ұ߽�
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
		printf("�߽�����δ����\n");
		goto fail;
	}


	//Դ����
	for (i = 1;; i++)
	{
		for (g = 0; g < G; g++)
		{
			//miu<0����������ɨ
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
			//miu>0����������ɨ
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

			//��������ͨ���ܶ�
			for (k = 0; k < ng; k++)
			{
				pphi[k][g] = 0;
				for (j = 0; j < m; j++)
					pphi[k][g] += omg[j] * phin[j][k][g];
			}
		}
		//�����µ���Ч��ֳ����
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

		//����k_eff�Ƿ�����
		err_k = fabs(k_eff - k0) / k0;
		if (err_k < 1e-7)
		{
			printf("Դ����������%d\n", i);
			printf("err_k=%e\n", err_k);
			break;
		}

		if (i > 1e4)
		{
			printf("������������ָ��ֵ��%d", i);
			goto fail;
		}

		//�����µ�Դ��
		switch (G)
		{
		case 2:
			for (k = 0; k < ng; k++)
			{
				for (g = 0; g < G; g++)
				{
					//ɢ��Դ
					x00 = 0;
					for (g2 = 0; g2 < G; g2++)
						x00 += pphi[k][g2] * grid_material->G2.scattering[g2][g];
					//�ѱ�Դ
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
					//ɢ��Դ
					x00 = 0;
					for (g2 = 0; g2 < G; g2++)
						x00 += pphi[k][g2] * grid_material->G7.scattering[g2][g];
					//�ѱ�Դ
					x01 = 0;
					for (g2 = 0; g2 < G; g2++)
						x01 += pphi[k][g2] * grid_material->G7.niu[g2] * grid_material->G7.fission[g2];
					Q[k][g] = x00 / 2 + x01 * grid_material->G7.chi[g] / (2 * k_eff);
				}
			}
			break;
		default:
			printf("��Ⱥ��Ŀ����ϲ���\n");
			goto fail;
		}

		//����߽�����
		//��߽�
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
			printf("�߽�����δ����\n");
			goto fail;
		}
		//�ұ߽�
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
			printf("�߽�����δ����\n");
			goto fail;
		}

		//Ϊ��һ�ε�����׼��
		mat_asm(pphi, pphi0, ng, G);
		_3Dmat_asm(phin, phi0, m, ng, G);
		k0 = k_eff;
	}

	//���
	printf("k_eff=%f\n", k_eff);
	//������ӽ�ͨ���ܶ���ָ���ļ�
	//_1D_transport_SN_GN_output(G, ng, m, phin, grid_centers, miu, nom);
	//�������ͨ���ܶ���ָ���ļ�
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