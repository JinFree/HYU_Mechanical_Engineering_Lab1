#include "elliptic.h"
namespace Elliptic
{	
	double w_Elliptic = 1.0;	//relaxation factor
	void main(void)
	{
		int iter=0;
		int scheme = schemeselector();
		char str[50] = "Error";
		NAME(str, scheme);
		double *T = (double *)calloc(sizeof(double), Gx*Gy);
		double *Tnew = (double *)calloc(sizeof(double), Gx*Gy);
		Init_Cond(T, Tnew);
		Boundary_Cond(T, Tnew);
		iter = compute(T, Tnew, scheme);
		FileWriter(str, T, scheme, iter);
		free(T);
		free(Tnew);
	}

	void Init_Cond(double *T, double *Tnew)
	{
		memset(T, 0, Gx*Gy*sizeof(double));
		memset(Tnew, 0, Gx*Gy*sizeof(double));
	}
	void Boundary_Cond(double *T, double *Tnew)
	{
		int i;
		for (i = 0; i < Gy; i++)
		{
			T[i*Gx] = Tnew[i*Gx] = 0;						//[0,y]
			T[(i + 1)*Gx - 1] = Tnew[(i + 1)*Gx - 1] = 0;	//[1,y]
		}
		for (i = 0; i < Gx; i++)
		{
			T[Gx*(Gy - 1) + i] = Tnew[Gx*(Gy - 1) + i] = 0;	//[x,1]
			T[i] = Tnew[i] = T0;							//[x,0]
		}
	}

	int schemeselector(void)
	{
		int scheme;
		printf("Jacobi: 1\n");
		printf("PGS: 2\n");
		printf("PSOR: 3\n");
		while (1)
		{
			printf("What scheme you want to solve?:");
			scanf("%d", &scheme);
			getchar();
			if (scheme >= 1 && scheme <= 3)
			{
				break;
			}
			else
				printf("Input proper value\n");
		}
		return scheme;
	}
	void get_w(void)
	{
		printf("input w, 1.0 ~ 1.9: ");
		scanf("%lf", &w_Elliptic);
		getchar();
	}
	void NAME(char str[50], int scheme_num)
	{
		switch (scheme_num)
		{
		case 1:
			strcpy(str, "Jacobi Method");
			break;
		case 2:
			strcpy(str, "PGS Method");
			break;
		case 3:
			strcpy(str, "PSOR Method");
			get_w();
			break;
		}
	}

	int compute(double *T, double *Tnew, int scheme)
	{
		double beta = GSx / GSy;
		int iter = 0;
		switch (scheme) //Switch를 이용하여 미리 설정한 Scheme에 알맞은 연산을 수행하도록 함
		{
		case 1:
			iter = CalJacobi(T, Tnew, beta);
			break;
		case 2:
			iter = CalPGS(T, Tnew, beta);
			break;
		case 3:
			iter = CalPSOR(T, Tnew, beta);
			break;
		}
		printf("\n");
		return iter;
	}

	double CalErr(double *T, double *Tnew)
	{
		double error = 0.0;
		int i, j;
		for (j = 1; j < Gy - 1; j++)
			for (i = 1; i < Gx - 1; i++)
				error = error + fabs(Tnew[i + j*Gx] - T[i + j*Gx]); //Error계산
		return error;
	}
	int CalJacobi(double *T, double *Tnew, double beta)
	{
		int iter = 0;
		double error = 0.0;
		do{
			Jacobi(T, Tnew, beta);
			error = CalErr(T, Tnew);
			memcpy(T, Tnew, Gx*Gy*sizeof(double)); //T=Tnew연산
			iter++;
			printf("iteration = %d\r", iter);
		} while (error >= ermax);
		return iter;
	}

	int CalPGS(double *T, double *Tnew, double beta)
	{
		int iter = 0;
		double error = 0.0;
		do{
			PGS(T, Tnew, beta);
			error = CalErr(T, Tnew);
			memcpy(T, Tnew, Gx*Gy*sizeof(double)); //T=Tnew연산
			iter++;
			printf("iteration = %d\r", iter);
		} while (error >= ermax);
		return iter;
	}

	int CalPSOR(double *T, double *Tnew, double beta)
	{
		int iter = 0;
		double error = 0.0;
		do{
			PSOR(T, Tnew, beta);
			error = CalErr(T, Tnew);
			memcpy(T, Tnew, Gx*Gy*sizeof(double)); //T=Tnew연산
			iter++;
			printf("iteration = %d\r", iter);
		} while (error >= ermax);
		return iter;
	}

	void Jacobi(double *T, double *Tnew, double beta)
	{
		int i, j;
		for (j = 1; j<Gy - 1; j++)
			for (i = 1; i<Gx - 1; i++)
				Tnew[j*Gx + i] = (T[j*Gx + (i + 1)] + T[j*Gx + (i - 1)] + beta*beta*(T[(j + 1)*Gx + i] + T[(j - 1)*Gx + i])) / (2 * (1 + beta*beta));
	}
	void PGS(double *T, double *Tnew, double beta)
	{
		int i, j;
		for (j = 1; j<Gy - 1; j++)
			for (i = 1; i<Gx - 1; i++)
				Tnew[j*Gx + i] = (T[j*Gx + (i + 1)] + Tnew[j*Gx + (i - 1)] + beta*beta*(T[(j + 1)*Gx + i] + Tnew[(j - 1)*Gx + i])) / (2 * (1 + beta*beta));
	}
	void PSOR(double *T, double *Tnew, double beta)
	{
		int i, j;
		for (j = 1; j<Gy - 1; j++)
			for (i = 1; i<Gx - 1; i++)
				Tnew[j*Gx + i] = (1 - w_Elliptic)*T[j*Gx + i] + (T[j*Gx + (i + 1)] + Tnew[j*Gx + (i - 1)] + beta*beta*(T[(j + 1)*Gx + i] + Tnew[(j - 1)*Gx + i]))*w_Elliptic / (2 * (1 + beta*beta));
	}
	void FileWriter(char str[50], double *T, int scheme, int iter)
	{
		FILE *fp;
		char name[50];
		int i, j;
		if (scheme <= 2)
			sprintf(name, "%s.dat", str);
		else
			sprintf(name, "%s, w=%.3f.dat", str, w_Elliptic);
		fp = fopen(name, "w"); //쓰기권한 획득
		fprintf(fp, "iteration=%d\n", iter);
		fprintf(fp, "zone i=%d,j=%d\n", Gx, Gy);
		for (j = 0; j < Gy; j++)
			for (i = 0; i < Gx; i++)
				fprintf(fp, "%f\t%f\t%f\n", double(i)*GSx, double(j)*GSy, T[i + j*Gx]);
		fclose(fp);
	}
}