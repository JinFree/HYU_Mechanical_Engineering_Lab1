#include "hyperbolic.h"
namespace Hyperbolic
{
	int scheme = 1;
	int N = 0; // number of grid
	int L = 0; // 각 case 별 domain length
	double dx = 0.1;
	double dt = 0.1;
	double end_t = 1.0;
	double c = 0.0;
	void main(void)
	{
		//scheme 정하면서 N 수정
		scheme = schemeselector();
		int CFL_Checker = 0;
		dt = dt * 2.0;
		for (CFL_Checker = 0; CFL_Checker < 3; CFL_Checker++)
		{
			dt = dt * 0.5;
			CalCFL(scheme);
			char str[50] = "Error";
			NAME(str, scheme);
			N = int(double(L) / dx);
			//N을 이용해서 U 및 Unew 생성
			double *U = (double *)calloc(sizeof(double), N);
			double *Unew = (double *)calloc(sizeof(double), N);
			//Init
			Init_Cond(U, Unew, scheme);
			//BC
			Boundary_Cond(U, Unew, scheme);
			//Filteopen
			FILE *fp;
			fp = FileOpener(str);
			//Filewrite
			double t = 0.0;
			//Calculate & Filewrite
			for (t = 0.0; t < end_t + dt; t += dt)
			{
				printf("\rt=%.3f", t);
				FileWriter(fp, U, t);
				SCHEMEOPEN(U, Unew, scheme);
				Time_Marching(U, Unew);
				Boundary_Cond(U, Unew, scheme);
			}
			printf("\n");
			fclose(fp);
			free(U);
			free(Unew);
		}
	}


	int schemeselector(void)
	{
		int scheme;
		printf("Linear - Explicit Upwind: 1\n");
		printf("Linear - Lax Method: 2\n");
		printf("Linear - Lax Wendroff: 3\n");
		printf("Nonlinear - Explicit Upwind: 4\n");
		printf("Nonlinear - Lax Method: 5\n");
		printf("Nonlinear - Lax Wendroff: 6\n");
		while (1)
		{
			printf("What scheme you want to solve?:");
			scanf("%d", &scheme);
			getchar();
			if (scheme >= 1 && scheme <= 3)
			{
				dx = 5.0;
				dt = 0.02;
				end_t = 1.0;
				L = 400;
				break;
			}
			else if (scheme >= 4 && scheme <= 6)
			{
				dx = 0.1;
				dt = 0.1;
				end_t = 1.8;
				L = 4;
				break;
			}
			else
				printf("Input proper value\n");
		}
		return scheme;
	}
	void CalCFL(int scheme)
	{
		if (scheme >= 1 && scheme <= 3)
		{
			c = alpha * dt / dx;
			printf("alpha = %.1f\tc=%.5f\ndx=%.5f\tdt=%.5f\tent_t=%.5f\n", alpha, c, dx, dt, end_t);
		}
		else if (scheme >= 4 && scheme <= 6)
		{
			c = dt / dx;
			printf("c=%.2f\ndx=%.2f\tdt=%.3f\tent_t=%.2f\n",c,  dx, dt, end_t);
		}
	}
	void NAME(char str[50], int scheme_num)
	{
		switch (scheme_num)
		{
		case 1:
			strcpy(str, "Linear Ex Upwind");
			break;
		case 2:
			strcpy(str, "Linear Lax Method");
			break;
		case 3:
			strcpy(str, "Linear Lax Wendroff");
			break;
		case 4:
			strcpy(str, "Nonlinear Ex Upwind");
			break;
		case 5:
			strcpy(str, "Nonlinear Lax Method");
			break;
		case 6:
			strcpy(str, "Nonlinear Lax Wendroff");
			break;
		}
	}

	FILE* FileOpener(char str[50])
	{
		FILE *fp;
		char name[100];
		sprintf(name, "%s, CFL = %.2f, dt = %.3f.dat", str, c, dt);
		fp = fopen(name, "w");
		return fp;
	}
	void FileWriter(FILE *fp, double *U, double t)
	{
		fprintf(fp, "VARIABLES = x,U\n");
		fprintf(fp, "zone T=\"t=%f\"\n", t);
		int i;
		for (i = 0; i<N; i++)
		{
			fprintf(fp, "%f\t%f\n", double(i*dx), U[i]);
		}
	}

	void Init_Cond(double *U, double *Unew, int scheme_num)
	{
		memset(U, 0, N*sizeof(double));
		memset(Unew, 0, N*sizeof(double));
		if (scheme >= 1 && scheme <= 3)
		{
			Init_Linear(U, Unew);
		}
		else if (scheme >= 4 && scheme <= 6)
		{
			Init_Non(U, Unew);
		}
	}
	void Init_Linear(double *U, double *Unew)
	{
		int i;
		for (i = int(50 / dx); i <= int(110 / dx); i++)
			U[i] = abs(100 * (sin(PI*(i*dx - 50) / 60)));
	}
	void Init_Non(double *U, double *Unew)
	{
		int i;
		for (i = 0; i <= int(2.0 / dx); i++)
			U[i] = 1.0;
	}

	void Boundary_Cond(double *U, double *Unew, int scheme_num)
	{
		if (scheme >= 1 && scheme <= 3)
		{
			Boundary_Cond_Linear(U, Unew);
		}
		else if (scheme >= 4 && scheme <= 6)
		{
			Boundary_Cond_Non(U, Unew);
		}
	}
	void Boundary_Cond_Linear(double *U, double *Unew)
	{
		U[0] = Unew[0] = 0;
		U[N - 1] = Unew[N - 1] = 0;
	}
	void Boundary_Cond_Non(double *U, double *Unew)
	{
		U[0] = Unew[0] = 1;
		U[N - 1] = Unew[N - 1] = 0;
	}

	void Time_Marching(double *U, double *Unew)
	{
		memcpy(U, Unew, N*sizeof(double));
	}

	void SCHEMEOPEN(double *U, double *Unew, int scheme)
	{
		switch (scheme)
		{
		case 1:
			Ex_Upwind(U, Unew);
			break;
		case 2:
			Lax_method(U, Unew);
			break;
		case 3:
			Lax_Wendroff(U, Unew);
			break;
		case 4:
			Non_Ex_Upwind(U, Unew);
			break;
		case 5:
			Non_Lax_Method(U, Unew);
			break;
		case 6:
			Non_Lax_Wendroff(U, Unew);
			break;
		}
	}

	void Ex_Upwind(double *U, double *Unew)
	{
		int i;
		for (i = 1; i< N-1; i++)
			Unew[i] = U[i] - c*(U[i] - U[i - 1]);
	}
	void Lax_method(double *U, double *Unew)
	{
		int i;
		for (i = 1; i < N-1; i++)
			Unew[i] = (U[i + 1] + U[i - 1]) / 2.0 - 0.5*c*(U[i + 1] - U[i - 1]);
	}
	void Lax_Wendroff(double *U, double *Unew)
	{
		int i;
		for (i = 1; i< N - 1; i++)
			Unew[i] = U[i] - 0.5*c*(U[i + 1] - U[i - 1]) + 0.5*c*c*(U[i + 1] - 2 * U[i] + U[i - 1]);
	}

	void Non_Ex_Upwind(double *U, double *Unew)
	{
		int i;
		double *E = (double *)calloc(sizeof(double), N);
		for (i = 0; i<N; i++)
			E[i] = pow(U[i], 2) / 2.0;
		for (i = 1; i<N-1; i++)
			Unew[i] = U[i] - c*(E[i] - E[i - 1]);
		free(E);
	}
	void Non_Lax_Method(double *U, double *Unew)
	{
		int i;
		double *E = (double *)calloc(sizeof(double), N);
		for (i = 0; i<N; i++)
			E[i] = pow(U[i], 2) / 2.0;
		for (i = 1; i<N - 1; i++)
			Unew[i] = 0.5*(U[i + 1] + U[i - 1]) - 0.5*c*(E[i + 1] - E[i - 1]);
		free(E);
	}
	void Non_Lax_Wendroff(double *U, double *Unew)
	{
		int i;
		double *E = (double *)calloc(sizeof(double), N);
		for (i = 0; i<N; i++)
			E[i] = pow(U[i], 2) / 2.0;
		for (i = 1; i<N - 1; i++)
			Unew[i] = U[i] - 0.5*c*(E[i + 1] - E[i - 1]) + 0.25*c*c*((U[i + 1] + U[i])*(E[i + 1] - E[i]) - (U[i] + U[i - 1])*(E[i] - E[i - 1]));
		free(E);
	}
}