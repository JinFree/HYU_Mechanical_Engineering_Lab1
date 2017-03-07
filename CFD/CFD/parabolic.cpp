#include "parabolic.h"
namespace Parabolic
{
	void main(void)
	{
		double dt = 0.000001, diffusion;
		int N = 41;
		double dx = 0.04 / 40.0;
		double dinamic_viscosity, t_end=0.0007;
		int U0 = 40;
		printf("dinamic_viscosity: ");
		scanf("%lf", &dinamic_viscosity);
		getchar();
		diffusion = ( dinamic_viscosity * dt ) / (dx * dx);
		printf("Diffusion = %.3f, dinamic_viscosity = %.3f\n", diffusion, dinamic_viscosity);
		if (diffusion <= 0.5)
			Explicit(diffusion, dx, dt, t_end, N, U0);
		else
		{
			printf("You Cannot Get Right Answer In Explicit When Diffusion Number Larger Than 0.5\n");
			Explicit(diffusion, dx, dt, t_end, N, U0);
		}			
		Implicit(diffusion, dx, dt, t_end, N, U0);
	}

	void Init_Cond(double *U, double *Unew, int N)
	{
		for (int i = 1; i < N; i++)
			U[i] = Unew[i] = 0;
	}
	void Boundary_Cond(double *U, double *Unew, int U0)
	{
		U[0] = Unew[0] = U0;
	}

	void Explicit(double diffusion, double dx, double dt, double t_end, int N, int U0)
	{
		double *U = (double *)malloc(sizeof(double)* N);
		double *Unew = (double *)malloc(sizeof(double)* N);
		Init_Cond(U, Unew, N);
		Boundary_Cond(U, Unew, U0);
		FILE *Parabolic_Explicit;
		char name[100];
		sprintf(name, "FTCS, diffusion=%.3f.dat", diffusion);
		Parabolic_Explicit = fopen(name, "w");
		fprintf(Parabolic_Explicit, "variables=u(m/sec) y(m)\n");
		int iter = (int)(t_end / dt);
		double t;
		for (int i = 0; i <= iter;i++)
		{
			t = i*dt;
			Explicit_Solver(U, Unew, diffusion, t, dx, N, Parabolic_Explicit, i);
		}
		fclose(Parabolic_Explicit);
		free(U);
		free(Unew);
	}
	void Explicit_Solver(double *U, double *Unew, double diffusion, double t, double dx, int N, FILE *Parabolic_Explicit, int iter)
	{
		if ( iter % 100 == 0)
			File_Write(U, N, dx, t, Parabolic_Explicit);
		Explicit_FTCS(U, Unew, diffusion, N);
		Time_Marching(U, Unew, N);
	}
	void Explicit_FTCS(double *U, double *Unew, double diffusion, int N)
	{
		int i;
		for (i = 1; i < N - 1; i++)
		{
			Unew[i] = diffusion*(U[i - 1] + U[i + 1]) + (1.0 - 2.0*diffusion)*U[i];
		}
	}
	
	void Implicit(double diffusion, double dx, double dt, double t_end, int N, int U0)
	{
		double *U = (double *)malloc(sizeof(double)* N);
		double *Unew = (double *)malloc(sizeof(double)* N);
		double *A = (double *)malloc(sizeof(double)* 3);
		Create_AMat_Lassonen(A, diffusion);

		Init_Cond(U, Unew, N);
		Boundary_Cond(U, Unew, U0);
		FILE *Parabolic_Implicit;
		char name[100];
		sprintf(name, "Lassonen, diffusion=%.3f.dat", diffusion);
		Parabolic_Implicit = fopen(name, "w");
		int iter = (int)(t_end / dt);
		double t;
		for (int i = 0; i <= iter; i++)
		{
			t = i*dt;
			Implicit_Solver(U, Unew, A, diffusion, t, dx, N, Parabolic_Implicit, i);
		}
		fclose(Parabolic_Implicit);
		free(U);
		free(Unew);
	}
	void Implicit_Solver(double *U, double *Unew, double *A, double diffusion, double t, double dx, int N, FILE *Parabolic_Implicit, int iter)
	{
		if (iter % 100 == 0)
			File_Write(U, N, dx, t, Parabolic_Implicit);
		matrix::TDMA(A, Unew, U, N);
		Time_Marching(U, Unew, N);
	}
	void Create_AMat_Lassonen(double *A, double diffusion)
	{
		A[0] = -diffusion;
		A[1] = 1.0 + 2.0*diffusion;
		A[2] = -diffusion;
	}

	void Time_Marching(double *U, double *Unew, int N)
	{
		memcpy(U, Unew, N*sizeof(double));
	}

	void File_Write(double *U,int N, double dx, double t, FILE *file)
	{
		fprintf(file, "variables = U(m/s) Y(m)\n");
		fprintf(file, "zone T=\"t=%.6f\"\n", t);
		for (int i = 0; i < N; i++)
			fprintf(file, "%f\t%f\n", U[i], i*dx);
	}
}