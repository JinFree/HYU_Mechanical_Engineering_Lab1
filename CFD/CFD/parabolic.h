#pragma once
#include "Main.h"
namespace Parabolic
{
	void main(void);
	
	void Init_Cond(double *U, double *Unew, int N);
	void Boundary_Cond(double *U, double *Unew, int U0);

	void Explicit(double diffusion, double dx, double dt, double t_end, int N, int U0);
	void Explicit_Solver(double *U, double *Unew, double diffusion, double t, double dx, int N, FILE *Parabolic_Explicit, int iter);
	void Explicit_FTCS(double *U, double *Unew, double diffusion, int N);

	void Implicit(double diffusion, double dx, double dt, double t_end, int N, int U0);
	void Implicit_Solver(double *U, double *Unew, double *A, double diffusion, double t, double dx, int N, FILE *Parabolic_Implicit, int iter);
	void Create_AMat_Lassonen(double *A, double diffusion);
	
	void Time_Marching(double *U, double *Unew, int N);

	void File_Write(double *U, int N, double dx, double t, FILE *file);
}