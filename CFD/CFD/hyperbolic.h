#pragma once
#include "Main.h"
namespace Hyperbolic
{
	#define PI					3.14159265359
	#define alpha				250.0
	void main(void);


	int schemeselector(void);
	void CalCFL(int scheme);
	void NAME(char str[50], int scheme_num);

	FILE* FileOpener(char str[50]);
	void FileWriter(FILE *fp, double *U, double t);

	void Init_Cond(double *U, double *Unew, int scheme_num);
	void Init_Linear(double *U, double *Unew);
	void Init_Non(double *U, double *Unew);

	void Boundary_Cond(double *U, double *Unew, int scheme_num);
	void Boundary_Cond_Linear(double *U, double *Unew);
	void Boundary_Cond_Non(double *U, double *Unew);

	void Time_Marching(double *U, double *Unew);

	void SCHEMEOPEN(double *U, double *Unew, int scheme);

	void Ex_Upwind(double *U, double *Unew);
	void Lax_method(double *U, double *Unew);
	void Lax_Wendroff(double *U, double *Unew);

	void Non_Ex_Upwind(double *U, double *Unew);
	void Non_Lax_Method(double *U, double *Unew);
	void Non_Lax_Wendroff(double *U, double *Unew);

}