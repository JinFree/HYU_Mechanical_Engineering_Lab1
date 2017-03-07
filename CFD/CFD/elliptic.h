#pragma once
#include "Main.h"
namespace Elliptic
{
	/////////////////////////////////////
	#define Gx		21		//x-grid
	#define Gy		41		//y-grid
	#define GSx		0.05	//x-gridsize
	#define	GSy		0.05	//y-gridsize
	#define T0		100		//Boundary Condition
	#define ermax	0.01	
	
	/////////////////////////////////////
	void main(void);

	void Init_Cond(double *T, double *Tnew);
	void Boundary_Cond(double *T, double *Tnew);

	int schemeselector(void);
	void get_w(void);
	void NAME(char str[50], int scheme_num);
	int compute(double *T, double *Tnew, int scheme);
	int CalJacobi(double *T, double *Tnew, double beta);
	int CalPGS(double *T, double *Tnew, double beta);
	int CalPSOR(double *T, double *Tnew, double beta);
	double CalErr(double *T, double *Tnew);
	void Jacobi(double *T, double *Tnew, double beta);
	void PGS(double *T, double *Tnew, double beta);
	void PSOR(double *T, double *Tnew, double beta);

	void FileWriter(char str[50], double *T, int scheme, int iter);
}