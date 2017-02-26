#include "TDMA.h"
namespace matrix
{
	void TDMA(double *A, double *x, double *b, int N) //A*x=b, NÀº matrix size factor
	{
		int i;
		double *P = (double *)calloc(N, sizeof(double));
		double *Q = (double *)calloc(N, sizeof(double));
		P[0] = 0.0;
		Q[0] = b[0];
		for (i = 1; i < N - 1; i++)
		{
			P[i] = -A[2] / ((A[0] * P[i - 1]) + A[1]);
			Q[i] = (b[i] - (A[0] * Q[i - 1])) / ((A[0] * P[i - 1]) + A[1]);
		}
		x[N - 1] = Q[N - 1];
		for (i = N - 2; i > -1; i--)
		{
			x[i] = P[i] * x[i + 1] + Q[i];
		}
		free(P);
		free(Q);
	}
}