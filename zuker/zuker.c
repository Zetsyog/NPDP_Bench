#include <mem.h>
#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>

#define CHECK_VALID 0

int **W;
int **V;
int **EFL;
int **EHF;
int **VMI;
int **VBI;
int **WZ;

int main(int argc, char *argv[]) {
	int i, j, k, m;

	W = mem();
	V = mem();
	EFL = mem();
	EHF = mem();
	VMI = mem();
	VBI = mem();
	WZ = mem();

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			W[i][j] = i * j;
			V[i][j] = i + 1;
			EHF[i][j] = i + 1;
			EFL[i][j] = i + 1;
		}
	}

	metrics_start_timer;

#pragma scop
	for (i = N - 1; i >= 0; i--) {
		for (j = i + 1; j < N; j++) {
			for (k = i + 1; k < j; k++) {
				for (m = k + 1; m < j; m++) {
					if (k - i + j - m > 2 && k - i + j - m < 30)
						V[i][j] = MIN(V[k][m] + EFL[i][j], V[i][j]);
				}
				W[i][j] += MIN(MIN(W[i][k], W[k + 1][j]), W[i][j]);
				if (k < j - 1)
					V[i][j] = MIN(W[i + 1][k] + W[k + 1][j - 1], V[i][j]);
			}
			V[i][j] = MIN(MIN(V[i + 1][j - 1], EHF[i][j]), V[i][j]);
			W[i][j] = MIN(MIN(MIN(W[i + 1][j], W[i][j - 1]), V[i][j]), W[i][j]);
		}
	}
#pragma endscop

	metrics_stop_timer;
	metrics_print_timer;

#ifdef RESULT_DUMP
	fprintf(stderr, "Matrix V\n");
	fprintf(stderr, "====================\n");
	array_print(stderr, V);
	fprintf(stderr, "Matrix W\n");
	fprintf(stderr, "====================\n");
	array_print(stderr, W);
#endif

	free_mem(W);
	free_mem(V);
	free_mem(EFL);
	free_mem(EHF);
	free_mem(VMI);
	free_mem(VBI);
	free_mem(WZ);

	return 0;
}
