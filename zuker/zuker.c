// Zuker-like RNA folding recurrence (see README.md).

#include <mem.h>
#include <metrics.h>

// n: dataset size (symbolic parameter).
// V[i][j]: energy of i..j when i and j are paired, W[i][j]: energy of i..j.
// EFL, EHF: energy terms of interior loops and hairpin loops.
static void kernel_zuker(int n, int W[NPDP_DIM][NPDP_DIM], int V[NPDP_DIM][NPDP_DIM], int EFL[NPDP_DIM][NPDP_DIM], int EHF[NPDP_DIM][NPDP_DIM]) {
	int i, j, k, m;

#pragma scop
	for (i = n - 1; i >= 0; i--) {
		for (j = i + 1; j < n; j++) {
			for (k = i + 1; k < j; k++) {
				for (m = k + 1; m < j; m++) {
					if (k - i + j - m > 2 && k - i + j - m < 30)
						V[i][j] = NPDP_MIN(V[k][m] + EFL[i][j], V[i][j]);
				}
				W[i][j] += NPDP_MIN(NPDP_MIN(W[i][k], W[k + 1][j]), W[i][j]);
				if (k < j - 1)
					V[i][j] = NPDP_MIN(W[i + 1][k] + W[k + 1][j - 1], V[i][j]);
			}
			V[i][j] = NPDP_MIN(NPDP_MIN(V[i + 1][j - 1], EHF[i][j]), V[i][j]);
			W[i][j] = NPDP_MIN(NPDP_MIN(NPDP_MIN(W[i + 1][j], W[i][j - 1]), V[i][j]), W[i][j]);
		}
	}
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	int(*W)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*V)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*EFL)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*EHF)[NPDP_DIM] = NPDP_ALLOC_2D(int);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			W[i][j] = i * j;
			V[i][j] = i + 1;
			EHF[i][j] = i + 1;
			EFL[i][j] = i + 1;
		}
	}

	NPDP_TIMER_START;
	kernel_zuker(N, W, V, EFL, EHF);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP
	fprintf(stderr, "Matrix V\n");
	fprintf(stderr, "====================\n");
	NPDP_DUMP_2D(stderr, V, N, " %lld ", long long);
	fprintf(stderr, "Matrix W\n");
	fprintf(stderr, "====================\n");
	NPDP_DUMP_2D(stderr, W, N, " %lld ", long long);
#endif

	free(W);
	free(V);
	free(EFL);
	free(EHF);
	return 0;
}
#endif // NPDP_NO_MAIN
