// Knuth's optimal binary search tree recurrence (see README.md).

#include <mem.h>
#include <metrics.h>

// n: dataset size (symbolic parameter).
// c: cost table, w: weights.
static void kernel_knuth(int n, int c[NPDP_DIM][NPDP_DIM], int w[NPDP_DIM][NPDP_DIM]) {
	int i, j, k;

#pragma scop
	for (i = n - 1; i >= 1; i--)
		for (j = i + 1; j <= n; j++)
			for (k = i + 1; k < j; k++)
				c[i][j] = NPDP_MIN(c[i][j], w[i][j] + c[i][k] + c[k][j]);
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	int(*c)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*w)[NPDP_DIM] = NPDP_ALLOC_2D(int);

	for (int i = 0; i < NPDP_DIM; i++)
		for (int j = 0; j < NPDP_DIM; j++) {
			c[i][j] = i + j;
			w[i][j] = i - j;
		}

	NPDP_TIMER_START;
	kernel_knuth(N, c, w);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP
	fprintf(stderr, "DUMPING RESULTS\n");
	fprintf(stderr, "Array Q\n");
	fprintf(stderr, "===============\n");
	NPDP_DUMP_2D(stderr, c, N, " %lld ", long long);
#endif

	free(c);
	free(w);
	return 0;
}
#endif // NPDP_NO_MAIN
