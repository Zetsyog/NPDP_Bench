// Needleman-Wunsch global sequence alignment with general gap penalties
// (see README.md).

#include <mem.h>
#include <metrics.h>

#include <limits.h>

// n: dataset size (symbolic parameter), both sequences have length n.
// F: scoring matrix, m1/m2: best score ending with a vertical/horizontal gap.
// W[k]: penalty of a gap of length k.
// sim[i][j]: similarity of a[i] and b[j] (+1 match, -1 mismatch).
static void kernel_nw(int n, int F[NPDP_DIM][NPDP_DIM], int m1[NPDP_DIM][NPDP_DIM], int m2[NPDP_DIM][NPDP_DIM], int W[NPDP_DIM],
					  int sim[NPDP_DIM][NPDP_DIM]) {
	int i, j, k;

#pragma scop
	for (i = 1; i <= n; i++)
		for (j = 1; j <= n; j++) {
			// Block S
			m1[i][j] = INT_MIN;
			for (k = 1; k <= i; k++)
				m1[i][j] = NPDP_MAX(m1[i][j], F[i - k][j] - W[k]);
			m2[i][j] = INT_MIN;
			for (k = 1; k <= j; k++)
				m2[i][j] = NPDP_MAX(m2[i][j], F[i][j - k] - W[k]);
			// Unlike Smith-Waterman (sw), no clamp to 0: global alignment
			F[i][j] = NPDP_MAX(F[i - 1][j - 1] + sim[i][j], NPDP_MAX(m1[i][j], m2[i][j]));
		}
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	srand(NPDP_RSEED);

	int(*F)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*m1)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*m2)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*sim)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int *W = NPDP_ALLOC_1D(int);
	char *a = NPDP_ALLOC_1D(char);
	char *b = NPDP_ALLOC_1D(char);

	// W is the gap alignment
	W[0] = 2;
	for (int i = 1; i <= N; i++)
		W[i] = i * W[0];

	// Global alignment: aligning a prefix against the empty sequence
	// costs a gap of that length.
	F[0][0] = 0;
	for (int i = 1; i <= N; i++) {
		F[i][0] = -W[i];
		F[0][i] = -W[i];
	}

	// Sequences are 1-indexed in the kernel (a[1..N])
	npdp_rand_seq(a + 1, N);
	npdp_rand_seq(b + 1, N);
	npdp_build_sim_table(a, b, 1, -1, sim);

	NPDP_TIMER_START;
	kernel_nw(N, F, m1, m2, W, sim);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP

	fprintf(stderr, "Array F\n");
	fprintf(stderr, "===============\n");
	NPDP_DUMP_2D(stderr, F, N, " %lld ", long long);

#endif

	free(F);
	free(m1);
	free(m2);
	free(sim);
	free(W);
	free(a);
	free(b);
	return 0;
}
#endif // NPDP_NO_MAIN
