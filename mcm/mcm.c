// Matrix chain multiplication order (see README.md and Cormen et al.).

#include <mem.h>
#include <metrics.h>

#include <limits.h>

// Matrix A_i has dimension p[i-1] x p[i] for i = 1..n-1 (n = length of p).
// m[i][j]: minimum number of scalar multiplications to compute A_i..A_j.
// s[i][j]: split point k of the optimal parenthesization of A_i..A_j.
// Row and column 0 of m and s are not used.
static void kernel_mcm(int n, int m[NPDP_DIM][NPDP_DIM], int s[NPDP_DIM][NPDP_DIM], int p[NPDP_DIM]) {
	int i, k, L;

#pragma scop
	// L is chain length.
	for (L = 2; L < n; L++) {
		for (i = 1; i < n - L + 1; i++) {
			m[i][i + L - 1] = INT_MAX;
			for (k = i; k <= i + L - 1 - 1; k++) {
				// Record the split before updating the minimum (same semantics as
				// the original minsq helper function, without a hidden side effect).
				s[i][i + L - 1] = (m[i][k] + m[k + 1][i + L - 1] + p[i - 1] * p[k] * p[i + L - 1] < m[i][i + L - 1]) ? k : s[i][i + L - 1];
				m[i][i + L - 1] = NPDP_MIN(m[i][i + L - 1], m[i][k] + m[k + 1][i + L - 1] + p[i - 1] * p[k] * p[i + L - 1]);
			}
		}
	}
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	srand(NPDP_RSEED);

	int *p = NPDP_ALLOC_1D(int);
	int(*m)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*s)[NPDP_DIM] = NPDP_ALLOC_2D(int);

	for (int i = 0; i < N; i++)
		p[i] = (rand() % 20) + 1;

	// cost is zero when multiplying one matrix (m is zero-initialized).
	NPDP_TIMER_START;
	kernel_mcm(N, m, s, p);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP
	fprintf(stderr, "Minimum number of multiplications is %d \n", N > 1 ? m[1][N - 1] : 0);
	fprintf(stderr, "Array s\n");
	fprintf(stderr, "===============\n");
	NPDP_DUMP_UPPER_2D(stderr, s, N, " %d ", int);
#endif

	free(p);
	free(m);
	free(s);
	return 0;
}
#endif // NPDP_NO_MAIN
