// Counting RNA secondary structures, corrected recurrence (see README.md).
// codes corrected by Piotr Blaszynski

// Values grow exponentially with N: use unsigned arithmetic so overflow
// wraps around (well defined) instead of being undefined behavior.
#define NPDP_INT_TYPE unsigned int

#include <mem.h>
#include <metrics.h>

// n: dataset size (symbolic parameter).
// NPDP_PAIRED(x, y) == 1 if nucleotides x and y can pair (function call by
// default, lookup table with -DNPDP_PAIR_TABLE, see include/mem.h).
static void kernel_counting_im(int n, NPDP_INT_TYPE c[NPDP_DIM][NPDP_DIM], NPDP_PAIR_PARAM) {
#pragma scop
	for (int i = n - 2; i >= 1; i--)
		for (int j = i + 2; j <= n; j++) {
			for (int k = i; k <= j - 1; k++)
				c[i][j] += NPDP_PAIRED(k, j) ? c[i][k - 1] + c[k + 1][j - 1] : 0;
			c[i][j] = c[i][j] + c[i][j - 1];
		}
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	srand(NPDP_RSEED);

	NPDP_INT_TYPE(*c)[NPDP_DIM] = NPDP_ALLOC_2D(NPDP_INT_TYPE);
	char *rna = NPDP_ALLOC_1D(char);

	for (int i = 0; i < NPDP_DIM; i++)
		for (int j = 0; j < NPDP_DIM; j++)
			c[i][j] = i + j;

	npdp_rand_seq(rna, NPDP_DIM);
	npdp_pair_t pair = npdp_pair_init(rna);

	NPDP_TIMER_START;
	kernel_counting_im(N, c, pair);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP
	fprintf(stderr, "DUMPING RESULTS\n");
	fprintf(stderr, "Array Q\n");
	fprintf(stderr, "===============\n");
	NPDP_DUMP_2D(stderr, c, N, " %lld ", long long);
#endif

	free(c);
	npdp_pair_free(pair);
	free(rna);
	return 0;
}
#endif // NPDP_NO_MAIN
