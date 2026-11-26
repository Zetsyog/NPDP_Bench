// Nussinov RNA folding: maximum number of base pairs (see README.md).

#define NPDP_FP_TYPE long double

#include <mem.h>
#include <metrics.h>

// n: dataset size (symbolic parameter).
// S[i][j]: maximum number of pairs in the subsequence i..j.
// NPDP_PAIRED(x, y) == 1 if nucleotides x and y can pair (function call by
// default, lookup table with -DNPDP_PAIR_TABLE, see include/mem.h).
static void kernel_nussinov(int n, NPDP_FP_TYPE S[NPDP_DIM][NPDP_DIM], NPDP_PAIR_PARAM) {
	int i, j, k;

#pragma scop
	for (i = n - 1; i >= 0; i--) {
		for (j = i + 1; j < n; j++) {
			for (k = 0; k < j - i; k++) {
				S[i][j] = NPDP_MAX(S[i][k + i] + S[k + i + 1][j], S[i][j]);
			}
			for (k = 0; k < 1; k++) {
				S[i][j] = NPDP_MAX(S[i][j], S[i + 1][j - 1] + NPDP_PAIRED(i, j));
			}
		}
	}
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	srand(NPDP_RSEED);

	NPDP_FP_TYPE(*S)[NPDP_DIM] = NPDP_ALLOC_2D(NPDP_FP_TYPE);
	char *rna = NPDP_ALLOC_1D(char);

	npdp_rand_seq(rna, N);
	npdp_pair_t pair = npdp_pair_init(rna);

	NPDP_TIMER_START;
	kernel_nussinov(N, S, pair);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP
	fprintf(stderr, "RNA sequence\n");
	fprintf(stderr, "==================\n");
	npdp_dump_seq(stderr, rna, N);
	fprintf(stderr, "\nArray S\n");
	fprintf(stderr, "==================\n");
	NPDP_DUMP_2D(stderr, S, N, " %.5f ", double);
#endif

	free(S);
	npdp_pair_free(pair);
	free(rna);
	return 0;
}
#endif // NPDP_NO_MAIN
