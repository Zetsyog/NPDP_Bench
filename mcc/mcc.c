// McCaskill partition function of an RNA sequence (see README.md).

#include <mem.h>
#include <metrics.h>

#include <math.h>

#include "mcc.h"

// n: dataset size, l: minimum loop length (symbolic parameters).
// ert: Boltzmann weight of a base pair, exp(-Ebp/RT).
// NPDP_PAIRED(x, y) == 1 if nucleotides x and y can pair (function call by
// default, lookup table with -DNPDP_PAIR_TABLE, see include/mem.h).
static void kernel_mcc(int n, int l, double ert, double Q[NPDP_DIM][NPDP_DIM], double Qbp[NPDP_DIM][NPDP_DIM], NPDP_PAIR_PARAM) {
	int i, j, k;

	//  compute the partition functions Q and Qbp
#pragma scop
	if (n >= 1 && l >= 0 && l <= 5) {
		for (i = n - 1; i >= 0; i--) {
			for (j = i + 1; j < n; j++) {
				Q[i][j] = Q[i][j - 1];
				for (k = 0; k < j - i - l; k++) {
					Qbp[k + i][j] = Q[k + i + 1][j - 1] * ert * NPDP_PAIRED(k + i, j - 1);
					Q[i][j] += Q[i][k + i] * Qbp[k + i][j];
				}
			}
		}
	}
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	srand(NPDP_RSEED);

	// computed in single precision, as in the original code
	float ert = exp((float)-NPDP_EBP / (float)NPDP_RT);

	char *rna = NPDP_ALLOC_1D(char);
	double(*Q)[NPDP_DIM] = NPDP_ALLOC_2D(double);
	double(*Qbp)[NPDP_DIM] = NPDP_ALLOC_2D(double);

	npdp_rand_seq(rna, N);
	npdp_pair_t pair = npdp_pair_init(rna);

	npdp_rna_array_init(Q, 0.4, 0.4);
	npdp_rna_array_init(Qbp, 0.5, 0.5);

	NPDP_TIMER_START;
	kernel_mcc(N, NPDP_MIN_LOOP, ert, Q, Qbp, pair);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP
	fprintf(stderr, "DUMPING RESULTS\n");
	fprintf(stderr, "Array Q\n");
	fprintf(stderr, "===============\n");
	NPDP_DUMP_UPPER_2D(stderr, Q, N, " %5.3f ", double);
	fprintf(stderr, "Array Qbp\n");
	fprintf(stderr, "===============\n");
	NPDP_DUMP_UPPER_2D(stderr, Qbp, N, " %5.3f ", double);
#endif

	free(Q);
	free(Qbp);
	npdp_pair_free(pair);
	free(rna);
	return 0;
}
#endif // NPDP_NO_MAIN
