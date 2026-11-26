// Maximum expected accuracy RNA folding (see README.md).

#include <mem.h>
#include <metrics.h>

#include <math.h>

#include "mea.h"

// n: dataset size, l: minimum loop length (symbolic parameters).
// ert: Boltzmann weight of a base pair, exp(-Ebp/RT).
// NPDP_PAIRED(x, y) == 1 if nucleotides x and y can pair (function call by
// default, lookup table with -DNPDP_PAIR_TABLE, see include/mem.h).
static void kernel_mea(int n, int l, double ert, double Q[NPDP_DIM][NPDP_DIM], double Qbp[NPDP_DIM][NPDP_DIM], double Pbp[NPDP_DIM][NPDP_DIM],
					   double Pu[NPDP_DIM][NPDP_DIM], double M[NPDP_DIM][NPDP_DIM], double Puu[NPDP_DIM], NPDP_PAIR_PARAM) {
	int i, j, k, p, q;

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

	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			Pbp[i][j] = (Q[0][i] * Q[j][n - 1] * Qbp[i][j]) / Q[0][n - 1];
			for (p = 0; p < i; p++) {
				for (q = j + 1; q < n; q++) {
					Pbp[i][j] += (Pbp[p][q] * ert * Q[p + 1][i] * Qbp[i][j] * Q[j + 1][q - 1]) / (Qbp[p][q] == 0 ? 1 : Qbp[p][q]);
				}
			}
		}
	}

	for (i = n - 1; i >= 0; i--) {
		for (j = i + 1; j < n; j++) {
			Pu[i][j] = (Q[0][i] * Q[j][n - 1] * 1) / Q[0][n - 1];
			for (p = 0; p < i; p++) {
				for (q = j + 1; q < n; q++) {
					Pu[i][j] += (Pbp[p][q] * ert * Q[p + 1][i] * 1 * Q[j + 1][q - 1]) / (Qbp[p][q] == 0 ? 1 : Qbp[p][q]);
				}
			}
		}
	}

	for (i = 0; i <= n; i++) {
		Puu[i] = 1;
		for (j = i + 1; j < n; j++) {
			Puu[i] += -1 * Pbp[i][j + 1];
		}
		for (k = 0; k < i; k++) {
			Puu[i] += -1 * Pbp[k][i + 1];
		}
	}

	for (i = n - 1; i >= 0; i--) {
		for (j = i + 1; j < n; j++) {
			for (k = 0; k < j - i - l; k++) {
				// k == 0 would read M[i][i-1], i.e. M[0][-1] when i == 0 (out of bounds in the
				// original code). The lower triangle is 0, so use 0 for that term.
				M[i][j] = NPDP_MAX(M[i][j], (k > 0 ? M[i][k + i - 1] : 0) + M[k + i + 1][j - 1] + NPDP_DELTA * Pbp[k + i][j]) * NPDP_PAIRED(k + i, j - 1);
			}
			M[i][j] = NPDP_MAX(M[i][j], M[i][j - 1] + Puu[j - 1]);
		}
	}
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	srand(NPDP_RSEED);

	// computed in single precision, as in the original code
	float ert = exp((float)-NPDP_EBP / (float)NPDP_RT);

#ifdef NPDP_DUMP
	fprintf(stderr, " -exp(Ebp/RT) = %5.3f\n", ert);
#endif

	char *rna = NPDP_ALLOC_1D(char);
	npdp_rand_seq(rna, N);

#ifdef NPDP_DUMP
	fprintf(stderr, "Sequence: ");
	for (int i = 0; i < N; i++) {
		fprintf(stderr, "%c", rna[i]);
	}
	fprintf(stderr, "\n\n");
#endif

	double (*Q)[NPDP_DIM] = NPDP_ALLOC_2D(double);
	double (*Qbp)[NPDP_DIM] = NPDP_ALLOC_2D(double);
	double (*Pbp)[NPDP_DIM] = NPDP_ALLOC_2D(double);
	double (*Pu)[NPDP_DIM] = NPDP_ALLOC_2D(double);
	double (*M)[NPDP_DIM] = NPDP_ALLOC_2D(double);
	double *Puu = NPDP_ALLOC_1D(double);

	npdp_pair_t pair = npdp_pair_init(rna);

	npdp_rna_array_init(Q, 1.0, 1.0);
	npdp_rna_array_init(Qbp, 0, 0);
	npdp_rna_array_init(Pbp, 0, 0);
	npdp_rna_array_init(Pu, 0, 0);
	npdp_rna_array_init(M, 0, 0);

	NPDP_TIMER_START;
	kernel_mea(N, NPDP_MIN_LOOP, ert, Q, Qbp, Pbp, Pu, M, Puu, pair);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP
	fprintf(stderr, "Q\n");
	NPDP_DUMP_UPPER_2D(stderr, Q, N, " %5.3f ", double);

	fprintf(stderr, "Qbp\n");
	NPDP_DUMP_UPPER_2D(stderr, Qbp, N, " %5.3f ", double);

	fprintf(stderr, "Pbp\n");
	NPDP_DUMP_UPPER_2D(stderr, Pbp, N, " %5.3f ", double);

	fprintf(stderr, "Pu\n");
	NPDP_DUMP_UPPER_2D(stderr, Pu, N, " %5.3f ", double);

	fprintf(stderr, "Puu\n");
	for (int i = 0; i < N - 1; i++) {
		fprintf(stderr, "%3.3f ", Puu[i]);
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "M\n");
	NPDP_DUMP_UPPER_2D(stderr, M, N, " %5.3f ", double);
#endif

	free(Q);
	free(Qbp);
	free(Pbp);
	free(Pu);
	free(M);
	free(Puu);
	npdp_pair_free(pair);
	free(rna);
	return 0;
}
#endif // NPDP_NO_MAIN
