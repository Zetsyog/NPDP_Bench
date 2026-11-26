#include <limits.h>
#include <math.h>
#include <mem.h>
#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>

#include "mea.h"

double **Q;
double **Qbp;
double **Pbp;
double **Pu;
double **M;
double *Puu;

float ERT;

char *RNA; // only ACGU

int paired(int i, int j) {
	char nt1 = RNA[i];
	char nt2 = RNA[j];
	if ((nt1 == 'A' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'A') || (nt1 == 'G' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'G') || (nt1 == 'G' && nt2 == 'U') ||
		(nt1 == 'U' && nt2 == 'G')) {
		return 1;
	} else {
		return 0;
	}
}

// n: dataset size, l: minimum loop length.
// Passed as arguments so that they are symbolic parameters in the scop.
static void kernel_mea(int n, int l) {
	int i, j, k, p, q;

	//  compute the partition functions Q and Qbp
#pragma scop
	if (n >= 1 && l >= 0 && l <= 5) {
		for (i = n - 1; i >= 0; i--) {
			for (j = i + 1; j < n; j++) {
				Q[i][j] = Q[i][j - 1];
				for (k = 0; k < j - i - l; k++) {
					Qbp[k + i][j] = Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
					Q[i][j] += Q[i][k + i] * Qbp[k + i][j];
				}
			}
		}
	}

	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			Pbp[i][j] = (Q[0][i] * Q[j][n - 1] * Qbp[i][j]) / Q[0][n - 1]; //  Pbp[i][j] = (Q[1][i]*Q[j+1][N]*Qbp[i][j])/Q[0][N-1];
			for (p = 0; p < i; p++) {
				for (q = j + 1; q < n; q++) {
					Pbp[i][j] += (Pbp[p][q] * ERT * Q[p + 1][i] * Qbp[i][j] * Q[j + 1][q - 1]) / (Qbp[p][q] == 0 ? 1 : Qbp[p][q]);
				}
			}
		}
	}

	for (i = n - 1; i >= 0; i--) {
		for (j = i + 1; j < n; j++) {
			Pu[i][j] = (Q[0][i] * Q[j][n - 1] * 1) / Q[0][n - 1];
			for (p = 0; p < i; p++) {
				for (q = j + 1; q < n; q++) {
					Pu[i][j] += (Pbp[p][q] * ERT * Q[p + 1][i] * 1 * Q[j + 1][q - 1]) / (Qbp[p][q] == 0 ? 1 : Qbp[p][q]);
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
				M[i][j] = MAX(M[i][j], (k > 0 ? M[i][k + i - 1] : 0) + M[k + i + 1][j - 1] + delta * Pbp[k + i][j]) * paired(k + i, j - 1);
			}
			M[i][j] = MAX(M[i][j], M[i][j - 1] + Puu[j - 1]);
		}
	}
#pragma endscop
}

int main(int argc, char *argv[]) {
	srand(RSEED);

	int i;

	ERT = exp((float)-Ebp / (float)RT);

#ifdef RESULT_DUMP
	fprintf(stderr, " -exp(Ebp/RT) = %5.3f\n", ERT);
#endif

	RNA = calloc(DIM, sizeof(char)); // read from FASTA file
	rand_seq(RNA, N);

#ifdef RESULT_DUMP
	fprintf(stderr, "Sequence: ");
	for (i = 0; i < N; i++)
		fprintf(stderr, "%c", RNA[i]);
	fprintf(stderr, "\n\n");
#endif

	Q = memd();
	Qbp = memd();
	Pbp = memd();
	Pu = memd();
	M = memd();
	Puu = malloc(DIM * sizeof(double));

	rna_array_init(Q, 1.0, 1.0);
	rna_array_init(Qbp, 0, 0);
	rna_array_init(Pbp, 0, 0);
	rna_array_init(Pu, 0, 0);
	rna_array_init(M, 0, 0);

	metrics_start_timer;
	kernel_mea(N, MIN_LOOP);
	metrics_stop_timer;
	metrics_print_timer;

#ifdef RESULT_DUMP
	fprintf(stderr, "Q\n");
	rna_array_print(stderr, Q);

	fprintf(stderr, "Qbp\n");
	rna_array_print(stderr, Qbp);

	fprintf(stderr, "Pbp\n");
	rna_array_print(stderr, Pbp);

	fprintf(stderr, "Pu\n");
	rna_array_print(stderr, Pu);

	fprintf(stderr, "Puu\n");
	for (i = 0; i < N - 1; i++)
		fprintf(stderr, "%3.3f ", Puu[i]);
	fprintf(stderr, "\n");

	fprintf(stderr, "M\n");
	rna_array_print(stderr, M);
#endif

	free(Puu);
	free_memd(Q);
	free_memd(Qbp);
	free_memd(Pbp);
	free_memd(Pu);
	free_memd(M);
	free(RNA);

	return 0;
}
