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

int main(int argc, char *argv[]) {
	srand(RSEED);

	int i, j, k, ll, p, q;

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
	double *Puu = (double *)malloc(DIM * sizeof(double));

	rna_array_init(Q, 1.0, 1.0);
	rna_array_init(Qbp, 0, 0);
	rna_array_init(Pbp, 0, 0);
	rna_array_init(Pu, 0, 0);
	rna_array_init(M, 0, 0);

	//  compute the partition functions Q and Qbp

	metrics_start_timer
#pragma scop
		if (N >= 1 && l >= 0 && l <= 5) for (i = N - 1; i >= 0; i--) {
		for (j = i + 1; j < N; j++) {
			Q[i][j] = Q[i][j - 1];
			for (k = 0; k < j - i - l; k++) {
				Qbp[k + i][j] = Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
				Q[i][j] += Q[i][k + i] * Qbp[k + i][j];
			}
		}
	}

	for (i = 0; i < N; i++) {
		for (j = i + 1; j < N; j++) {
			Pbp[i][j] = (Q[0][i] * Q[j][N - 1] * Qbp[i][j]) / Q[0][N - 1]; //  Pbp[i][j] = (Q[1][i]*Q[j+1][N]*Qbp[i][j])/Q[0][N-1];
			for (p = 0; p < i; p++) {
				for (q = j + 1; q < N; q++) {
					Pbp[i][j] += (Pbp[p][q] * ERT * Q[p + 1][i] * Qbp[i][j] * Q[j + 1][q - 1]) / (Qbp[p][q] == 0 ? 1 : Qbp[p][q]);
				}
			}
		}
	}

	for (i = N - 1; i >= 0; i--) {
		for (j = i + 1; j < N; j++) {
			Pu[i][j] = (Q[0][i] * Q[j][N - 1] * 1) / Q[0][N - 1];
			for (p = 0; p < i; p++) {
				for (q = j + 1; q < N; q++) {
					Pu[i][j] += (Pbp[p][q] * ERT * Q[p + 1][i] * 1 * Q[j + 1][q - 1]) / (Qbp[p][q] == 0 ? 1 : Qbp[p][q]);
				}
			}
		}
	}

	for (i = 0; i <= N; i++) {
		Puu[i] = 1;
		for (j = i + 1; j < N; j++) {
			Puu[i] += -1 * Pbp[i][j + 1];
		}
		for (k = 0; k < i; k++) {
			Puu[i] += -1 * Pbp[k][i + 1];
		}
	}

	for (i = N - 1; i >= 0; i--) {
		for (j = i + 1; j < N; j++) {
			for (k = 0; k < j - i - l; k++) {
				M[i][j] = MAX(M[i][j], M[i][k + i - 1] + M[k + i + 1][j - 1] + delta * Pbp[k + i][j]) * paired(k + i, j - 1);
			}
			M[i][j] = MAX(M[i][j], M[i][j - 1] + Puu[j - 1]);
		}
	}
#pragma endscop

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
