#include <math.h>
#include <mem.h>
#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>

#include "mcc.h"

double **Q;
double **Qbp;

float ERT;

unsigned char *RNA; // only ACGU

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

	RNA = calloc(DIM, sizeof(unsigned char)); // read from FASTA file
	rand_seq(RNA, N);

	Q = memd();
	Qbp = memd();

	rna_array_init(Q, 0.4, 0.4);
	rna_array_init(Qbp, 0.5, 0.5);

	metrics_start_timer;
	//  compute the partition functions Q and Qbp
#pragma scop
	if (N >= 1 && l >= 0 && l <= 5) {
		for (i = N - 1; i >= 0; i--) {
			for (j = i + 1; j < N; j++) {
				Q[i][j] = Q[i][j - 1];
				for (k = 0; k < j - i - l; k++) {
					Qbp[k + i][j] = Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
					Q[i][j] += Q[i][k + i] * Qbp[k + i][j];
				}
			}
		}
	}
#pragma endscop
	metrics_stop_timer;
	metrics_print_timer;

#ifdef RESULT_DUMP
	fprintf(stderr, "DUMPING RESULTS\n");
	fprintf(stderr, "Array Q\n");
	fprintf(stderr, "===============\n");
	rna_array_print(stderr, Q);
	fprintf(stderr, "Array Qbp\n");
	fprintf(stderr, "===============\n");
	rna_array_print(stderr, Qbp);
#endif

	//    FILE *f = fopen("client.data", "wb");
	//    fwrite(Q, sizeof(double), sizeof(double)*DIM*DIM, f);

	free_memd(Q);
	free_memd(Qbp);
	free(RNA);

	return 0;
}
