// codes corrected by Piotr Blaszynski

#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>

int **c;
char *RNA; // only ACGU

#include "mem.h"

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

	int i, j, k, l = 0;

	// printf(" -exp(Ebp/RT) = %5.3f\n", ERT);

	c = mem();

	for (i = 0; i < DIM; i++) {
		for (j = 0; j < DIM; j++) {
			c[i][j] = i + j;
		}
	}

	RNA = (char *)malloc(DIM * sizeof(char *)); // read from FASTA file
	rand_seq(RNA, DIM);

	metrics_start_timer;

#pragma scop
	for (int i = N - 2; i >= 1; i--)
		for (int j = i + 2; j <= N; j++) {
			for (int k = i; k <= j - 1; k++)
				c[i][j] += paired(k, j) ? c[i][k - 1] + c[k + 1][j - 1] : 0;
			c[i][j] = c[i][j] + c[i][j - 1];
		}
#pragma endscop

	metrics_stop_timer;
	metrics_print_timer;

#ifdef RESULT_DUMP
	fprintf(stderr, "DUMPING RESULTS\n");
	fprintf(stderr, "Array Q\n");
	fprintf(stderr, "===============\n");
	array_print(stderr, c);
#endif

	free_mem(c);
	free(RNA);

	return 0;
}
