// set check to 0 to measure time, set to 1 to check validity

#include <mem.h>
#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>

int **c;
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

int main() {
	srand(RSEED);

	int i, j, k, l = 0;

	// printf(" -exp(Ebp/RT) = %5.3f\n", ERT);

	c = mem();

	for (i = 0; i < DIM; i++) {
		for (j = 0; j < DIM; j++) {
			c[i][j] = i + j;
		}
	}

	RNA = malloc(DIM * sizeof(char)); // read from FASTA file
	rand_seq(RNA, DIM);

	metrics_start_timer;
	//  compute the partition functions Q and Qbp
#pragma scop
	for (i = N - 2; i >= 1; i--) {
		for (j = i + 2; j <= N; j++) {
			for (k = i; k <= j - l; k++) {
				c[i][j] += paired(k, j) ? c[i][k - 1] + c[k + 1][j - 1] : 0;
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
	array_print(stderr, c);
#endif

	free_mem(c);
	free(RNA);
	return 0;
}
