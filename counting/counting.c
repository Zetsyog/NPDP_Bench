// Values grow exponentially with N: use unsigned arithmetic so overflow
// wraps around (well defined) instead of being undefined behavior.
#define _bench_int_type unsigned int

#include <mem.h>
#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>

_bench_int_type **c;
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
static void kernel_counting(int n, int l) {
	int i, j, k;

#pragma scop
	for (i = n - 2; i >= 1; i--) {
		for (j = i + 2; j <= n; j++) {
			for (k = i; k <= j - l; k++) {
				c[i][j] += paired(k, j) ? c[i][k - 1] + c[k + 1][j - 1] : 0;
			}
		}
	}
#pragma endscop
}

int main() {
	srand(RSEED);

	int i, j;

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
	kernel_counting(N, 0);
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
