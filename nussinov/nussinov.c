// nussinov

#define _bench_fp_type long double

#include <math.h>
#include <mem.h>
#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

_bench_fp_type **S;
char *RNA;

int can_pair(char *input, int a, int b) {

	return ((((input[a] == 'A' && input[b] == 'U') || (input[a] == 'U' && input[b] == 'A') || (input[a] == 'G' && input[b] == 'C') ||
			  (input[a] == 'C' && input[b] == 'G') || (input[a] == 'G' && input[b] == 'U') || (input[a] == 'U' && input[b] == 'G')))) /* && (a < b - 1))) */
			   ? 1
			   : 0;
}

// n: dataset size, passed as an argument so that it is a symbolic
// parameter in the scop.
static void kernel_nussinov(int n) {
	int i, j, k;

#pragma scop
	for (i = n - 1; i >= 0; i--) {

		for (j = i + 1; j < n; j++) {
			for (k = 0; k < j - i; k++) {
				S[i][j] = MAX(S[i][k + i] + S[k + i + 1][j], S[i][j]);
			}
			for (k = 0; k < 1; k++) {
				S[i][j] = MAX(S[i][j], S[i + 1][j - 1] + can_pair(RNA, i, j));
			}
		}
	}
#pragma endscop
}

int main(int argc, char *argv[]) {

	S = memd();
	RNA = rna_seq_alloc();

	srand(RSEED);
	rand_seq(RNA, N);

	metrics_start_timer;
	kernel_nussinov(N);
	metrics_stop_timer;
	metrics_print_timer;

#ifdef RESULT_DUMP
	fprintf(stderr, "RNA sequence\n");
	fprintf(stderr, "==================\n");
	rna_seq_print(stderr, RNA);
	fprintf(stderr, "\nArray S\n");
	fprintf(stderr, "==================\n");
	darray_print(stderr, S);
#endif

	rna_seq_free(RNA);
	free_memd(S);

	return 0;
}
