// use knut_fix.cpp

#include <mem.h>
#include <metrics.h>
#include <stdlib.h>

int **c;
int **w;

// n: dataset size, passed as an argument so that it is a symbolic
// parameter in the scop.
static void kernel_knuth(int n) {
	int i, j, k;

#pragma scop
	for (i = n - 1; i >= 1; i--)
		for (j = i + 1; j <= n; j++)
			for (k = i + 1; k < j; k++)
				c[i][j] = MIN(c[i][j], w[i][j] + c[i][k] + c[k][j]);
#pragma endscop
}

int main(int argc, char *argv[]) {

	int i, j;

	c = mem();
	w = mem();

	for (i = 0; i < DIM; i++)
		for (j = 0; j < DIM; j++) {
			c[i][j] = i + j;
			w[i][j] = i - j;
		}

	metrics_start_timer;
	kernel_knuth(N);
	metrics_stop_timer;
	metrics_print_timer;

#ifdef RESULT_DUMP
	fprintf(stderr, "DUMPING RESULTS\n");
	fprintf(stderr, "Array Q\n");
	fprintf(stderr, "===============\n");
	array_print(stderr, c);
#endif
	free_mem(c);
	free_mem(w);

	return 0;
}
