#include <limits.h>
#include <mem.h>
#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>

int **F;
int **m1;
int **m2;

int *W;
char *a;
char *b;

// Similarity score of the elements that constituted the two sequences
int s(char x, char z) { return (x == z) ? 1 : -1; }

int main(int argc, char *argv[]) {
	srand(RSEED);

	int i, j, k;

	// H is the scoring matrix
	F = mem();
	m1 = mem();
	m2 = mem();

	W = calloc(DIM, sizeof(int));

	a = rna_seq_alloc();
	b = rna_seq_alloc();

	// W is the gap alignment
	W[0] = 2;
	for (i = 1; i <= N; i++)
		W[i] = i * W[0];

	// Global alignment: aligning a prefix against the empty sequence
	// costs a gap of that length.
	F[0][0] = 0;
	for (i = 1; i <= N; i++) {
		F[i][0] = -W[i];
		F[0][i] = -W[i];
	}

	// Sequences are 1-indexed in the kernel (a[1..N])
	rand_seq(a + 1, N);
	rand_seq(b + 1, N);

	metrics_start_timer;

#pragma scop
	for (i = 1; i <= N; i++)
		for (j = 1; j <= N; j++) {
			// Block S
			m1[i][j] = INT_MIN;
			for (k = 1; k <= i; k++)
				m1[i][j] = MAX(m1[i][j], F[i - k][j] - W[k]);
			m2[i][j] = INT_MIN;
			for (k = 1; k <= j; k++)
				m2[i][j] = MAX(m2[i][j], F[i][j - k] - W[k]);
			// Unlike Smith-Waterman (sw), no clamp to 0: global alignment
			F[i][j] = MAX(F[i - 1][j - 1] + s(a[i], b[j]), MAX(m1[i][j], m2[i][j]));
		}
#pragma endscop

	metrics_stop_timer;
	metrics_print_timer;

#ifdef RESULT_DUMP

	fprintf(stderr, "Array F\n");
	fprintf(stderr, "===============\n");
	array_print(stderr, F);

#endif

	free_mem(F);
	free_mem(m1);
	free_mem(m2);
	free(W);
	rna_seq_free(a);
	rna_seq_free(b);

	return 0;
}
