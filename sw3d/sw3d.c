#include <limits.h>
#include <math.h>
#include <mem.h>
#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>

int ***H;

int ***m1;
int ***m2;
int ***m3;
int ***m4;
int ***m5;
int ***m6;

int *W;
unsigned char *a;
unsigned char *b;
unsigned char *c;

// Similarity score of the elements that constituted the two sequences
int s(unsigned char x, unsigned char z) { return (x == z) ? 3 : -3; }

int main(int argc, char *argv[]) {
	srand(RSEED);
	int i, j, k, l;

	// H is the scoring matrix
	H = mem3d();

	m1 = mem3d();
	m2 = mem3d();
	m3 = mem3d();
	m4 = mem3d();
	m5 = mem3d();
	m6 = mem3d();

	W = (int *)malloc(DIM * sizeof(int));
	a = rna_seq_alloc();
	b = rna_seq_alloc();
	c = rna_seq_alloc();

	for (i = 0; i <= N; i++) {
		H[i][0][0] = 0;
		H[0][i][0] = 0;
		H[0][0][i] = 0;
	}

	// W is the gap alignment
	W[0] = 2;
	for (i = 1; i <= N; i++)
		W[i] = i * W[0];

	rand_seq(a, N);
	rand_seq(b, N);
	rand_seq(c, N);

#pragma scop
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			for (l = 1; l <= N; l++) {
				// Block S
				m1[i][j][l] = INT_MIN;
				for (k = 1; k <= i; k++) {
					m1[i][j][l] = MAX(m1[i][j][l], H[i - k][j][l] - 2 * W[k]);
				}

				m2[i][j][l] = INT_MIN;
				for (k = 1; k <= j; k++) {
					m2[i][j][l] = MAX(m2[i][j][l], H[i][j - k][l] - 2 * W[k]);
				}

				m3[i][j][l] = INT_MIN;
				for (k = 1; k <= l; k++) {
					m3[i][j][l] = MAX(m3[i][j][l], H[i][j][l - k] - 2 * W[k]);
				}

				m4[i][j][l] = INT_MIN;
				for (k = 1; k <= min(i, j); k++) {
					m4[i][j][l] = MAX(m4[i][j][l], H[i - k][j - k][l] - W[k] + s(a[i], b[j]));
				}

				m5[i][j][l] = INT_MIN;
				for (k = 1; k <= min(j, l); k++) {
					m5[i][j][l] = MAX(m5[i][j][l], H[i][j - k][l - k] - W[k] + s(b[j], c[l]));
				}

				m6[i][j][l] = INT_MIN;
				for (k = 1; k <= min(i, l); k++) {
					m6[i][j][l] = MAX(m6[i][j][l], H[i - k][j][l - k] - W[k] + s(a[i], c[l]));
				}

				H[i][j][l] = MAX(0, MAX(H[i - 1][j - 1][l - 1] + s(a[i], b[j]) + s(a[i], c[l]) + s(b[j], c[l]),
										MAX(m1[i][j][l], MAX(m2[i][j][l], MAX(m3[i][j][l], MAX(m4[i][j][l], MAX(m5[i][j][l], m6[i][j][l])))))));
			}
		}
	}
#pragma endscop

#ifdef RESULT_DUMP

	fprintf(stderr, "Matrix H\n");
	fprintf(stderr, "==================\n");
	array3d_print(stderr, H);
#endif

	free_mem3d(H);
	free_mem3d(m1);
	free_mem3d(m2);
	free_mem3d(m3);
	free_mem3d(m4);
	free_mem3d(m5);
	free_mem3d(m6);

	free(W);

	rna_seq_free(a);
	rna_seq_free(b);
	rna_seq_free(c);

	return 0;
}
