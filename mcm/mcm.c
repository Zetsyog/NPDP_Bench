// See the Cormen book for details of the following algorithm
#include <limits.h>
#include <math.h>
#include <mem.h>
#include <metrics.h>
#include <stdio.h>

int **s;

int minsq(int a, int b, int i, int j, int k) {
	if (a < b) {
		s[i][j] = k;
		return a;
	} else
		return b;
}

// Matrix Ai has dimension p[i-1] x p[i] for i = 1..n
int MatrixChainOrder(int p[], int n) {

	/* For simplicity of the program, one extra row and one
	   extra column are allocated in m[][].  0th row and 0th
	   column of m[][] are not used */
	int **m = mem();
	s = mem();

	int i, j, k, L, q;

	/* m[i,j] = Minimum number of scalar multiplications needed
	   to compute the matrix A[i]A[i+1]...A[j] = A[i..j] where
	   dimension of A[i] is p[i-1] x p[i] */

	// cost is zero when multiplying one matrix.
	for (i = 1; i < n; i++)
		m[i][i] = 0;

	metrics_start_timer;

#pragma scop
	// L is chain length.
	for (L = 2; L < n; L++) {
		for (i = 1; i < n - L + 1; i++) {
			m[i][i + L - 1] = INT_MAX;
			for (k = i; k <= i + L - 1 - 1; k++) {
				m[i][i + L - 1] = minsq(m[i][k] + m[k + 1][i + L - 1] + p[i - 1] * p[k] * p[i + L - 1], m[i][i + L - 1], i, i + L - 1, k);
			}
		}
	}
#pragma endscop

	metrics_stop_timer;
	metrics_print_timer;

	return m[1][n - 1];
}

int main(int argc, char *argv[]) {
	srand(RSEED);
	int *p = calloc(DIM, sizeof(int));
	for (int i = 0; i < N; i++) {
		// TODO: max matrix size as macro ?
		p[i] = (rand() % 20) + 1;
	}

	// int p[] = {1, 2, 3, 4};
	// N = sizeof(p)/sizeof(p[0]);

	fprintf(stderr, "Minimum number of multiplications is %d \n", MatrixChainOrder(p, N));

	return 0;
}
