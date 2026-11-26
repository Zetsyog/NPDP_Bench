// Smith-Waterman local alignment of three sequences (see README.md).

#include <mem.h>
#include <metrics.h>

#include <limits.h>

// n: dataset size (symbolic parameter), the three sequences have length n.
// H: scoring cube, m1..m6: best scores ending with each kind of gap.
// W[k]: penalty of a gap of length k.
// sab[i][j], sbc[j][l], sac[i][l]: pairwise similarities of a[i], b[j], c[l]
// (+3 match, -3 mismatch).
// Loop bounds like k <= min(i, j) are written k <= i && k <= j so that they
// contain no function or macro call.
static void kernel_sw3d(int n, int H[NPDP_DIM][NPDP_DIM][NPDP_DIM], int m1[NPDP_DIM][NPDP_DIM][NPDP_DIM], int m2[NPDP_DIM][NPDP_DIM][NPDP_DIM],
						int m3[NPDP_DIM][NPDP_DIM][NPDP_DIM], int m4[NPDP_DIM][NPDP_DIM][NPDP_DIM], int m5[NPDP_DIM][NPDP_DIM][NPDP_DIM],
						int m6[NPDP_DIM][NPDP_DIM][NPDP_DIM], int W[NPDP_DIM], int sab[NPDP_DIM][NPDP_DIM], int sbc[NPDP_DIM][NPDP_DIM],
						int sac[NPDP_DIM][NPDP_DIM]) {
	int i, j, k, l;

#pragma scop
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) {
			for (l = 1; l <= n; l++) {
				// Block S
				m1[i][j][l] = INT_MIN;
				for (k = 1; k <= i; k++) {
					m1[i][j][l] = NPDP_MAX(m1[i][j][l], H[i - k][j][l] - 2 * W[k]);
				}

				m2[i][j][l] = INT_MIN;
				for (k = 1; k <= j; k++) {
					m2[i][j][l] = NPDP_MAX(m2[i][j][l], H[i][j - k][l] - 2 * W[k]);
				}

				m3[i][j][l] = INT_MIN;
				for (k = 1; k <= l; k++) {
					m3[i][j][l] = NPDP_MAX(m3[i][j][l], H[i][j][l - k] - 2 * W[k]);
				}

				m4[i][j][l] = INT_MIN;
				for (k = 1; k <= i && k <= j; k++) {
					m4[i][j][l] = NPDP_MAX(m4[i][j][l], H[i - k][j - k][l] - W[k] + sab[i][j]);
				}

				m5[i][j][l] = INT_MIN;
				for (k = 1; k <= j && k <= l; k++) {
					m5[i][j][l] = NPDP_MAX(m5[i][j][l], H[i][j - k][l - k] - W[k] + sbc[j][l]);
				}

				m6[i][j][l] = INT_MIN;
				for (k = 1; k <= i && k <= l; k++) {
					m6[i][j][l] = NPDP_MAX(m6[i][j][l], H[i - k][j][l - k] - W[k] + sac[i][l]);
				}

				H[i][j][l] = NPDP_MAX(
					0,
					NPDP_MAX(H[i - 1][j - 1][l - 1] + sab[i][j] + sac[i][l] + sbc[j][l],
							 NPDP_MAX(m1[i][j][l], NPDP_MAX(m2[i][j][l], NPDP_MAX(m3[i][j][l], NPDP_MAX(m4[i][j][l], NPDP_MAX(m5[i][j][l], m6[i][j][l])))))));
			}
		}
	}
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	srand(NPDP_RSEED);

	// H is the scoring cube. Local alignment: its three border faces
	// (index 0 on any axis) are 0 (zero-initialized).
	int (*H)[NPDP_DIM][NPDP_DIM] = NPDP_ALLOC_3D(int);
	int (*m1)[NPDP_DIM][NPDP_DIM] = NPDP_ALLOC_3D(int);
	int (*m2)[NPDP_DIM][NPDP_DIM] = NPDP_ALLOC_3D(int);
	int (*m3)[NPDP_DIM][NPDP_DIM] = NPDP_ALLOC_3D(int);
	int (*m4)[NPDP_DIM][NPDP_DIM] = NPDP_ALLOC_3D(int);
	int (*m5)[NPDP_DIM][NPDP_DIM] = NPDP_ALLOC_3D(int);
	int (*m6)[NPDP_DIM][NPDP_DIM] = NPDP_ALLOC_3D(int);
	int (*sab)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int (*sbc)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int (*sac)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int *W = NPDP_ALLOC_1D(int);
	char *a = NPDP_ALLOC_1D(char);
	char *b = NPDP_ALLOC_1D(char);
	char *c = NPDP_ALLOC_1D(char);

	// W is the gap alignment
	W[0] = 2;
	for (int i = 1; i <= N; i++) {
		W[i] = i * W[0];
	}

	// Sequences are 1-indexed in the kernel (a[1..N])
	npdp_rand_seq(a + 1, N);
	npdp_rand_seq(b + 1, N);
	npdp_rand_seq(c + 1, N);
	npdp_build_sim_table(a, b, 3, -3, sab);
	npdp_build_sim_table(b, c, 3, -3, sbc);
	npdp_build_sim_table(a, c, 3, -3, sac);

	NPDP_TIMER_START;
	kernel_sw3d(N, H, m1, m2, m3, m4, m5, m6, W, sab, sbc, sac);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP
	fprintf(stderr, "Matrix H\n");
	fprintf(stderr, "==================\n");
	NPDP_DUMP_3D(stderr, H, N, " %d ", int);
#endif

	free(H);
	free(m1);
	free(m2);
	free(m3);
	free(m4);
	free(m5);
	free(m6);
	free(sab);
	free(sbc);
	free(sac);
	free(W);
	free(a);
	free(b);
	free(c);
	return 0;
}
#endif // NPDP_NO_MAIN
