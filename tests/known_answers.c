// Known-answer tests: run a benchmark kernel on a small hand-checkable input.
//
// Compiled once per benchmark with -DTEST_<BENCH> (see run_tests.sh), e.g.
//   gcc -I include -DN=16 -DNPDP_NO_MAIN -DTEST_MCM tests/known_answers.c -lm
// The benchmark source is included with NPDP_NO_MAIN so that its static
// kernel function can be called directly.

#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, ...)                                                                                                                                       \
	do {                                                                                                                                                       \
		if (!(cond)) {                                                                                                                                         \
			fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__);                                                                                               \
			fprintf(stderr, __VA_ARGS__);                                                                                                                      \
			fprintf(stderr, "\n");                                                                                                                             \
			failures++;                                                                                                                                        \
		}                                                                                                                                                      \
	} while (0)

#if defined(TEST_MCM)
#include "../mcm/mcm.c"

// CLRS (Cormen et al.), section 15.2: 6 matrices of dimensions
// 30x35, 35x15, 15x5, 5x10, 10x20, 20x25 -> 15125 multiplications,
// optimal parenthesization ((A1 (A2 A3)) ((A4 A5) A6)), i.e. s[1][6] = 3.
static void run(void) {
	static const int dims[] = {30, 35, 15, 5, 10, 20, 25};
	int n = 7;
	int *p = NPDP_ALLOC_1D(int);
	int(*m)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*s)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	memcpy(p, dims, sizeof dims);

	kernel_mcm(n, m, s, p);
	CHECK(m[1][6] == 15125, "mcm: m[1][6] = %d, expected 15125", m[1][6]);
	CHECK(s[1][6] == 3, "mcm: s[1][6] = %d, expected 3", s[1][6]);
	CHECK(m[2][5] == 7125, "mcm: m[2][5] = %d, expected 7125", m[2][5]);

	free(p);
	free(m);
	free(s);
}

#elif defined(TEST_NUSSINOV)
#include "../nussinov/nussinov.c"

static NPDP_FP_TYPE fold(const char *seq) {
	int n = (int)strlen(seq);
	char *rna = NPDP_ALLOC_1D(char);
	NPDP_FP_TYPE(*S)[NPDP_DIM] = NPDP_ALLOC_2D(NPDP_FP_TYPE);
	memcpy(rna, seq, (size_t)n);
	npdp_pair_t pair = npdp_pair_init(rna); // table only with -DNPDP_PAIR_TABLE
	kernel_nussinov(n, S, pair);
	NPDP_FP_TYPE res = S[0][n - 1];
	npdp_pair_free(pair);
	free(rna);
	free(S);
	return res;
}

static void run(void) {
	// three G-C pairs can be nested
	CHECK(fold("GGGAAACCC") == 3, "nussinov: GGGAAACCC = %Lf, expected 3", (long double)fold("GGGAAACCC"));
	// A cannot pair with A
	CHECK(fold("AAAA") == 0, "nussinov: AAAA = %Lf, expected 0", (long double)fold("AAAA"));
	// G-U wobble pairs count too
	CHECK(fold("GGUU") == 2, "nussinov: GGUU = %Lf, expected 2", (long double)fold("GGUU"));
}

#elif defined(TEST_SW) || defined(TEST_NW)
#if defined(TEST_SW)
#include "../sw/sw.c"
#define ALIGN kernel_sw
#else
#include "../nw/nw.c"
#define ALIGN kernel_nw
#endif

// Align a and b (same length), return the bottom-right score.
static int align(const char *sa, const char *sb) {
	int n = (int)strlen(sa);
	int(*F)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*m1)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*m2)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int(*sim)[NPDP_DIM] = NPDP_ALLOC_2D(int);
	int *W = NPDP_ALLOC_1D(int);
	char *a = NPDP_ALLOC_1D(char);
	char *b = NPDP_ALLOC_1D(char);
	W[0] = 2;
	for (int i = 1; i <= n; i++)
		W[i] = i * W[0];
#if defined(TEST_NW)
	for (int i = 1; i <= n; i++) {
		F[i][0] = -W[i];
		F[0][i] = -W[i];
	}
#endif
	memcpy(a + 1, sa, (size_t)n);
	memcpy(b + 1, sb, (size_t)n);
	npdp_build_sim_table(a, b, 1, -1, sim);
	ALIGN(n, F, m1, m2, W, sim);
	int res = F[n][n];
	free(F);
	free(m1);
	free(m2);
	free(sim);
	free(W);
	free(a);
	free(b);
	return res;
}

static void run(void) {
	// identical sequences: n matches
	CHECK(align("ACGUACGU", "ACGUACGU") == 8, "align: identical = %d, expected 8", align("ACGUACGU", "ACGUACGU"));
#if defined(TEST_NW)
	// global alignment of fully different sequences: 4 mismatches (-4) beat any gap (W[k] = 2k)
	CHECK(align("AAAA", "UUUU") == -4, "nw: AAAA/UUUU = %d, expected -4", align("AAAA", "UUUU"));
#else
	// local alignment never goes below 0
	CHECK(align("AAAA", "UUUU") == 0, "sw: AAAA/UUUU = %d, expected 0", align("AAAA", "UUUU"));
#endif
}

#elif defined(TEST_TRIANG)
#include "../triang/triang.c"

// Square of side 10: two triangles of perimeter 10 + 10 + 10*sqrt(2).
static void run(void) {
	static const int pts[4][2] = {{0, 0}, {10, 0}, {10, 10}, {0, 10}};
	double(*table)[NPDP_DIM] = NPDP_ALLOC_2D(double);
	double(*D)[NPDP_DIM] = NPDP_ALLOC_2D(double);
	for (int x = 0; x < 4; x++)
		for (int y = 0; y < 4; y++) {
			int dx = pts[x][0] - pts[y][0];
			int dy = pts[x][1] - pts[y][1];
			D[x][y] = sqrt(dx * dx + dy * dy);
		}
	kernel_triang(4, table, D);
	double expected = 40.0 + 20.0 * sqrt(2.0);
	CHECK(fabs(table[0][3] - expected) < 1e-9, "triang: square = %.12f, expected %.12f", table[0][3], expected);
	free(table);
	free(D);
}

#else
#error "define one of TEST_MCM, TEST_NUSSINOV, TEST_SW, TEST_NW, TEST_TRIANG"
#endif

int main(void) {
	run();
	if (failures == 0)
		printf("ok\n");
	return failures == 0 ? 0 : 1;
}
