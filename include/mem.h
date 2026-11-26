// Shared helpers for the NPDP benchmarks.
//
// Every shared macro and function is prefixed with NPDP_ / npdp_.
// Configuration macros (set with -D at compile time):
//   N            dataset size (default 8)
//   NPDP_DIM     allocated extent of each array dimension (default N + 10)
//   NPDP_RSEED   seed for the random input generators (default time(NULL))
//   NPDP_DUMP    dump the results on stderr (see each benchmark)
//   NPDP_TIME    print the kernel execution time on stdout (see metrics.h)
//
// Conventions for polyhedral compilers:
//   - Kernels never use N or NPDP_DIM inside a scop: main() passes the size
//     to kernel_<name>(int n, ...) so that it is a symbolic parameter, even
//     for tools that run the preprocessor first.
//   - Arrays are contiguous, statically sized (T A[NPDP_DIM][NPDP_DIM]) and
//     passed as kernel parameters, so that every access is affine.
//   - Scops contain no function calls: per-element functions are replaced by
//     lookup tables filled before the timed region (npdp_build_sim_table, ...).
//     Exception: base pairing in the RNA kernels is a function call by default
//     (as in the original codes); -DNPDP_PAIR_TABLE replaces it by a table
//     (see NPDP_PAIRED below).
//   NPDP_PAIR_TABLE  use a lookup table for base pairing (RNA kernels)
#ifndef NPDP_MEM_H
#define NPDP_MEM_H

// clock_gettime (metrics.h) needs POSIX; this header is included first.
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef N
#define N 8
#endif

#ifndef NPDP_DIM
#define NPDP_DIM (N + 10)
#endif

#ifndef NPDP_RSEED
#define NPDP_RSEED time(NULL)
#endif

#define NPDP_MIN(x, y) ((x) < (y) ? (x) : (y))
#define NPDP_MAX(x, y) ((x) > (y) ? (x) : (y))

#ifndef NPDP_FP_TYPE
#define NPDP_FP_TYPE double
#endif

#ifndef NPDP_INT_TYPE
#define NPDP_INT_TYPE int
#endif

// ---------------------------------------------------------------------------
// Allocation: zero-initialized, contiguous arrays. Release with free().
// ---------------------------------------------------------------------------

static inline void *npdp_xcalloc(size_t nmemb, size_t size) {
	void *p = calloc(nmemb, size);
	if (p == NULL) {
		fprintf(stderr, "npdp: cannot allocate %zu x %zu bytes\n", nmemb, size);
		exit(EXIT_FAILURE);
	}
	return p;
}

// T *          with NPDP_DIM elements
#define NPDP_ALLOC_1D(T) ((T *)npdp_xcalloc((size_t)NPDP_DIM, sizeof(T)))
// T (*)[NPDP_DIM]   usable as A[i][j]
#define NPDP_ALLOC_2D(T) ((T(*)[NPDP_DIM])npdp_xcalloc((size_t)NPDP_DIM * NPDP_DIM, sizeof(T)))
// T (*)[NPDP_DIM][NPDP_DIM] usable as A[i][j][k]
#define NPDP_ALLOC_3D(T) ((T(*)[NPDP_DIM][NPDP_DIM])npdp_xcalloc((size_t)NPDP_DIM * NPDP_DIM * NPDP_DIM, sizeof(T)))

// ---------------------------------------------------------------------------
// Result dumps. fmt is a printf format for one element after the cast.
// ---------------------------------------------------------------------------

// Full n x n matrix.
#define NPDP_DUMP_2D(fp, A, n, fmt, cast)                                                                                                                      \
	do {                                                                                                                                                       \
		for (int npdp_i_ = 0; npdp_i_ < (n); npdp_i_++) {                                                                                                      \
			for (int npdp_j_ = 0; npdp_j_ < (n); npdp_j_++)                                                                                                    \
				fprintf((fp), (fmt), (cast)(A)[npdp_i_][npdp_j_]);                                                                                             \
			fprintf((fp), "\n");                                                                                                                               \
		}                                                                                                                                                      \
		fprintf((fp), "\n");                                                                                                                                   \
	} while (0)

// Upper triangle of an n x n matrix (lower triangle printed as blanks).
#define NPDP_DUMP_UPPER_2D(fp, A, n, fmt, cast)                                                                                                                \
	do {                                                                                                                                                       \
		for (int npdp_i_ = 0; npdp_i_ < (n); npdp_i_++) {                                                                                                      \
			for (int npdp_j_ = 0; npdp_j_ < (n); npdp_j_++) {                                                                                                  \
				if (npdp_i_ > npdp_j_)                                                                                                                         \
					fprintf((fp), "       ");                                                                                                                  \
				else                                                                                                                                           \
					fprintf((fp), (fmt), (cast)(A)[npdp_i_][npdp_j_]);                                                                                         \
			}                                                                                                                                                  \
			fprintf((fp), "\n");                                                                                                                               \
		}                                                                                                                                                      \
		fprintf((fp), "\n");                                                                                                                                   \
	} while (0)

// n x n x n array, one n x n slice per block.
#define NPDP_DUMP_3D(fp, A, n, fmt, cast)                                                                                                                      \
	do {                                                                                                                                                       \
		for (int npdp_i_ = 0; npdp_i_ < (n); npdp_i_++) {                                                                                                      \
			for (int npdp_j_ = 0; npdp_j_ < (n); npdp_j_++) {                                                                                                  \
				for (int npdp_k_ = 0; npdp_k_ < (n); npdp_k_++)                                                                                                \
					fprintf((fp), (fmt), (cast)(A)[npdp_i_][npdp_j_][npdp_k_]);                                                                                \
				fprintf((fp), "\n");                                                                                                                           \
			}                                                                                                                                                  \
			fprintf((fp), "\n");                                                                                                                               \
		}                                                                                                                                                      \
	} while (0)

// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------

// Fill a[0..n-1] with random nucleotides from {A, G, C, U}.
static inline void npdp_rand_seq(char *a, int n) {
	static const char nt[4] = {'A', 'G', 'C', 'U'};
	for (int i = 0; i < n; i++) {
		a[i] = nt[rand() % 4];
	}
}

// Print the first n characters of a sequence.
static inline void npdp_dump_seq(FILE *fp, const char *a, int n) {
	for (int i = 0; i < n; i++) {
		fprintf(fp, " %c ", a[i]);
	}
	fprintf(fp, "\n");
}

// RNA matrices: S[i][j] = def on the diagonal and first row, def2 elsewhere.
static inline void npdp_rna_array_init(double S[][NPDP_DIM], double def, double def2) {
	for (int i = 0; i <= N + 5; i++) {
		for (int j = 0; j <= N + 5; j++) {
			S[i][j] = (i == j || i == 0) ? def : def2;
		}
	}
}

// 1 if nucleotides x and y can form a base pair (Watson-Crick or G-U wobble).
static inline int npdp_can_pair(char x, char y) {
	return (x == 'A' && y == 'U') || (x == 'U' && y == 'A') || (x == 'G' && y == 'C') || (x == 'C' && y == 'G') || (x == 'G' && y == 'U') ||
		   (x == 'U' && y == 'G');
}

// pair[x][y] = npdp_can_pair(seq[x], seq[y]) for the whole NPDP_DIM extent.
// Replaces a function call in the kernels by an affine array read.
static inline void npdp_build_pair_table(const char *seq, int pair[][NPDP_DIM]) {
	for (int x = 0; x < NPDP_DIM; x++) {
		for (int y = 0; y < NPDP_DIM; y++) {
			pair[x][y] = npdp_can_pair(seq[x], seq[y]);
		}
	}
}

// 1 if nucleotides seq[x] and seq[y] can pair: the paired()/can_pair()
// function of the original codes.
static inline int npdp_paired(const char *seq, int x, int y) { return npdp_can_pair(seq[x], seq[y]); }

// Base pairing in the RNA kernels (counting, counting_im, mcc, mea, nussinov).
// In a kernel:  static void kernel_x(int n, ..., NPDP_PAIR_PARAM)
//               ... NPDP_PAIRED(k, j) ...            (inside the scop)
// In main:      npdp_pair_t pair = npdp_pair_init(rna);
//               kernel_x(N, ..., pair);
//               npdp_pair_free(pair);
//
// Default: NPDP_PAIRED(x, y) calls npdp_paired(seq, x, y), as the original
// codes did.
// -DNPDP_PAIR_TABLE: NPDP_PAIRED(x, y) reads pair[x][y], a table built before
// the timed region, so that the scop contains no function call once
// preprocessed (Pluto parses the unpreprocessed source and sees a call-like
// NPDP_PAIRED(x, y) in both modes). Results are identical in both modes.
#ifdef NPDP_PAIR_TABLE
#define NPDP_PAIR_PARAM int pair[NPDP_DIM][NPDP_DIM]
#define NPDP_PAIRED(x, y) pair[x][y]
typedef int (*npdp_pair_t)[NPDP_DIM];
static inline npdp_pair_t npdp_pair_init(const char *seq) {
	npdp_pair_t pair = NPDP_ALLOC_2D(int);
	npdp_build_pair_table(seq, pair);
	return pair;
}
static inline void npdp_pair_free(npdp_pair_t pair) { free(pair); }
#else
#define NPDP_PAIR_PARAM const char seq[NPDP_DIM]
#define NPDP_PAIRED(x, y) npdp_paired(seq, x, y)
typedef const char *npdp_pair_t;
static inline npdp_pair_t npdp_pair_init(const char *seq) { return seq; }
static inline void npdp_pair_free(npdp_pair_t pair) { (void)pair; }
#endif

// sim[x][y] = match if a[x] == b[y], mismatch otherwise (alignment scores).
static inline void npdp_build_sim_table(const char *a, const char *b, int match, int mismatch, int sim[][NPDP_DIM]) {
	for (int x = 0; x < NPDP_DIM; x++) {
		for (int y = 0; y < NPDP_DIM; y++) {
			sim[x][y] = (a[x] == b[y]) ? match : mismatch;
		}
	}
}

#endif // NPDP_MEM_H
