// Minimum cost polygon triangulation (see README.md and Cormen et al.).
// Accelerating Minimum Cost Polygon Triangulation Code with the TRACO Compiler,
// Palkowski, Bielecki, FedCSIS 2018: https://annals-csis.org/Volume_17/drp/pdf/8.pdf

#include <mem.h>
#include <metrics.h>

#include <limits.h>
#include <math.h>

// Minimum cost triangulation of a convex polygon of n vertices.
// n: dataset size (symbolic parameter).
// table[i][j]: minimum cost to triangulate vertices i..j.
// D[x][y]: euclidean distance between vertices x and y. The cost of the
// triangle (i, j, k) is its perimeter D[i][j] + D[j][k] + D[k][i], written
// inline so that the scop contains no function call.
static void kernel_triang(int n, double table[NPDP_DIM][NPDP_DIM], double D[NPDP_DIM][NPDP_DIM]) {
	int j, k, gap;

#pragma scop
	for (gap = 0; gap < n; gap++) {
		for (j = gap; j < n; j++) // i = j - gap
		{
			if (gap < 2)
				table[j - gap][j] = 0.0;
			else {
				table[j - gap][j] = INT_MAX;
				for (k = j - gap + 1; k < j; k++) {
					table[j - gap][j] = NPDP_MIN(table[j - gap][j], table[j - gap][k] + table[k][j] + (D[j - gap][j] + D[j][k] + D[k][j - gap]));
				}
			}
		}
	}
#pragma endscop
}

#ifndef NPDP_NO_MAIN
int main(void) {
	int(*points)[2] = npdp_xcalloc(NPDP_DIM, sizeof *points);
	double(*table)[NPDP_DIM] = NPDP_ALLOC_2D(double);
	double(*D)[NPDP_DIM] = NPDP_ALLOC_2D(double);

	// Random points
	srand(NPDP_RSEED);
	for (int i = 0; i < NPDP_DIM; i++) {
		points[i][0] = rand() % 100;
		points[i][1] = rand() % 100;
	}

	for (int x = 0; x < NPDP_DIM; x++)
		for (int y = 0; y < NPDP_DIM; y++) {
			int dx = points[x][0] - points[y][0];
			int dy = points[x][1] - points[y][1];
			D[x][y] = sqrt(dx * dx + dy * dy);
		}

	NPDP_TIMER_START;
	kernel_triang(N, table, D);
	NPDP_TIMER_STOP;
	NPDP_TIMER_PRINT;

#ifdef NPDP_DUMP
	fprintf(stderr, "%f\n", N > 0 ? table[0][N - 1] : 0.0);
#endif

	free(points);
	free(table);
	free(D);
	return 0;
}
#endif // NPDP_NO_MAIN
