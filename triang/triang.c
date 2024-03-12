
// See the Cormen book for details of the following algorithm

// Accelerating Minimum Cost Polygon Triangulation Code with the TRACO Compiler, Palkowski Bielecki FedCsis 2018
// https://annals-csis.org/Volume_17/drp/pdf/8.pdf

#include <limits.h>
#include <math.h>
#include <mem.h>
#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>

int **points;
double **table;

// A utility function to find distance between two points in a plane
double dist(int *p1, int *p2) { return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1])); }

// A utility function to find cost of a triangle. The cost is considered
// as perimeter (sum of lengths of all edges) of the triangle
double cost(int i, int j, int k) {
	int *p1 = points[i];
	int *p2 = points[j];
	int *p3 = points[k];
	return dist(p1, p2) + dist(p2, p3) + dist(p3, p1);
}

// Matrix Ai has dimension p[i-1] x p[i] for i = 1..n
double mcTDP() {

	table = memd();

	int i, j, k, gap, q;

	metrics_start_timer;

#pragma scop
	for (gap = 0; gap < N; gap++) {
		for (j = gap; j < N; j++) // i = j - gap
		{
			if (gap < 2)
				table[j - gap][j] = 0.0;
			else {
				table[j - gap][j] = INT_MAX;
				for (k = j - gap + 1; k < j; k++) {
					table[j - gap][j] = MIN(table[j - gap][j], table[j - gap][k] + table[k][j] + cost(j - gap, j, k));
				}
			}
		}
	}
#pragma endscop

	metrics_stop_timer;
	metrics_print_timer;

	return table[0][N - 1];
}

int main(int argc, char *argv[]) {
	points = malloc(DIM * sizeof(int *));

	for (int i = 0; i < DIM; i++)
		points[i] = malloc(2 * sizeof(int));

	// Random points
	srand(RSEED);
	for (int i = 0; i < DIM; i++) {
		points[i][0] = rand() % 100;
		points[i][1] = rand() % 100;
	}

	double res = mcTDP();

#ifdef RESULT_DUMP
	fprintf(stderr, "%f\n", res);
#endif

	for (int i = 0; i < DIM; i++)
		free(points[i]);
	free(points);

	free_memd(table);

	return 0;
}
