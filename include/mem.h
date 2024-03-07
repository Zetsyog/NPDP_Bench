#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef N
#define N 8
#else
#define DIM (N + 10)
#endif

#ifndef DIM
#define DIM (N + 10)
#endif

#ifndef RSEED
#define RSEED time(NULL)
#endif

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

int **mem() {
	int i;
	int **S;
	S = malloc(DIM * sizeof(int *));

	for (i = 0; i < DIM; i++) {
		S[i] = malloc(DIM * sizeof(int));
	}

	return S;
}

double **memd() {
	int i;
	double **S = malloc(DIM * sizeof(double *));

	for (i = 0; i < DIM; i++) {
		S[i] = malloc(DIM * sizeof(double));
	}

	return S;
}

void free_mem(int **S) {
	int i;
	for (i = 0; i < DIM; i++)
		free(S[i]);
	free(S);
}

void free_memd(double **S) {
	int i;
	for (i = 0; i < DIM; i++)
		free(S[i]);
	free(S);
}

void rand_seq(char *a, int _N) {
	int i, tmp;

	for (i = 0; i < _N; i++) {
		tmp = rand() % 4;

		switch (tmp) {
		case 0:
			a[i] = 'A';
			break;
		case 1:
			a[i] = 'G';
			break;
		case 2:
			a[i] = 'C';
			break;
		case 3:
			a[i] = 'U';
			break;
		}
	}
}

void rna_array_init(double **S, double def, double def2) {

	int i, j;

	for (i = 0; i <= N + 5; i++) {
		for (j = 0; j <= N + 5; j++) {
			if (i == j || i == 0) {
				S[i][j] = def;
			} else {
				S[i][j] = def2;
			}
		}
	}
}

void rna_array_print(FILE *fp, double **S) {
	int i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i > j) {
				fprintf(fp, "       ");
			} else {
				fprintf(fp, " %5.3f ", S[i][j]);
			}
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}

void array_print(FILE *fp, int **S) {
	int i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			fprintf(fp, " %d ", S[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}
