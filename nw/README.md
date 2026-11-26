# nw: Needleman-Wunsch global alignment

## Problem

Needleman-Wunsch global alignment of two sequences of length `n`, with a general gap penalty `W[k] = 2k`. It is the same recurrence as [sw](../sw/README.md) without the clamp to 0. The borders hold the cost of aligning a prefix against the empty sequence:

```
F[0][0] = 0,  F[i][0] = F[0][i] = -W[i]
m1[i][j] = max_{1≤k≤i} F[i-k][j] - W[k]
m2[i][j] = max_{1≤k≤j} F[i][j-k] - W[k]
F[i][j]  = max(F[i-1][j-1] + sim(a[i], b[j]), m1[i][j], m2[i][j])
```

## Kernel

```c
static void kernel_nw(int n, int F[NPDP_DIM][NPDP_DIM], int m1[NPDP_DIM][NPDP_DIM], int m2[NPDP_DIM][NPDP_DIM], int W[NPDP_DIM], int sim[NPDP_DIM][NPDP_DIM]) {
#pragma scop
	for (i = 1; i <= n; i++)
		for (j = 1; j <= n; j++) {
			m1[i][j] = INT_MIN;
			for (k = 1; k <= i; k++)
				m1[i][j] = NPDP_MAX(m1[i][j], F[i - k][j] - W[k]);
			m2[i][j] = INT_MIN;
			for (k = 1; k <= j; k++)
				m2[i][j] = NPDP_MAX(m2[i][j], F[i][j - k] - W[k]);
			F[i][j] = NPDP_MAX(F[i - 1][j - 1] + sim[i][j], NPDP_MAX(m1[i][j], m2[i][j]));
		}
#pragma endscop
}
```

The arrays and the iteration space are the same as in [sw](../sw/README.md).

## Characteristics

The loop structure, dependences and cost are the same as sw: 3 loops, 5 statements, O(n³) time and 4 arrays of `(n+10)²` ints. Only the values differ (negative scores, no clamp), so the two kernels cover the same iteration space with different data.

## Inputs

- Two random sequences over {A, C, G, U} (`NPDP_RSEED`), stored at indices `1..n`.
- `W[k] = 2k`.
- Borders set to `-W[i]`.

## Output (`NPDP_DUMP`)

The `n × n` matrix `F`.

## Differences from the original

- **Algorithm:** upstream's `nw` was a copy of sw: it kept the clamp to 0 and had zero borders, so it computed a local alignment. Its `nw.c` did not compile. This version is a real global alignment.
- The same fixes as sw: `b[j]` instead of `b[i]`, 1-based sequences, and the `sim` table instead of a function call.

## Validation

`tests/known_answers.c` checks that identical sequences of length 8 give 8, and that `AAAA` vs `UUUU` gives -4 (four mismatches cost less than any gap).

## References

- S. B. Needleman, C. D. Wunsch, *A general method applicable to the search for similarities in the amino acid sequence of two proteins*, Journal of Molecular Biology 48(3), 443–453, 1970.
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/nw.cpp`).
