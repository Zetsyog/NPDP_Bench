# sw: Smith-Waterman local alignment

## Problem

Smith-Waterman local alignment of two sequences `a` and `b` of length `n`, with a general (arbitrary length) gap penalty `W[k]` (here `W[k] = 2k`). Because of the general gap penalty, every cell scans its whole row and column, which makes the algorithm cubic, and the dependences non-uniform:

```
m1[i][j] = max_{1≤k≤i} H[i-k][j] - W[k]                 (gap in b)
m2[i][j] = max_{1≤k≤j} H[i][j-k] - W[k]                 (gap in a)
H[i][j]  = max(0, H[i-1][j-1] + sim(a[i], b[j]), m1[i][j], m2[i][j])
```

`sim` is +1 for a match and -1 for a mismatch. The borders `H[i][0]` and `H[0][j]` are 0. The clamp to 0 makes the alignment local.

## Kernel

```c
static void kernel_sw(int n, int H[NPDP_DIM][NPDP_DIM], int m1[NPDP_DIM][NPDP_DIM], int m2[NPDP_DIM][NPDP_DIM], int W[NPDP_DIM], int sim[NPDP_DIM][NPDP_DIM]) {
#pragma scop
	for (i = 1; i <= n; i++)
		for (j = 1; j <= n; j++) {
			m1[i][j] = INT_MIN;
			for (k = 1; k <= i; k++)
				m1[i][j] = NPDP_MAX(m1[i][j], H[i - k][j] - W[k]);
			m2[i][j] = INT_MIN;
			for (k = 1; k <= j; k++)
				m2[i][j] = NPDP_MAX(m2[i][j], H[i][j - k] - W[k]);
			H[i][j] = NPDP_MAX(0, NPDP_MAX(H[i - 1][j - 1] + sim[i][j], NPDP_MAX(m1[i][j], m2[i][j])));
		}
#pragma endscop
}
```

| Array | Meaning |
|---|---|
| `H` | scoring matrix |
| `m1`, `m2` | best score ending with a vertical or horizontal gap |
| `W[k]` | penalty of a gap of length `k` |
| `sim[i][j]` | similarity of `a[i]` and `b[j]`, built by `npdp_build_sim_table` |

## Characteristics

- **Loop nest:** 3 loops (two sibling `k` loops), 5 statements, at depths 2 and 3.
- **Cost:** O(n³) time, 4 arrays of `(n+10)²` ints.
- **Parameters:** `n`.
- **Dependences:** `H[i-k][j]` (whole column above) and `H[i][j-k]` (whole row to the left) give non-uniform dependences. `H[i-1][j-1]` is a uniform one. Rows and columns are both carried, which leads to wavefront parallelism.

## Inputs

- Two random sequences over {A, C, G, U} (`NPDP_RSEED`), stored at indices `1..n`.
- `W[k] = 2k`.

## Output (`NPDP_DUMP`)

The `n × n` matrix `H`. Row and column `n` are not printed.

## Differences from the original

- **Score index:** the original scored `s(a[i], b[i])` instead of `s(a[i], b[j])`.
- **Gap sign:** upstream's `sw_oryg.h` added the gap penalty (`+ W[k]`), which made gaps a bonus.
- **Sequence indexing:** sequences were filled at `0..n-1` while the kernel reads `1..n`, so `a[n]` was uninitialized.
- **No function call:** `s(a[i], b[j])` is replaced by the `sim` table, so the scop contains no call.

## Validation

`tests/known_answers.c` checks that identical sequences of length 8 give 8, and that `AAAA` vs `UUUU` gives 0.

## References

- T. F. Smith, M. S. Waterman, *Identification of common molecular subsequences*, Journal of Molecular Biology 147(1), 195–197, 1981.
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/sw.cpp`).
