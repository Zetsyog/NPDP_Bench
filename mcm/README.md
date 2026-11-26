# mcm: matrix chain multiplication

## Problem

Given matrices `A_1 … A_{n-1}`, where `A_i` is `p[i-1] × p[i]`, find the parenthesization of the product that minimizes the number of scalar multiplications (Cormen et al., *Introduction to Algorithms*, section 15.2):

```
m[i][i] = 0
m[i][j] = min_{i ≤ k < j} ( m[i][k] + m[k+1][j] + p[i-1] · p[k] · p[j] )
s[i][j] = the k that reaches the minimum
```

## Kernel

The table is filled by increasing chain length `L` (`j = i + L - 1`):

```c
static void kernel_mcm(int n, int m[NPDP_DIM][NPDP_DIM], int s[NPDP_DIM][NPDP_DIM], int p[NPDP_DIM]) {
#pragma scop
	for (L = 2; L < n; L++)
		for (i = 1; i < n - L + 1; i++) {
			m[i][i + L - 1] = INT_MAX;
			for (k = i; k <= i + L - 1 - 1; k++) {
				s[i][i + L - 1] = (m[i][k] + m[k + 1][i + L - 1] + p[i - 1] * p[k] * p[i + L - 1] < m[i][i + L - 1]) ? k : s[i][i + L - 1];
				m[i][i + L - 1] = NPDP_MIN(m[i][i + L - 1], m[i][k] + m[k + 1][i + L - 1] + p[i - 1] * p[k] * p[i + L - 1]);
			}
		}
#pragma endscop
}
```

| Array | Meaning |
|---|---|
| `m` | minimum cost of `A_i … A_j` (row and column 0 unused) |
| `s` | optimal split point |
| `p` | matrix dimensions (length `n`) |

## Characteristics

- **Loop nest:** 3 loops, 3 statements: the initialization at depth 2, then 2 statements at depth 3.
- **Cost:** O(n³) time, 2 arrays of `(n+10)²` ints.
- **Parameters:** `n`, the length of `p` (`n - 1` matrices).
- **Iteration order:** the loops iterate over diagonals (`L`), not rows. The access functions `m[i][i+L-1]` are therefore skewed, which is a different starting point for schedulers than the row-wise kernels.
- **Dependences:** `m[i][k]` and `m[k+1][j]`, non-uniform.

## Inputs

- `p[i]` is random in 1..20 (`NPDP_RSEED`).
- `m` and `s` start at 0.

## Output (`NPDP_DUMP`)

- The optimal cost `m[1][n-1]`.
- The upper triangle of the split table `s`.

## Differences from the original

- **Split recording:** the original recorded the split through `minsq(a, b, i, j, k)`, a function called inside the scop that wrote `s[i][j]` as a hidden side effect, invisible to dependence analysis. Here the split is recorded by an explicit statement, with exactly the same semantics: `s` is updated before `m`, when the new cost is strictly smaller.
- `s` is now part of the dump.
- Memory leaks are fixed.

## Validation

`tests/known_answers.c` checks the CLRS example: `p = {30, 35, 15, 5, 10, 20, 25}` gives `m[1][6] = 15125`, `s[1][6] = 3` and `m[2][5] = 7125`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest, C. Stein, *Introduction to Algorithms*, section 15.2 "Matrix-chain multiplication".
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`mcm/mmc.c`).
