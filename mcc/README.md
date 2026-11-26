# mcc: McCaskill RNA partition function

## Problem

McCaskill's algorithm computes the partition function of an RNA sequence: the sum, over all secondary structures, of their Boltzmann weights. From it, base-pair probabilities can be derived. This kernel is the simplified core recurrence used in the NPDP suite. Each base pair has the same weight, `ert = exp(-Ebp / RT)`:

```
Qbp[k][j] = Q[k+1][j-1] · ert · pair(k, j-1)        (k and j-1 paired)
Q[i][j]   = Q[i][j-1] + Σ_{k=i}^{j-l-1} Q[i][k] · Qbp[k][j]
```

`Q[i][j]` is the partition function of the subsequence `i..j` and `Qbp[k][j]` the contribution of structures closed by a pair. The initial values (`Q = 0.4`, `Qbp = 0.5`) are synthetic.

## Kernel

```c
static void kernel_mcc(int n, int l, double ert, double Q[NPDP_DIM][NPDP_DIM], double Qbp[NPDP_DIM][NPDP_DIM], NPDP_PAIR_PARAM) {
#pragma scop
	if (n >= 1 && l >= 0 && l <= 5) {
		for (i = n - 1; i >= 0; i--)
			for (j = i + 1; j < n; j++) {
				Q[i][j] = Q[i][j - 1];
				for (k = 0; k < j - i - l; k++) {
					Qbp[k + i][j] = Q[k + i + 1][j - 1] * ert * NPDP_PAIRED(k + i, j - 1);
					Q[i][j] += Q[i][k + i] * Qbp[k + i][j];
				}
			}
	}
#pragma endscop
}
```

| Array | Meaning |
|---|---|
| `Q` | partition function |
| `Qbp` | contribution of structures closed by a base pair |
| `NPDP_PAIRED(x, y)` | 1 if nucleotides `x` and `y` can pair (A-U, G-C, G-U): an `npdp_paired(seq, x, y)` call by default, as in the original, or a read of the precomputed table `pair[x][y]` with `-DNPDP_PAIR_TABLE` |

## Characteristics

- **Loop nest:** 3 loops, 3 statements (1 at depth 2, 2 at depth 3), inside a guard on the parameters.
- **Cost:** O(n³) time, 2 `double` arrays of `(n+10)²` elements, plus the `int` pair table with `-DNPDP_PAIR_TABLE`.
- **Parameters:**
  - `n`;
  - `l`, the minimum hairpin loop length (`NPDP_MIN_LOOP`, 0–5). A pair `(x, y)` is considered only if `y - x ≥ l`;
  - `ert`, a read-only scalar.
- **Dependences:** two statements in the `k` loop, linked by `Qbp[k+i][j]` (produced and immediately consumed). `Q[i][j]` depends non-uniformly on `Q[i][k+i]` and `Q[k+i+1][j-1]`.
- **Floating point:** `Q[i][j]` accumulates a sum. Compile with `-ffp-contract=off` to compare results exactly across loop orders.

## Inputs and macros

- A random sequence over {A, C, G, U} (`NPDP_RSEED`).
- `NPDP_MIN_LOOP` (default 0).
- `NPDP_EBP` (default 0) and `NPDP_RT` (default 1). `ert` is computed in single precision, as in the original.

## Output (`NPDP_DUMP`)

The upper triangles of `Q` and `Qbp` (`%5.3f`).

## Differences from the original

- **Base pairing:** the original function call is kept by default (`NPDP_PAIRED(x, y)` expands to `npdp_paired(seq, x, y)`), for backward compatibility. With `-DNPDP_PAIR_TABLE`, it reads a lookup table `pair[x][y]` built before the timer, so the preprocessed scop contains no call. Results are identical in both modes.
- **Parameters:** `ERT` and `l` are kernel parameters instead of a global and a macro.
- **Macro names:** `Ebp`, `RT` and `l` are renamed `NPDP_EBP`, `NPDP_RT` and `NPDP_MIN_LOOP`.

## References

- J. S. McCaskill, *The equilibrium partition function and base pair binding probabilities for RNA secondary structure*, Biopolymers 29(6-7), 1105–1119, 1990.
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/mcc.cpp`).
