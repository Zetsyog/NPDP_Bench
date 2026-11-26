# nussinov: Nussinov RNA folding

## Problem

The Nussinov algorithm predicts an RNA secondary structure by maximizing the number of base pairs. Pairs must be nested (no pseudoknots). `S[i][j]` is the maximum number of pairs in the subsequence `i..j`:

```
S[i][j] = max( max_{i ≤ k < j} S[i][k] + S[k+1][j],      (bifurcation)
               S[i+1][j-1] + pair(i, j) )                 (i pairs with j)
```

No minimum hairpin loop length is enforced, so adjacent nucleotides may pair.

## Kernel

```c
static void kernel_nussinov(int n, NPDP_FP_TYPE S[NPDP_DIM][NPDP_DIM], NPDP_PAIR_PARAM) {
#pragma scop
	for (i = n - 1; i >= 0; i--)
		for (j = i + 1; j < n; j++) {
			for (k = 0; k < j - i; k++)
				S[i][j] = NPDP_MAX(S[i][k + i] + S[k + i + 1][j], S[i][j]);
			for (k = 0; k < 1; k++)
				S[i][j] = NPDP_MAX(S[i][j], S[i + 1][j - 1] + NPDP_PAIRED(i, j));
		}
#pragma endscop
}
```

| Array | Meaning |
|---|---|
| `S` | DP table, `long double` (as in the original NPDP code) |
| `NPDP_PAIRED(x, y)` | 1 if nucleotides `x` and `y` can pair (A-U, G-C, G-U): an `npdp_paired(seq, x, y)` call by default, as in the original, or a read of the precomputed table `pair[x][y]` with `-DNPDP_PAIR_TABLE` |

## Characteristics

- **Loop nest:** 3 loops, 2 statements. The second statement is in a single-iteration loop (`k < 1`) kept from the original, so both statements sit at depth 3.
- **Cost:** O(n³) time. Memory is `(n+10)²` `long double`, plus the `int` pair table with `-DNPDP_PAIR_TABLE`.
- **Parameters:** `n`.
- **Dependences:** `S[i][k+i]` (same row) and `S[k+i+1][j]` (same column), non-uniform. This is the most studied NPDP kernel in the polyhedral tiling literature.

## Inputs

A random sequence over {A, C, G, U} (`NPDP_RSEED`). `S` starts at 0.

## Output (`NPDP_DUMP`)

The sequence, then the `n × n` matrix `S` (`%.5f`).

## Differences from the original

- **Initialization:** upstream never filled the RNA sequence (the random generator was commented out), so the kernel read uninitialized memory. A random sequence is used here.
- **Dump format:** the dump printed `long double` values with `%f`, which is undefined behavior and printed zeros. Values are now cast to `double`.
- **Base pairing:** the original function call is kept by default (`NPDP_PAIRED(x, y)` expands to `npdp_paired(seq, x, y)`), for backward compatibility. With `-DNPDP_PAIR_TABLE`, it reads a lookup table `pair[x][y]` built before the timer, so the preprocessed scop contains no call. Results are identical in both modes. (The original function was `can_pair(RNA, i, j)`.)

## Validation

`tests/known_answers.c` checks three cases: `GGGAAACCC` gives 3, `AAAA` gives 0 and `GGUU` gives 2 (G-U wobble pairs).

## References

- R. Nussinov, G. Pieczenik, J. R. Griggs, D. J. Kleitman, *Algorithms for loop matchings*, SIAM Journal on Applied Mathematics 35(1), 68–82, 1978.
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/nussinov.cpp`).
