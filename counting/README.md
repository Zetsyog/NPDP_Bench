# counting: RNA structure counting recurrence

## Problem

This kernel is derived from the classic dynamic program that counts the secondary structures of an RNA sequence. In that program, the count for a subsequence `i..j` is obtained by trying every position `k` that can pair with `j`.

This version keeps the dependence pattern of that recurrence but not its arithmetic. It sums instead of multiplying, it starts from synthetic values (`c[i][j] = i + j`), and it has no `c[i][j-1]` term:

```
c[i][j] += Σ_{k = i}^{j - l} [pair(k, j)] · (c[i][k-1] + c[k+1][j-1])
```

The results are therefore not structure counts; the kernel is a benchmark of the loop structure. [counting_im](../counting_im/README.md) is the "corrected" variant, which adds `c[i][j-1]`.

## Kernel

```c
static void kernel_counting(int n, int l, NPDP_INT_TYPE c[NPDP_DIM][NPDP_DIM], NPDP_PAIR_PARAM) {
#pragma scop
	for (i = n - 2; i >= 1; i--)
		for (j = i + 2; j <= n; j++)
			for (k = i; k <= j - l; k++)
				c[i][j] += NPDP_PAIRED(k, j) ? c[i][k - 1] + c[k + 1][j - 1] : 0;
#pragma endscop
}
```

| Array | Meaning |
|---|---|
| `c` | DP table, `unsigned int` (values grow exponentially; unsigned overflow wraps, which is well defined) |
| `NPDP_PAIRED(x, y)` | 1 if nucleotides `x` and `y` can pair (A-U, G-C, G-U): an `npdp_paired(seq, x, y)` call by default, as in the original, or a read of the precomputed table `pair[x][y]` with `-DNPDP_PAIR_TABLE` |

## Characteristics

- **Loop nest:** 3 loops, 1 statement.
- **Cost:** O(n³) time, 1 array of `(n+10)²` elements, plus the `int` pair table with `-DNPDP_PAIR_TABLE`.
- **Parameters:** `n`, and `l`, the minimum loop length (always 0 in `main`).
- **Dependences:** `c[i][j]` reads `c[i][k-1]` (same row, left) and `c[k+1][j-1]` (a lower row). The distance depends on `k`, so the dependences are non-uniform. Rows are computed bottom-up (`i` decreasing).

## Inputs

- A random sequence over {A, C, G, U}, seeded with `NPDP_RSEED`, filled over the whole `NPDP_DIM` extent.
- `c[i][j] = i + j`.

## Output (`NPDP_DUMP`)

The `n × n` matrix `c`.

## Differences from the original

- **Recurrence:** the original was `c[i][j] += c[i][j-1] + paired(k,j) ? ... : 0`. Because of operator precedence it tested `c[i][j-1] + paired(k,j)`; this version keeps only the conditional term (see counting_im for the corrected recurrence).
- **Initialization:** the original initialization loop was broken (`for i ... for i`, with `j` never set).
- **Overflow:** the arithmetic is unsigned, so overflow is well defined.
- **Base pairing:** the original function call is kept by default (`NPDP_PAIRED(x, y)` expands to `npdp_paired(seq, x, y)`), for backward compatibility. With `-DNPDP_PAIR_TABLE`, it reads a lookup table `pair[x][y]` built before the timer, so the preprocessed scop contains no call. Results are identical in both modes.

## References

- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/counting.cpp`).
