# counting_im: RNA structure counting recurrence (corrected)

## Problem

This is the version of [counting](../counting/README.md) corrected by Piotr Blaszynski. After summing the contributions of every position `k` that pairs with `j`, it adds the count of the subsequence without `j`:

```
c[i][j] += Σ_{k = i}^{j - 1} [pair(k, j)] · (c[i][k-1] + c[k+1][j-1])
c[i][j] += c[i][j-1]
```

Like counting, it uses sums and synthetic initial values (`c[i][j] = i + j`). It reproduces the dependence pattern of the RNA structure counting recurrence, not actual counts.

## Kernel

```c
static void kernel_counting_im(int n, NPDP_INT_TYPE c[NPDP_DIM][NPDP_DIM], NPDP_PAIR_PARAM) {
#pragma scop
	for (int i = n - 2; i >= 1; i--)
		for (int j = i + 2; j <= n; j++) {
			for (int k = i; k <= j - 1; k++)
				c[i][j] += NPDP_PAIRED(k, j) ? c[i][k - 1] + c[k + 1][j - 1] : 0;
			c[i][j] = c[i][j] + c[i][j - 1];
		}
#pragma endscop
}
```

| Array | Meaning |
|---|---|
| `c` | DP table, `unsigned int` (wrapping overflow) |
| `NPDP_PAIRED(x, y)` | 1 if nucleotides `x` and `y` can pair (A-U, G-C, G-U): an `npdp_paired(seq, x, y)` call by default, as in the original, or a read of the precomputed table `pair[x][y]` with `-DNPDP_PAIR_TABLE` |

## Characteristics

- **Loop nest:** 3 loops, 2 statements. The second statement sits at depth 2, after the `k` loop, so the loop nest is imperfect.
- **Cost:** O(n³) time, 1 array of `(n+10)²` elements, plus the `int` pair table with `-DNPDP_PAIR_TABLE`.
- **Parameters:** `n`.
- **Dependences:** non-uniform, as in counting, plus `c[i][j-1]`, a uniform dependence carried by `j`.

## Inputs

- A random sequence over {A, C, G, U} (`NPDP_RSEED`).
- `c[i][j] = i + j`.

## Output (`NPDP_DUMP`)

The `n × n` matrix `c`.

## Differences from the original

- **Overflow:** the arithmetic is unsigned. Signed overflow was undefined behavior from N≈40.
- **Base pairing:** the original function call is kept by default (`NPDP_PAIRED(x, y)` expands to `npdp_paired(seq, x, y)`), for backward compatibility. With `-DNPDP_PAIR_TABLE`, it reads a lookup table `pair[x][y]` built before the timer, so the preprocessed scop contains no call. Results are identical in both modes.
- **Allocation:** the original allocated the sequence with `sizeof(char *)` per element.

## References

- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`counting/counting_im.cpp`).
