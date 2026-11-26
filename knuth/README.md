# knuth: optimal binary search tree recurrence

## Problem

This is the dynamic program for optimal binary search trees (Knuth 1971). It is the canonical "interval" recurrence: the optimal cost of the interval `(i, j)` is the best split `k` of the interval into `(i, k)` and `(k, j)`:

```
c[i][j] = min( c[i][j], min_{i<k<j} ( w[i][j] + c[i][k] + c[k][j] ) )
```

The kernel uses the cubic form, not Knuth's O(n²) speed-up that restricts the range of `k`. The weights and initial costs are synthetic (`w[i][j] = i - j`, `c[i][j] = i + j`).

## Kernel

```c
static void kernel_knuth(int n, int c[NPDP_DIM][NPDP_DIM], int w[NPDP_DIM][NPDP_DIM]) {
#pragma scop
	for (i = n - 1; i >= 1; i--)
		for (j = i + 1; j <= n; j++)
			for (k = i + 1; k < j; k++)
				c[i][j] = NPDP_MIN(c[i][j], w[i][j] + c[i][k] + c[k][j]);
#pragma endscop
}
```

| Array | Meaning |
|---|---|
| `c` | cost table |
| `w` | weight table (read only) |

## Characteristics

- **Loop nest:** 3 loops, 1 statement.
- **Cost:** O(n³) time, 2 arrays of `(n+10)²` ints.
- **Parameters:** `n`.
- **Dependences:** `c[i][j]` depends on `c[i][k]` (same row) and `c[k][j]` (same column, a lower row), with distances that depend on `k`. This is the typical NPDP pattern, and it is also the simplest kernel of the suite.

## Inputs

- `c[i][j] = i + j`, `w[i][j] = i - j`.
- There is no random input, so `NPDP_RSEED` has no effect.

## Output (`NPDP_DUMP`)

The `n × n` matrix `c`.

## Differences from the original

- Arrays are contiguous and passed as parameters.
- Otherwise identical to `knuth/knut_fix.cpp` upstream. Upstream's `knut.cpp` read an uninitialized array.

## References

- D. E. Knuth, *Optimum binary search trees*, Acta Informatica 1, 14–25, 1971.
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/knut.cpp`).
