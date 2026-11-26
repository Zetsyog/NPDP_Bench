# sw3d: Smith-Waterman alignment of three sequences

## Problem

This is a local alignment of three sequences `a`, `b`, `c` of length `n`, generalizing [sw](../sw/README.md) to a 3D scoring cube `H[i][j][l]`. Each cell considers:

- `m1`, `m2`, `m3`: a gap of length `k` along a single axis (penalty `2·W[k]`);
- `m4`, `m5`, `m6`: a diagonal step of length `k` in one of the three planes (penalty `W[k]`, plus the pairwise similarity);
- the full diagonal `H[i-1][j-1][l-1]` plus the three pairwise similarities;
- 0 (local alignment).

Similarity is +3 for a match and -3 for a mismatch.

## Kernel (abridged)

```c
static void kernel_sw3d(int n, int H[][NPDP_DIM][NPDP_DIM], int m1..m6[][NPDP_DIM][NPDP_DIM], int W[NPDP_DIM],
                        int sab[][NPDP_DIM], int sbc[][NPDP_DIM], int sac[][NPDP_DIM]) {
#pragma scop
	for (i = 1; i <= n; i++)
		for (j = 1; j <= n; j++)
			for (l = 1; l <= n; l++) {
				m1[i][j][l] = INT_MIN;
				for (k = 1; k <= i; k++) m1[i][j][l] = NPDP_MAX(m1[i][j][l], H[i - k][j][l] - 2 * W[k]);
				/* m2 (j axis), m3 (l axis): same pattern */
				m4[i][j][l] = INT_MIN;
				for (k = 1; k <= i && k <= j; k++) m4[i][j][l] = NPDP_MAX(m4[i][j][l], H[i - k][j - k][l] - W[k] + sab[i][j]);
				/* m5 (j, l plane), m6 (i, l plane): same pattern */
				H[i][j][l] = NPDP_MAX(0, NPDP_MAX(H[i-1][j-1][l-1] + sab[i][j] + sac[i][l] + sbc[j][l], max(m1..m6)));
			}
#pragma endscop
}
```

See [sw3d.c](sw3d.c) for the exact code. Bounds such as `k <= min(i, j)` are written `k <= i && k <= j`, so that they contain no function or macro call.

## Characteristics

- **Loop nest:** 4 loops (six sibling `k` loops), 13 statements, at depths 3 and 4.
- **Cost:** O(n⁴) time. Memory is 7 arrays of `(n+10)³` ints: about 340 MB at n = 220 and 3.4 GB at n = 500.
- **Parameters:** `n`.
- **Dependences:** non-uniform along every axis and every plane diagonal, plus the uniform `(1, 1, 1)` dependence. Parallelism is a 3D wavefront.

## Inputs

- Three random sequences over {A, C, G, U} (`NPDP_RSEED`), stored at indices `1..n`.
- `W[k] = 2k`.
- `H` starts at 0, including its three border faces.

## Output (`NPDP_DUMP`)

The `n × n × n` cube `H`, one `n × n` slice per block.

## Differences from the original

- **Border initialization:** only the three edges `H[i][0][0]`, `H[0][i][0]` and `H[0][0][i]` were initialized, but the kernel reads the three border faces (`H[i][j][0]`, ...). Valgrind reported uninitialized reads. Everything is now zero-initialized.
- **Sequence indexing:** sequences are 1-based, matching the kernel.
- **No function calls:** the three `s(x, y)` calls are replaced by the pairwise tables `sab`, `sbc` and `sac`. `min(i, j)` in loop bounds is written as a conjunction.
- **Cleanup:** the unused header `mem3d.h` is removed.

## References

- T. F. Smith, M. S. Waterman, Journal of Molecular Biology 147(1), 195–197, 1981 (pairwise algorithm).
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/sw3d.cpp`).
