# triang: minimum cost polygon triangulation

## Problem

Triangulate a polygon with vertices `0 … n-1` so that the sum of the triangle costs is minimal. The cost of a triangle is its perimeter. `table[i][j]` is the optimal cost for the sub-polygon `i..j`:

```
table[i][j] = 0                                                          if j - i < 2
table[i][j] = min_{i<k<j} table[i][k] + table[k][j] + cost(i, j, k)      otherwise
cost(i, j, k) = D[i][j] + D[j][k] + D[k][i]                              (perimeter)
```

The table is filled by increasing `gap = j - i`.

## Kernel

```c
static void kernel_triang(int n, double table[NPDP_DIM][NPDP_DIM], double D[NPDP_DIM][NPDP_DIM]) {
#pragma scop
	// initialization
	for (gap = 0; gap < n; gap++)
		for (j = gap; j < n; j++)              // i = j - gap
			if (gap < 2)
				table[j - gap][j] = 0.0;
			else
				table[j - gap][j] = INT_MAX;
	// computation: perfect nest, single statement
	for (gap = 2; gap < n; gap++)
		for (j = gap; j < n; j++)
			for (k = j - gap + 1; k < j; k++)
				table[j - gap][j] = NPDP_MIN(table[j - gap][j], table[j - gap][k] + table[k][j] + (D[j - gap][j] + D[j][k] + D[k][j - gap]));
#pragma endscop
}
```

| Array | Meaning |
|---|---|
| `table` | DP table |
| `D[x][y]` | Euclidean distance between vertices `x` and `y`, precomputed |

## Characteristics

- **Loop nests:** 2. A 2D initialization nest (2 statements under an affine conditional on `gap`), then a perfect 3D nest with a single statement.
- **Cost:** O(n³) time, 2 `double` arrays of `(n+10)²` elements.
- **Parameters:** `n`.
- **Iteration order:** the loops iterate over diagonals (`gap`), and the access functions are skewed (`table[j-gap][...]`).
- **Dependences:** non-uniform, through `table[j-gap][k]` and `table[k][j]`.

## Inputs

- Random integer points in `[0, 100)²` (`NPDP_RSEED`).
- The points do not form a convex polygon; the recurrence is evaluated anyway, as a benchmark.

## Output (`NPDP_DUMP`)

The optimal cost `table[0][n-1]`.

## Differences from the original

- **Point initialization:** upstream never initialized the points (valgrind reported uninitialized reads).
- **Initialization split:** the original initialized `table[j-gap][j]` inside the `gap`/`j` loops, just before the `k` loop. This made the nest imperfect and hard for polyhedral compilers. The initialization is now a separate nest, and the computation is a perfect nest with a single statement. Upstream used this same form to generate the tiled codes in the appendix. Results are bit-identical, because the first `MIN(INT_MAX, x)` returns `x`.
- **No function call:** `cost(i, j, k)` called `dist()` three times inside the scop. A 3D cost table would be too large, so the cost is written inline as three reads of the precomputed distance matrix `D`. The summation order is unchanged, so results are bit-identical.

## Validation

`tests/known_answers.c` checks a square of side 10: two triangles, total cost `40 + 20√2`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest, *Introduction to Algorithms* (1st ed.), section 16.4 "Optimal polygon triangulation".
- M. Palkowski, W. Bielecki, *Accelerating Minimum Cost Polygon Triangulation Code with the TRACO Compiler*, FedCSIS 2018, [PDF](https://annals-csis.org/Volume_17/drp/pdf/8.pdf).
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/triang.cpp`).

## Appendix: reference tiled codes

These tiled versions of the original kernel come from the upstream `info.txt`. They use the original `cost()` function and the macro `N`:

- **R+:** TRACO, tiles based on the transitive closure.
- **TS:** TRACO, tile skewing.
- **PLT:** Pluto.

The kernel now has the same split structure, but they are kept for reference only: they were not tested against this version.

<details>
<summary>Initialization split and tiled codes</summary>

```c
// initialization split from the computation
for (gap = 0; gap < N; gap++){
   for (j = gap; j < N; j++){    // i = j - gap
      if (gap < 2)  // polygon vi...vj has fewer than 3 vertices,
         table[j-gap][j] = 0;
      else
         table[j-gap][j] = INT_MAX;
   }
}

// R+
for( c1 = 2; c1 < N; c1 += 1)
  #pragma omp parallel for
  for( c3 = 0; c3 <= (N - c1 - 1) / 16; c3 += 1)
    for( c5 = 0; c5 <= (c1 - 2) / 16; c5 += 1)
      for( c9 = c1 + 16 * c3; c9 <= min(N - 1, c1 + 16 * c3 + 15); c9 += 1)
        for( c11 = -c1 + 16 * c5 + c9 + 1; c11 <= min(c9 - 1, -c1 + 16 * c5 + c9 + 16); c11 += 1)
          table[c9-c1][c9] = MIN(table[c9-c1][c9], table[c9-c1][c11] + table[c11][c9] + cost(c9-c1,c9,c11));

// TS
for( c0 = 0; c0 <= floord(N - 1, 16); c0 += 1)
  #pragma omp parallel for
  for( c1 = max(0, c0 - (N + 13) / 16 + 1); c1 <= c0; c1 += 1)
    for( c3 = max(16 * c0 + 16 * c1, 16 * c0 - 16 * c1 + 4); c3 <= min(min(2 * N - 16 * c0 + 16 * c1 - 2, N + 16 * c1 + 14), 16 * c0 + 16 * c1 + 45); c3 += 1)
      for( c4 = max(c0 - c1, -2 * c1 + (c3 + 3) / 16 - 2); c4 <= min(min((N - 2) / 16, -c1 + (c3 - 1) / 16), -c1 + (16 * c0 + 16 * c1 + c3 + 13) / 32); c4 += 1)
        for( c6 = max(max(max(max(2, 16 * c1), -N + c3 + 1), -8 * c0 + 8 * c1 + c3 / 2 - 7), -8 * c4 + (c3 + 1) / 2 - 7); c6 <= min(min(16 * c1 + 15, c3 - 16 * c4 - 1), -8 * c0 + 8 * c1 + c3 / 2); c6 += 1)
          for( c10 = max(16 * c4, c3 - 2 * c6 + 1); c10 <= min(16 * c4 + 15, c3 - c6 - 1); c10 += 1)
            table[(c3-c6)-c6][(c3-c6)] = MIN(table[(c3-c6)-c6][(c3-c6)], table[(c3-c6)-c6][c10] + table[c10][(c3-c6)] + cost((c3-c6)-c6,(c3-c6),c10));

// PLT (Pluto)
if (N >= 3) {
  for (t1=0;t1<=floord(N-1,8);t1++) {
    lbp=ceild(t1,2);
    ubp=min(floord(N-1,16),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(2,16*t1-16*t2);t3<=min(N-1,16*t1-16*t2+15);t3++) {
        for (t4=max(16*t2,t3);t4<=min(N-1,16*t2+15);t4++) {
          for (t5=-t3+t4+1;t5<=t4-1;t5++) {
            table[t4-t3][t4] = MIN(table[t4-t3][t4], table[t4-t3][t5] + table[t5][t4] + cost(t4-t3,t4,t5));;
          }
        }
      }
    }
  }
}
```

</details>
