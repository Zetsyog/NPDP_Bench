# zuker: Zuker-like RNA folding

## Problem

This kernel has the loop structure of Zuker's minimum free energy folding algorithm:

- `V[i][j]` is the minimum energy of `i..j` when `i` and `j` are paired. It considers hairpin loops (`EHF`), interior loops closed by an inner pair `(k, m)` whose size `k - i + j - m` is limited to the window (2, 30), and multi-branch decompositions `W[i+1][k] + W[k+1][j-1]`.
- `W[i][j]` is the minimum energy of `i..j`, with bifurcations `W[i][k]`, `W[k+1][j]`.

The energy tables are synthetic (`EFL[i][j] = EHF[i][j] = i + 1`, initial `W[i][j] = i·j`, `V[i][j] = i + 1`). The kernel reproduces the dependence structure of Zuker's algorithm, not real free energies.

## Kernel

```c
static void kernel_zuker(int n, int W[NPDP_DIM][NPDP_DIM], int V[NPDP_DIM][NPDP_DIM], int EFL[NPDP_DIM][NPDP_DIM], int EHF[NPDP_DIM][NPDP_DIM]) {
#pragma scop
	for (i = n - 1; i >= 0; i--)
		for (j = i + 1; j < n; j++) {
			for (k = i + 1; k < j; k++) {
				for (m = k + 1; m < j; m++)
					if (k - i + j - m > 2 && k - i + j - m < 30)
						V[i][j] = NPDP_MIN(V[k][m] + EFL[i][j], V[i][j]);
				W[i][j] += NPDP_MIN(NPDP_MIN(W[i][k], W[k + 1][j]), W[i][j]);
				if (k < j - 1)
					V[i][j] = NPDP_MIN(W[i + 1][k] + W[k + 1][j - 1], V[i][j]);
			}
			V[i][j] = NPDP_MIN(NPDP_MIN(V[i + 1][j - 1], EHF[i][j]), V[i][j]);
			W[i][j] = NPDP_MIN(NPDP_MIN(NPDP_MIN(W[i + 1][j], W[i][j - 1]), V[i][j]), W[i][j]);
		}
#pragma endscop
}
```

## Characteristics

- **Loop nest:** 4 loops, 5 statements, at depths 4, 3, 3, 2 and 2.
- **Cost:** O(n⁴) loop iterations. The interior-loop statement only executes inside the window of the affine guard `2 < k - i + j - m < 30`, but the loops still traverse the full 4D space. Memory is 4 arrays of `(n+10)²` ints.
- **Parameters:** `n`.
- **Dependences:** `V` and `W` depend on each other. Accesses include `V[k][m]` (a 2D inner interval), `W[i][k]`, `W[k+1][j]` and `W[i+1][k]`, so the dependences are non-uniform. Affine conditionals restrict the iteration domain.

## Inputs

Synthetic tables (see above). There is no random input, so `NPDP_RSEED` has no effect.

## Output (`NPDP_DUMP`)

The `n × n` matrices `V` and `W`.

## Differences from the original

- Arrays are contiguous and passed as parameters.
- Three unused arrays (`VMI`, `VBI`, `WZ`) are removed.
- Upstream did not compile with recent gcc (`floor`/`ceil` used without `<math.h>`).

## References

- M. Zuker, P. Stiegler, *Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information*, Nucleic Acids Research 9(1), 133–148, 1981.
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/zuker.cpp`).
