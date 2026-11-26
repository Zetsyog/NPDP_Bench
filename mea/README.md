# mea: maximum expected accuracy RNA folding

## Problem

MEA folding predicts the RNA secondary structure that maximizes the expected number of correctly predicted base pairs and unpaired nucleotides. Each pair is weighted by its probability. The pipeline has five stages, each a separate loop nest in a single scop:

1. **Partition function `Q`, `Qbp`:** as in [mcc](../mcc/README.md).
2. **Base-pair probabilities `Pbp[i][j]`:** a 4-deep nest. Every enclosing pair `(p, q)` with `p < i < j < q` contributes.
3. **Probabilities `Pu[i][j]`:** a 4-deep nest, same structure as stage 2.
4. **Unpaired probabilities `Puu[i] = 1 - Σ_j Pbp[i][j]`:** 2 loops.
5. **MEA dynamic program `M[i][j]`:** pairs weighted by `NPDP_DELTA · Pbp`, unpaired nucleotides by `Puu`.

The probability formulas are simplified versions of the real ones. The kernel is a benchmark of the loop structure and its dependences, not a validated folding tool.

## Kernel (abridged)

```c
static void kernel_mea(int n, int l, double ert, double Q[][], double Qbp[][], double Pbp[][], double Pu[][], double M[][], double Puu[], NPDP_PAIR_PARAM) {
#pragma scop
	/* 1 */ if (n >= 1 && l >= 0 && l <= 5) for i (n-1..0) for j (i+1..n-1) { Q[i][j] = Q[i][j-1]; for k (0..j-i-l-1) { Qbp[k+i][j] = ...; Q[i][j] += ...; } }
	/* 2 */ for i (0..n-1) for j (i+1..n-1) { Pbp[i][j] = Q[0][i] * Q[j][n-1] * Qbp[i][j] / Q[0][n-1];
	                                       for p (0..i-1) for q (j+1..n-1) Pbp[i][j] += Pbp[p][q] * ert * Q[p+1][i] * Qbp[i][j] * Q[j+1][q-1] / (Qbp[p][q] == 0 ? 1 : Qbp[p][q]); }
	/* 3 */ for i (n-1..0) for j (i+1..n-1) { Pu[i][j] = ...; for p (0..i-1) for q (j+1..n-1) Pu[i][j] += ...; }
	/* 4 */ for i (0..n) { Puu[i] = 1; for j (i+1..n-1) Puu[i] -= Pbp[i][j+1]; for k (0..i-1) Puu[i] -= Pbp[k][i+1]; }
	/* 5 */ for i (n-1..0) for j (i+1..n-1) { for k (0..j-i-l-1)
	                                            M[i][j] = MAX(M[i][j], (k > 0 ? M[i][k+i-1] : 0) + M[k+i+1][j-1] + NPDP_DELTA * Pbp[k+i][j]) * NPDP_PAIRED(k+i, j-1);
	                                         M[i][j] = MAX(M[i][j], M[i][j-1] + Puu[j-1]); }
#pragma endscop
}
```

See [mea.c](mea.c) for the exact code.

## Characteristics

- **Loop nests:** 5 nests in one scop. Their depths are 3, 4, 4, 2 and 3, with 12 statements in total. It is the largest scop of the suite.
- **Cost:** O(n⁴) time (stages 2 and 3). Memory is 5 `double` arrays of `(n+10)²` elements, plus the `int` pair table with `-DNPDP_PAIR_TABLE`.
- **Parameters:** `n`, `l` (`NPDP_MIN_LOOP`) and `ert`.
- **Dependences:**
  - stage 2 reads `Pbp[p][q]` of enclosing pairs, so it is computed from the outside in (`i` increasing, `q > j`);
  - the stages are chained through `Q`, `Qbp`, `Pbp` and `Puu`, so fusion and inter-nest scheduling matter.
- **Data-dependent conditionals:** the expressions `Qbp[p][q] == 0 ? 1 : ...` and `k > 0 ? ... : 0` are conditional, but the control flow stays static.
- **Floating point:** compile with `-ffp-contract=off` to compare results exactly.

## Inputs and macros

- A random sequence over {A, C, G, U} (`NPDP_RSEED`).
- `NPDP_MIN_LOOP` (0), `NPDP_EBP` (0), `NPDP_RT` (1), `NPDP_DELTA` (1).

## Output (`NPDP_DUMP`)

- `ert` and the sequence.
- The upper triangles of `Q`, `Qbp`, `Pbp`, `Pu` and `M`.
- `Puu`.

## Differences from the original

- **Out-of-bounds read fixed:** when `i = k = 0`, `M[i][k+i-1]` read `M[0][-1]`. The term is now 0, the value of the (zero) lower triangle used for every other `i`.
- **Single scop:** the five nests are in one scop. Upstream they were five separate scops; compilers can now optimize across stages.
- **Base pairing:** the original function call is kept by default (`NPDP_PAIRED(x, y)` expands to `npdp_paired(seq, x, y)`), for backward compatibility. With `-DNPDP_PAIR_TABLE`, it reads a lookup table `pair[x][y]` built before the timer, so the preprocessed scop contains no call. Results are identical in both modes.
- **Parameters:** `ERT` and `l` are kernel parameters instead of a global and a macro.
- **Macro names:** `Ebp`, `RT`, `l` and `delta` are renamed `NPDP_EBP`, `NPDP_RT`, `NPDP_MIN_LOOP` and `NPDP_DELTA`.
- **Input:** upstream overrode the random sequence with the fixed string `"GGUCCAC"` and did not compile as C (K&R `paired(i, j)`, missing `../mem.h`).

## References

- Z. J. Lu, J. W. Gloor, D. H. Mathews, *Improved RNA secondary structure prediction by maximizing expected pair accuracy*, RNA 15(10), 1805–1813, 2009.
- J. S. McCaskill, Biopolymers 29, 1105–1119, 1990 (partition function).
- Original code: [markpal/NPDP_Bench](https://github.com/markpal/NPDP_Bench) (`input/mea.cpp`).
