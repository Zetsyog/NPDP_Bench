# NPDP Bench

Non-serial polyadic dynamic programming (NPDP) benchmarks, for evaluating automatic source-to-source optimizing compilers (Pluto, PPCG, Polly, TRACO, ...) and manual loop transformations. The kernels are written in C as static control parts (`#pragma scop`) that polyhedral tools can analyze. OpenMP is supported.

NPDP recurrences, such as `c[i][j] = min_k(c[i][k] + c[k][j])`, have **non-uniform dependences**: their distance depends on the loop iterators. Tiling and parallelizing them is hard, which makes them a demanding test case for polyhedral compilers.

This repository is a cleaned-up rewrite of the [original NPDP Bench](https://github.com/markpal/NPDP_Bench). Many bugs of the original codes are fixed (see [Differences from the original](#differences-from-the-original-npdp-bench)), and correctness checks, reference results and a benchmarking driver are added.

## Benchmarks

| Benchmark | Domain | Loop depth | Time | Memory | Extra parameters |
|---|---|---|---|---|---|
| [counting](counting/README.md) | RNA structure counting | 3 | O(n³) | n² (+n²) | |
| [counting_im](counting_im/README.md) | RNA structure counting (corrected) | 3 | O(n³) | n² (+n²) | |
| [knuth](knuth/README.md) | optimal binary search tree | 3 | O(n³) | 2 n² | |
| [mcc](mcc/README.md) | RNA partition function (McCaskill) | 3 | O(n³) | 2 n² (+n²) | `l`, `ert` |
| [mcm](mcm/README.md) | matrix chain multiplication | 3 | O(n³) | 2 n² | |
| [mea](mea/README.md) | RNA maximum expected accuracy | 4 | O(n⁴) | 5 n² (+n²) | `l`, `ert` |
| [nussinov](nussinov/README.md) | RNA folding (Nussinov) | 3 | O(n³) | n² (+n²) | |
| [nw](nw/README.md) | global alignment (Needleman-Wunsch) | 3 | O(n³) | 4 n² | |
| [sw](sw/README.md) | local alignment (Smith-Waterman) | 3 | O(n³) | 4 n² | |
| [sw3d](sw3d/README.md) | 3-sequence local alignment | 4 | O(n⁴) | 7 n³ | |
| [triang](triang/README.md) | minimum cost polygon triangulation | 3 | O(n³) | 2 n² | |
| [zuker](zuker/README.md) | RNA folding (Zuker-like) | 4 | O(n⁴) | 4 n² | |

`n` is the dataset size (`-DN=...`). Memory counts arrays of `NPDP_DIM² = (n + 10)²` elements (4 or 8 bytes each). "(+n²)" is the pair lookup table, allocated only with `-DNPDP_PAIR_TABLE`.

## Quick start

```bash
nix develop                     # or any Linux with gcc, make, bash, awk (valgrind optional)
make test                       # build everything and run the test suite
cd knuth && make knuth EXTRA_FLAGS="-DN=2000 -DNPDP_TIME" && ./knuth
export PLUTO=/path/to/polycc
make knuth_pluto EXTRA_FLAGS="-DN=2000 -DNPDP_TIME" && OMP_NUM_THREADS=4 ./knuth_pluto
```

## Building

Each benchmark directory has a Makefile (see [common.mk](common.mk)) with these targets:

- `<bench>`: the original sequential program.
- `<bench>_pluto`: the Pluto-optimized program. `<bench>.pluto.c` is generated with `$PLUTO` (`--tile --parallel`, plus any per-benchmark `PLUTO_EXTRA_FLAGS`).
- `clean`: remove the binaries and generated sources.

Variables: `EXTRA_FLAGS` (macros below), `CFLAGS` (default `-march=native -O3 -fopenmp -ffp-contract=off`), `PLUTO`, `PLUTO_FLAGS`, `PLUTO_EXTRA_FLAGS`. The top-level Makefile provides `all` (every original program), `test`, `golden`, `format`, `format-check` and `clean`.

Keep `-ffp-contract=off`. Without it, the compiler fuses `a*b+c` into FMA instructions wherever it can. That changes the floating-point results of `mcc` and `mea` depending on how the loops are written, so original and transformed results would not match.

### Nix

`nix develop` provides gcc (OpenMP), make, bc, clang-format, valgrind and gdb. It keeps `-march=native` and disables the Nix compiler-wrapper hardening flags, so builds match a regular toolchain. Pluto is not packaged in nixpkgs: set `PLUTO` yourself.

## Configuration macros

Set at compile time, for example `make knuth EXTRA_FLAGS="-DN=1000 -DNPDP_TIME"`:

| Macro | Meaning | Default |
|---|---|---|
| `N` | dataset size | 8 |
| `NPDP_TIME` | print the kernel execution time (seconds) on stdout | off |
| `NPDP_DUMP` | print the results on stderr (to compare variants) | off |
| `NPDP_RSEED` | seed of the random inputs | `time(NULL)` |
| `NPDP_FLUSH_CACHE` | flush the caches before starting the timer | off |
| `NPDP_FLUSH_CACHE_SIZE` | bytes used to flush the caches | 64 MiB |
| `NPDP_DIM` | allocated extent of each array dimension | `N + 10` |
| `NPDP_MIN_LOOP` | (mcc, mea) minimum loop length 0–5 | 0 |
| `NPDP_EBP`, `NPDP_RT` | (mcc, mea) base pair energy and temperature | 0, 1 |
| `NPDP_DELTA` | (mea) base pair weighting 1–5 | 1 |
| `NPDP_PAIR_TABLE` | (counting, counting_im, mcc, mea, nussinov) read base pairing from a lookup table instead of calling a function | off |

The default `N` is only meant for quick checks. [tools/run_all.sh](tools/run_all.sh) lists sizes that give sequential runs of a few seconds.

## Conventions for polyhedral compilers

Each kernel follows rules that make it analyzable by any polyhedral tool, including tools that preprocess the source first:

- **Symbolic sizes.** The scop lives in `static void kernel_<name>(int n, ...)`, called from `main` with `N`. `N`, `NPDP_DIM` and the other configuration macros never appear inside a scop, so tools see the parameters `n` and `l`, not constants.
- **Contiguous, statically sized arrays.** Arrays are allocated as single blocks (`NPDP_ALLOC_2D`, `NPDP_ALLOC_3D`) and passed as `T A[NPDP_DIM][NPDP_DIM]` parameters. Every access is affine, and no pointer chasing is involved.
- **No function calls.** Per-element functions of the original codes (`s`, `cost`, `minsq`) are replaced by lookup tables filled before the timer starts: `sim`, `sab`/`sbc`/`sac`, and the distance matrix `D`. The `NPDP_MIN`/`NPDP_MAX` macros expand to conditional expressions.
- **Base pairing (optional table).** For backward compatibility, the RNA kernels (counting, counting_im, mcc, mea, nussinov) keep the original function call by default: `NPDP_PAIRED(x, y)` expands to `npdp_paired(seq, x, y)`. With `-DNPDP_PAIR_TABLE`, it reads a lookup table `pair[x][y]` built before the timer starts, and the preprocessed scop contains no call. Results are identical in both modes. Pluto parses the unpreprocessed source, so it sees `NPDP_PAIRED(x, y)` as a call either way. Tools that preprocess first (PPCG/pet, ...) see the array read.
- **No hidden side effects.** For example, mcm's split table `s` is now written by an explicit statement.

Shared helpers are in [include/mem.h](include/mem.h) and [include/metrics.h](include/metrics.h). Every shared name is prefixed with `NPDP_`/`npdp_`.

## Testing

`make test` (or `nix flake check`, which CI runs on every push) checks:

1. **Warnings:** every benchmark builds with `-Wall -Wextra -Werror`.
2. **Sanitizers:** AddressSanitizer, UndefinedBehaviorSanitizer and LeakSanitizer report nothing.
3. **Golden checksums:** the result dumps at `-O0` and at `-O3 -march=native` match the reference checksums in [tests/golden/](tests/golden/), which are computed at N=64 with seed 1.
4. **Known answers** ([tests/known_answers.c](tests/known_answers.c)) on hand-checkable inputs:
   - mcm: the CLRS example gives 15125.
   - nussinov: small sequences.
   - sw and nw: identical and fully different sequences.
   - triang: a square.
5. **Pair table:** the RNA kernels built with `-DNPDP_PAIR_TABLE` match the same golden checksums and known answers.
6. **Valgrind:** no invalid or uninitialized memory accesses (skipped if valgrind is not installed).

`make golden` regenerates the reference checksums. Only run it when a change of results is intended, and review the diff.

Sources follow [.clang-format](.clang-format): `make format` reformats them, `make format-check` (also run by `nix flake check`) fails on any unformatted file.

## Running experiments

[tools/run_all.sh](tools/run_all.sh) runs every benchmark × variant × thread count and writes one CSV row per measurement:

```bash
PLUTO=/path/to/polycc THREADS="1 2 4 8" tools/run_all.sh -o results.csv
BENCHES="knuth nussinov" SIZES="knuth=3000 nussinov=2500" NB_RUNS=10 tools/run_all.sh -o small.csv
```

For each benchmark, the script:

1. checks every variant's result against the original program (`correct` column, using [tools/compare_dumps.sh](tools/compare_dumps.sh));
2. does one warm-up run;
3. reports the median and standard deviation of the central `NB_MEDIAN_RUNS` of `NB_RUNS` timed runs.

Threads are pinned by default (`OMP_PROC_BIND=close`, `OMP_PLACES=cores`). Each row also records:

- the git commit and whether the tree had uncommitted changes;
- the host, CPU model and frequency governor;
- the gcc and Pluto versions, the flags and the seed.

Everything is built in a temporary copy of the sources. See the header of the script for every option.

Timing methodology:

- The timer uses `clock_gettime(CLOCK_MONOTONIC)` and measures only the kernel call. Input generation and lookup-table construction happen before the timer starts.
- `-DNPDP_FLUSH_CACHE` (`FLUSH_CACHE=1` for the driver) evicts the caches before each run.
- For stable results, use the `performance` CPU governor, an otherwise idle machine and several runs. The driver warns when the governor is not `performance`.

Other tools:

- [tools/get_chksum.sh](tools/get_chksum.sh), run inside a benchmark directory, prints the sha256 of the result dump of one make target. It uses seed 42 unless `-DNPDP_RSEED` is given.

  ```bash
  [ "$(../tools/get_chksum.sh "-DN=500" counting)" = "$(../tools/get_chksum.sh "-DN=500" counting_pluto)" ] && echo OK
  ```

- [tools/compare_dumps.sh](tools/compare_dumps.sh) (`compare_dumps.sh ref cand [rtol] [atol]`) compares two dumps with a numerical tolerance. This is useful for tools that reassociate floating-point operations.
- [tools/finetune_pluto.sh](tools/finetune_pluto.sh) runs an exhaustive search over Pluto tile sizes. It is configured through environment variables; see the script.

## Differences from the original NPDP Bench

The original codes are kept unmodified in [input/](input/). Compared to them:

- **Bug fixes:**
  - All arrays are zero-initialized. Several originals read uninitialized memory: the `nussinov` RNA sequence, the `triang` points and the `sw3d` borders.
  - `sw3d`: all three border faces of `H` are initialized, not only the edges.
  - `sw`, `nw`, `sw3d`: sequences are stored at indices `1..n` to match the 1-indexed kernels.
  - `sw`, `nw`: they score `a[i]` against `b[j]` instead of `b[i]`.
  - `sw`: it subtracts the gap penalty instead of adding it.
  - `nw` is a real Needleman-Wunsch global alignment (no clamp to 0, borders initialized to `-W[i]`). In the original it was a copy of `sw` and did not compile.
  - `mea`: the out-of-bounds read `M[0][-1]` (when `i = k = 0`) is replaced by 0.
  - `counting`, `counting_im`: use `unsigned int`. Their values grow exponentially, and signed overflow is undefined behavior.
  - `counting` uses `c[i][j] += pair[k][j] ? c[i][k-1] + c[k+1][j-1] : 0`. The original had an operator precedence bug; `counting_im` is the corrected version, which adds the `c[i][j-1]` term.
- **Structure:** each benchmark now follows the [polyhedral conventions](#conventions-for-polyhedral-compilers) above: kernel functions, contiguous arrays, lookup tables instead of function calls. Base pairing is the exception: its table is optional (`-DNPDP_PAIR_TABLE`). Results are bit-identical to the fixed codes with function calls.
- **mcm** dumps the split table `s` in addition to the optimal cost.

Each benchmark README has the details.

## Citation and license

If you use these benchmarks in your research, cite this repository ([CITATION.cff](CITATION.cff); GitHub shows a "Cite this repository" button) and the original NPDP Bench by Marek Palkowski. Record the version ([VERSION](VERSION)) and the commit you used. `tools/run_all.sh` writes both into every CSV row.

MIT license ([LICENSE](LICENSE)), the same as the original NPDP Bench, whose copyright notice is kept.
