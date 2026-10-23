# Input-boundary fixes and performance, 2026-10-22

The five review findings are fixed without adding checks to the refinement,
automorphism, or adjacency-permutation inner loops. Automorphism callbacks now
borrow immutable slices. Graph6 validates input before allocation, and Python
rejects multigraph objects before their edges are converted to adjacency bits.

Raw adjacency and mutable `DenseGraph` inputs must have padding checked on entry:
their contents can change between calls. Only the last word of each row is
examined; a graph whose order is divisible by 64 needs no padding scan. For
single-word rows, the validation uses a vectorizable OR reduction.

For repeated searches, `ValidatedGraph::new(&adjacency)` checks storage once
without allocating or copying. Its immutable borrow prevents invalidation.
`canonize_validated` and `canonize_validated_with_colors` skip subsequent
padding scans and dimension inference. These views preserve the raw API's
directedness convention: enable `with_digraph()` for directed or looped input.
They do not validate symmetry or attach colors.

Mutable raw inputs are checked on every call; the validation-once API avoids
repeating that cost.
Existing benchmark groups still time the checked raw API; the new
`input_validation` Criterion group separately measures checked and validated
calls, with validation outside the latter's timed region.

## Measurements

Baseline: commit `63f8b32`, before these changes. Intel Core i7-13700, CPU 2,
Rust 1.99.1, release/LTO, one codegen unit, warm manager buffers. The baseline,
checked API, validated API and C nauty 2.8.9 run in one process, with 24 rotated
timing rounds per workload. Three separate processes rotate the allocation
order of the three Rust managers as well. Graph construction is outside timing;
validation is included for checked calls and excluded for validated calls.
No callbacks are registered in these timing workloads.

Times below are medians of the three process medians, in microseconds.
The change column is the median of each process's paired validated/baseline
ratios, so it need not equal the ratio of the two displayed median times.
The C library is the existing native build used in the earlier review; no
compiler tuning was applied to C for this comparison.

With both Rust versions built using `-C target-cpu=native`:

| Workload | Baseline µs | Checked µs | Validated µs | C nauty µs | Validated change |
|---|---:|---:|---:|---:|---:|
| G(8, .5) | 0.1115 | 0.1191 | 0.1147 | 0.1128 | +3.2% |
| G(16, .5) | 0.2898 | 0.3055 | 0.3014 | 0.3282 | +4.0% |
| G(36, .5) | 0.6320 | 0.6430 | 0.6382 | 0.8054 | +0.7% |
| G(63, .5) | 1.6886 | 1.6993 | 1.6942 | 2.6531 | −1.1% |
| G(64, .5) | 1.7439 | 1.7259 | 1.7360 | 2.7435 | +0.9% |
| G(65, .5) | 1.8432 | 1.8204 | 1.7957 | 3.7424 | −3.1% |
| G(127, .5) | 7.0177 | 7.0105 | 6.9356 | 13.2075 | −0.9% |
| G(128, .5) | 6.8412 | 6.8113 | 6.7868 | 12.9499 | −0.5% |
| G(129, .5) | 7.5496 | 7.6563 | 7.5144 | 13.4708 | −0.4% |
| G(193, .5) | 14.8063 | 14.6876 | 14.5596 | 31.6782 | −2.8% |
| G(512, .5) | 97.6624 | 97.7598 | 96.6124 | 259.0212 | −1.1% |
| G(64, .1) | 2.3509 | 2.3861 | 2.3875 | 2.2430 | +1.4% |
| G(512, .1) | 53.0053 | 53.3357 | 51.7346 | 85.7542 | −2.7% |
| Path, 63 | 8.8011 | 8.8655 | 8.8510 | 9.9339 | +0.5% |
| Cycle, 65 | 20.9936 | 21.0958 | 20.8755 | 25.6425 | −0.6% |
| Complete, 64 | 323.7414 | 331.5584 | 329.3277 | 486.2264 | +2.2% |
| Colored G(36, .1), 3 colors | 1.3313 | 1.3176 | 1.3091 | 1.3290 | −1.2% |
| Colored G(65, .1), 3 colors | 1.8385 | 1.8343 | 1.8126 | 2.1161 | −1.5% |
| CFI, 200 vertices | 619.3998 | 637.1492 | 610.5468 | 738.6483 | −1.6% |

The geometric mean of the 19 native validated/baseline ratios is 0.9973:
essentially unchanged overall. Individual medians range from 3.1% faster to
4.0% slower. In the portable build, the geometric mean is 0.9801, with individual
medians ranging from 9.0% faster to 2.0% slower. These are local observations,
not proof of a universal zero-cost change. In particular, native G(64, .1)
is slower than C both before and after this change, and G(8, .5) is slightly
slower after the change, so canonaut is not faster than nauty on every graph
or build configuration.

An initial portable G(512, .1) run showed a 24% apparent regression. Rotating
allocation order moved the baseline timing from about 53 to 66 µs while the
new path stayed around 65 µs. The final median paired change was −1.5%.
Callgrind, using separate baseline and new runs with identical input, measured
exactly 1,102,678,000 refinement instructions and 739,280,000 permutation
instructions in each version over 2,000 calls. Total instruction counts were
1,898,476,189 before and 1,898,414,200 after. The allocation-order sensitivity
is consistent with cache-placement effects; instruction counts do not establish
wall-clock equivalence, and individual allocation layouts can still differ.

Graph6 parsing, including the added validation, measured:

| Order | Before µs | After µs |
|---|---:|---:|
| 8 | 0.0414 | 0.0413 |
| 36 | 0.6656 | 0.6594 |
| 65 | 2.1358 | 2.1256 |
| 193 | 18.9928 | 18.6895 |
| 512 | 223.1331 | 218.6857 |

These are 24 alternating rounds on generated G(n, .1), seed 42, portable Rust
builds. Checking each encoded byte once outside the edge loop offsets the new
length/header checks in these samples. Python adapter overhead was not timed.

## Verification and reproduction

- Rust debug and release tests, including boundary rejection, immutable-view
  borrow checking, callback lengths after manager shrinkage, and zero warmed
  allocations for checked and validated entry points.
- Python suite: 29 passing tests, with NetworkX installed.
- C-nauty comparisons and 28,800 structured-graph runs passed after the fixes.
- The independent exhaustive audit again passed 1,966,080 runs over all
  four-vertex digraphs, 15 color partitions, and both Schreier settings,
  checking canonical invariance, orbits, orders, and complete generator closure.
- Valgrind on the new input-boundary integration tests: zero errors.

Audit harnesses and original timing data are retained locally under
`/scr/kaldra/muhlir/canonaut-safety-fixes/`. `bench/src/main.rs` contains the
four-way benchmark; `CANONAUT_BENCH_ORDER=0`, `1`, and `2` rotate allocation
order. `bench-order*.csv` and `native-order*.csv` contain its final data;
`summarize.py` computes the aggregate values. `bench/src/bin/parser.rs` and
`parser.csv` provide the parser measurements. `verify/` contains the expanded
correctness audit; `oracles.log` records its passing results.

With a system nauty installation, run the committed comparison group with:

```sh
cargo bench --features c-nauty-bench --bench bench -- input_validation
```

Set `NAUTY_INCLUDE_DIR` and `NAUTY_LIB_DIR` for nonstandard installations.
