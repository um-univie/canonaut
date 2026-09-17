# Correctness and hotspot work, 2026-09-16

The baseline is the working tree before this change, **including** the existing
three-line `gca_first` backtracking correction. Measurements are local to an
Intel Core i7-13700, pinned to CPU 2, Rust release with LTO and one codegen unit.
The baseline and changed library ran in the same process, with 30 alternating
A/B rounds and warm buffers. Graph generation was outside the timed region.
Speedups below are medians of paired ratios; timings are means in microseconds.
They are observations on these instances, not guarantees across graph families
or CPUs. CFI's parameter 20 produces a 200-vertex graph.

| Workload | Before µs | After µs | Paired speedup |
|---|---:|---:|---:|
| Random G(64, .1) | 2.704 | 2.542 | 1.065× |
| Random G(64, .5) | 3.612 | 1.660 | 2.148× |
| Random G(512, .1) | 92.230 | 55.241 | 1.641× |
| Random G(512, .5) | 313.874 | 113.983 | 2.760× |
| Cycle, 128 vertices | 89.991 | 84.744 | 1.059× |
| Complete, 64 vertices | 483.742 | 323.837 | 1.501× |
| Colored G(36, .1), 3 colors | 1.477 | 1.337 | 1.111× |
| Colored G(64, .1), 3 colors | 2.526 | 2.255 | 1.124× |
| CFI, base order 20 | 720.171 | 726.149 | 0.999× |

The two real colored molecular fixtures in `tests/manager_reuse.rs` were also
benchmarked separately: 36 vertices went from 1.559 to 1.548 µs (1.007×), and
91 vertices from 22.012 to 21.867 µs (1.010×). Treat these as essentially
unchanged. The large dense-graph gains do not imply comparable gains on small,
sparse, heavily colored molecules.

With both libraries compiled using `-C target-cpu=native`:

| Workload | Before µs | After µs | Paired speedup |
|---|---:|---:|---:|
| Random G(64, .5) | 3.366 | 1.614 | 2.090× |
| Random G(512, .5) | 250.681 | 108.873 | 2.313× |
| Complete, 64 vertices | 406.737 | 334.399 | 1.207× |
| Colored G(64, .1), 3 colors | 2.426 | 2.160 | 1.128× |
| CFI, base order 20 | 633.261 | 650.851 | 0.972× |

The native CFI sample is about 3% slower; cycle and CFI results generally stay
near the baseline. The improvement is largest when permuting many adjacency
bits or checking many automorphism edges. No CPU-specific flags were added to
the distributed crate's default build.

## Hotspots and changes

Callgrind on 100 baseline calls to G(512, .1) attributed 44% of all instructions
to lines in multiword refinement and 32% to adjacency permutation, excluding
additional instructions attributed to inlined standard-library code. Complete
K64 put substantial work in automorphism checking as well as refinement.

Adjacency permutation and automorphism checking now iterate set bits using
`trailing_zeros` and `word &= word - 1`, converting each bit position back to
MSB-0 indexing. The previous loop computed a shift from `leading_zeros` before
it could clear the bit and advance. Separating bit clearing from index lookup
shortens the loop's dependency chain. Bit order is irrelevant to permutation
output and to the automorphism predicate. This does **not** change refinement
ordering, search target selection, or canonical-label comparison ordering.
Instruction count alone understates the improvement: Callgrind's total for the
random workload went from 104.4 million to 101.8 million instructions, while
native wall-clock time improved substantially.

Colored setup now sorts the manager's labeling buffer in place, with an explicit
vertex-index tie-breaker. It no longer allocates an order Vec or stable-sort
scratch space on each call. Known dimensions bypass inference, and colored
calls no longer enter setup twice. Raw input dimensions are inferred using an
integer square root instead of trial division. Canonical output is initially
allocated for n*m words. The allocation regression verifies **zero heap
requests** through repeated, warmed colored/uncolored size transitions.

Python result conversion fetched `labeling` and `orbits` through PyO3 getters
once per vertex, materializing entire lists repeatedly. Fetching each once
removes that quadratic copying. An isolated conversion benchmark, using the
same extension result on both sides, measured median times:

| n | Before µs | After µs | Speedup |
|---|---:|---:|---:|
| 64 | 40.75 | 9.94 | 4.10× |
| 512 | 3917.17 | 1080.52 | 3.63× |

These are result-conversion timings, **not** end-to-end Python call speedups.
Conversion of canonical edges still has cost proportional to the edge count.

## Correctness

- The active bitset's logical length now follows m when the graph shrinks,
  retaining allocation capacity. Stale active bits previously caused unchecked
  out-of-bounds reads and lost automorphism generators.
- Empty graphs, explicit active bitsets and missing canonical-output validation
  are covered by regressions. DenseGraph checks edge endpoints and uses full
  automorphism checking for loops.
- Python enables directed search, distinguishes graph kinds, validates endpoints,
  and gives fixed-color keys a deterministic ordering and complete encoding.
- Internal kernels/options are private. The subsequent Schreier repair exposes
  optional stabilizer pruning through `CanonautManager::with_schreier()`.
- Python canonical-byte keys have a new versioned format; stored keys must be
  regenerated. Supported color labels are stable built-in values, including
  tuples/frozensets. Custom objects must first be mapped to stable color labels.

Validation: all 136 library tests including the full reference corpora, nine
manager/graph regression tests, the allocation test, and 28 Python tests pass.
The optional C test additionally checks 68,600 graphs against nauty 2.8.9,
comparing canonical adjacency, orbit arrays and group orders with manager reuse.
It covers all labeled 4-vertex simple undirected graphs, all 4-vertex digraphs
including loops, and 3,000 mixed-size random cases with and without colors.
The C differential test passes in both portable and CPU-native release builds.
Workspace Clippy and documentation checks are also run.

Miri was attempted but could not start: the installed
`nightly-2025-06-05` toolchain has component metadata but lacks its Cargo binary.
The reuse regression passes normally in debug mode, where its original
unchecked out-of-bounds access caused an abort. This is not a substitute for a
completed Miri run.

## Reproducing checks and profiles

Use an artifact root appropriate to your machine; this session used
`/scr/kaldra/muhlir/canonaut-fixes`. Set `CARGO_TARGET_DIR` there for all commands.
For a nonstandard C installation, set `NAUTY_INCLUDE_DIR`, `NAUTY_LIB_DIR` and
the platform's dynamic-library search path.

```sh
cargo test --release --lib --tests -- --include-ignored
cargo test --release --features c-nauty-bench --test c_reference
cargo clippy --workspace --all-targets -- -D warnings
cargo test --doc
cargo build --profile flamegraph --bin profile
```

The profiling binary now accepts `profile <family> <n> [iterations]`.
Supported families include `random` (p=.5), `random10`, `colored` (p=.1, three
colors), `complete`, `cycle`, `path`, `star`, `cfi`, `hadamard` and `paley`.
Use its path under `$CARGO_TARGET_DIR/flamegraph` with a sampling profiler or
Callgrind; direct profile outputs to your artifact root. Example workloads are
`random10 512 100`, `complete 64 100`, and `colored 64 10000`.

The baseline source snapshot, paired benchmark harness, logs and Callgrind
outputs for this session are retained in
`/scr/kaldra/muhlir/canonaut-fixes/`.

## Schreier stabilizer repair

Completed the optional path using nauty 2.8.9 `filterschreier`, `expandschreier`
and `pruneset` as the reference. The chain now has an unfixed terminal level,
identity root markers, inverse-power transversal paths, and a base aligned with
all currently fixed search vertices. Filtering preserves essential generators;
derived permutations enter the random walk and are reclaimed when no longer
referenced after a base change. Reset retains storage and resizes logical buffers
when the manager changes graph order. Search passes only the active permutation
slice. The former rejection guard is removed; the public opt-in is
`CanonautManager::with_schreier()`.

As in nauty, expansion stops after ten consecutive unsuccessful random words.
This is not a certificate that every stabilizer element has been discovered.
Every discovered element belongs to the pointwise stabilizer, so incomplete
expansion can only miss pruning opportunities, never justify invalid pruning.

Validation: 138 library tests (including ignored reference corpora), nine manager
regressions and the allocation regression pass. The C comparison now runs its
68,600 cases in both modes (137,200 comparisons), checking canonical adjacency,
orbits and group order, plus 45 symmetric/colored cases covering complete,
cycle, star, CFI, Hadamard and Paley graphs. Small-group unit tests explicitly
enumerate group closures and check every removed vertex against the pointwise
stabilizer, including incompatible base changes, generator additions and reset.
A separate S6 test verifies actual pruning with two fixed vertices. Permutation
power tests cover all three implementation paths across a word boundary.
C reference calls are serialized because the reference build has shared scratch
storage. Workspace Clippy passes with warnings denied.

Seven alternating pairs, pinned CPU 2, portable release build; median ns/graph:

| Family | Order | Disabled | Enabled |
|---|---:|---:|---:|
| Random, p=0.5 | 64 | 1,650 | 1,660 |
| Complete | 64 | 324,783 | 331,266 |
| Complete | 193 | 8,480,808 | 8,542,552 |
| Cycle | 64 | 19,002 | 19,341 |
| CFI | 200 | 709,676 | 715,747 |
| Hadamard | 64 | 52,843 | 54,740 |
| Paley | 61 | 8,652 | 9,140 |

Node counts were unchanged in these samples; overhead was approximately 1–6%.
Keep the option disabled by default. These measurements do not establish a
speedup on other symmetric families. Reproduce with
`profile <family> <order> <iterations> [--schreier]`; CFI and Hadamard take base
and matrix order respectively, as documented in that binary.
