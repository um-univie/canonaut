# Entry-point and setup optimization, 2026-10-23

This pass compares against the working tree **after the input-safety fixes**,
before these optimizations. It does not compare against the pre-fix commit.
The snapshot is retained under
`/scr/kaldra/muhlir/canonaut-perf-20261023/baseline/`.

## Retained changes

- `ValidatedGraph::new` recognizes zero-to-64-word inputs directly: their order
  equals their length. General dimension inference is kept out of this fast
  path. Padding validation is still performed when needed.
- The private manager runner takes a `ValidatedGraph`, so it does not repeat
  the dimension check already established at the public boundary.
- Output preparation checks the common logical length of the private labeling,
  partition and orbit buffers once. Public canonical output remains independently
  checked because callers can resize it or take ownership of it.
- Algorithm-buffer preparation uses the jointly maintained vertex-buffer size
  and exact active-set word count to skip repeated capacity checks. Buffers
  still grow as needed, and changing graph width still updates the active set.
  Debug assertions check the private invariants on the fast path.

No unchecked public entry point was added. The refinement implementation is
unchanged. No additional public API change is needed for these optimizations.

## Measurements

Intel Core i7-13700, CPU 2, release/LTO, one codegen unit. Both baseline and
candidate are in the same benchmark binary, with preallocated managers and
identical graphs. Checked calls include validation; validated calls exclude
view construction. Four processes rotate allocation order of all four Rust
managers. Timed execution order rotates as well. Graph construction is excluded.

The table reports geometric means of per-workload ratios. Each workload ratio
is the median of four process-level paired timing ratios. A negative value
means lower elapsed time. Native uses `-C target-cpu=native`; portable uses
the compiler's default CPU target.

| Workloads | Build | Checked change | Validated change |
|---|---|---:|---:|
| Original 19 instances, seed 42 where random | Native | −0.15% | −0.92% |
| 80 additional random instances, seeds 0–7 | Native | −1.78% | −2.50% |
| Original 19 instances | Portable | −2.01% | −3.84% |
| 80 additional random instances | Portable | −4.76% | −4.88% |

The additional instances cover orders 8, 16, 36, 63 and 64 at edge probabilities
0.1 and 0.5. The original set also includes multiword graphs, colored graphs,
paths, cycles, complete graphs and CFI. Original cases use 24 timing rounds;
the additional sweep uses 12. Every case checks canonical output against both
the baseline and C nauty.

For example, on the native seed-42 cases, checked G(8, .5) changed from about
110.0 to 106.2 ns, and checked G(64, .1) from 2.087 to 2.018 µs. These are
medians of process medians; their ratios need not equal the paired ratios above.

This is a modest aggregate improvement, not a promise that every input gets
faster. On the original native set, checked per-case changes range from about
−3.7% to +3.0%; validated changes range from −2.6% to +2.6%. Timings include
compiler/code-placement and allocator-placement effects. In particular, the
large portable G(512, .1) gain is layout-sensitive, as in the earlier profiling;
it should not be attributed entirely to fewer setup checks. The diverse small
graph sweep is a separate check on the entry-point optimizations.

Callgrind provides independent evidence that the retained changes reduce work.
With 2,000 warmed calls per case, removing the fixed 13 loop-setup instructions
but retaining equal per-call harness overhead gives:

| Native instance | Baseline checked instructions/call | Optimized checked | Difference |
|---|---:|---:|---:|
| G(8, .5) | 3,170 | 3,121 | −49 |
| G(64, .1) | 50,586 | 50,534 | −52 |

Instruction counts are not latency estimates. The immutable-view path executes
3,038 and 50,491 instructions respectively, excluding initial validation.

## Refinement alternatives examined

The hot scalar partition scan was tested with native-width indices, bounded
SSE2 probes, a maintained bitset of cell boundaries, and four-way scalar
unrolling. Some improved selected seed-42 cases, but broader results did not
justify keeping them. The bounded SIMD variant slowed the 63-vertex path by
about 22% in its exploratory run. The boundary bitset also regressed several
families. Four-way unrolling regressed the additional 80-instance native sweep
by about 3% for checked calls and 2.4% for validated calls.

All of those kernel changes were reverted. This avoids selecting an optimization
solely because it improves the single graph that initially motivated profiling.
The existing scalar loop remains a possible future optimization target, but
none of these replacements demonstrated a sufficiently broad benefit.

## Verification

- Rust debug and release suites: 134 passing tests including doctests; 28
  corpus-dependent tests ignored. Warmed allocation-count test still passes.
- Added a reuse regression covering public pre-sizing, growing/shrinking graph
  orders, arbitrary canonical-output resizing, and taking the output buffer.
- Independent exhaustive oracle: 1,966,080 manager runs on all four-vertex
  digraphs, all 15 color partitions, and both Schreier settings.
- Expanded C comparisons and structured/relabeling/group-generator checks pass.
- Valgrind on the new reuse regression and input-boundary tests: zero errors.
- Clippy on all default-feature targets with warnings denied: passes.

## Reproduction artifacts

The directory `/scr/kaldra/muhlir/canonaut-perf-20261023/` contains:

- `bench/src/main.rs`: five-way benchmark (baseline checked/validated,
  candidate checked/validated, C).
- `bench/src/bin/sweep.rs`: additional graph-seed sweep.
- `capacity-native-*.csv`, `capacity-portable-*.csv`: retained-version results
  for the original cases.
- `capacity-sweep-native-*.csv`, `capacity-sweep-portable-*.csv`: retained-version
  results for the 80 additional instances.
- `summarize.py`, `kept-summary.txt`: aggregation script and output.
- `kept-*.callgrind`: warmed instruction profiles.
- `tests-kept-*.log`, `oracles-kept.log`, `valgrind-*.log`, `clippy-kept.log`:
  verification logs.

Other prefixes retain exploratory/rejected variants; they are not the final
implementation. In particular, the earlier `final-native-*` files belong to
the intermediate scalar-index experiment. Use the `capacity-*` files above.

Build the standalone benchmark with its manifest, a task-local
`CARGO_TARGET_DIR`, and optionally `RUSTFLAGS='-C target-cpu=native'`. Run with
`CANONAUT_BENCH_ORDER=0`, `1`, `2`, and `3`, pinned to the same CPU. The harness
uses the existing local C-nauty archive through `bench/build.rs`. Run aggregation
with `uv run python summarize.py 'capacity-native-*.csv'` from this artifact
directory (the script resolves its stored artifact root).
