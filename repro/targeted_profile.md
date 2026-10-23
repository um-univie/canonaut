# Targeted input-boundary and search profiling, 2026-10-23

The two graph instances where the earlier native comparison favored C were
profiled separately: G(8, .5) and G(64, .1), both seed 42. Production code was
not changed during this profiling pass. Baseline is the pre-fix source at
63f8b32; checked and validated refer to the pending fixes described in
[input_validation_review.md](input_validation_review.md).

## Where CPU time goes

Linux process-CPU-time timer sampling, pinned to CPU 2 of the i7-13700,
collected approximately 11,000–13,000 instruction-pointer samples per run.
The timer interval was 1 ms. The profiles cover warmed canonicalization calls;
graph generation, manager construction and initial validation of the reusable
view are outside the sampled region.

| Instance | Old Rust: `refine1` | Checked Rust: `refine1` | Validated Rust: `refine1` | C: `refine1` |
|---|---:|---:|---:|---:|
| G(8, .5) | 44.2% | 44.1% | 44.2% | 57.4% |
| G(64, .1) | 83.9% | 82.9% | 83.1% | 83.7% |

These are shares of each implementation's own CPU samples, not relative
speeds. Inlining can attribute some samples to library helpers. The validated
64-vertex run spends another approximately 10% in adjacency permutation and
2% directly in canonical-graph updating. At eight vertices, setup, dispatch
and canonical output occupy a much larger fraction of the call. Both Rust
cases finish with one search node, one canonical update and no generators;
automorphism callbacks are not exercised by these particular graphs.

Within Rust `refine1`, the singleton-pivot partition scan and split loop
(`src/refine.rs`, approximately lines 616–638) are prominent in both CPU
samples and instruction counts. C's existing native build contains AVX2
partition-end scans (`vpminsd`, `vpcmpeqd`, `vptest`) with scalar handling
around them. Rust's corresponding scan remains scalar, including repeated
32-to-64-bit index conversions. This is an observed code-generation difference
in a hot path and a useful optimization target; no controlled replacement
experiment has yet assigned a fraction of the wall-clock gap to it.

The same C source uses bit-mask-table loads where Rust generates bit-index
manipulation and `bt` instructions. Neither strategy is universally faster.
The existing C makefile specifies GCC with `-O4 -mpopcnt -march=native`.
Rust used release/LTO, one codegen unit, `-C target-cpu=native`, and line debug
information for profiling. The original uninstrumented comparison remains
the timing reference: sampler timings are perturbed and code layout differs.

## What the fixes cost

Callgrind collected only the warmed benchmark loop, using 20,000 calls for
eight vertices and 2,000 for 64 vertices. Counts below remove the fixed 13
loop-setup instructions, but retain equal per-call harness overhead.

| Instance | Old Rust instructions/call | Checked Rust | Validated Rust | C |
|---|---:|---:|---:|---:|
| G(8, .5) | 3,133 | 3,196 | 3,099 | 2,563 |
| G(64, .1) | 50,586 | 50,612 | 50,552 | 40,856 |

The raw checked API adds 63 instructions/call (+2.0%) at eight vertices and
26 (+0.051%) at 64 vertices compared with the old raw API. There is no padding
scan at 64 vertices, but entry checks and wrapper work still exist. Compared
with the new validated entry point, the checked path costs 97 and 60 extra
instructions respectively. The reusable view skips old dimension inference
as well as new padding validation, explaining its lower count than baseline.
This excludes initial view construction; a newly imported graph still needs
validation. Instruction deltas are not cycle or latency deltas.

The old and new Rust refinement and canonical-permutation routines execute
identical instruction counts on both inputs. In the 64-vertex case, Rust
`refine1` including its inlined integer helpers executes 43,046 instructions
per call versus C's 33,957. C's lower work count in the dominant kernel is
much larger than the added entry checks. These profiles do not measure the
graph6 or Python boundary checks; neither is on this canonicalization path.

## Private API boundary assessment

The core modules already have crate-private visibility. That is the right
place for kernels that rely on invariants, but the safe public entry points
must establish those invariants before reaching them:

- `canonize(&[usize])` accepts arbitrary caller-owned words. Making the callee
  private cannot guarantee that padding in this slice is zero.
- `DenseGraph` publicly exposes `adjacency`, `number_of_vertices` and
  `words_per_vertex`. Callers can invalidate its storage after construction.
- The old public callback setter handed safe callbacks mutable permutation
  and orbit slices. Private search state therefore remained externally
  mutable during the search. Immutable callback slices close this hole
  without a validation scan or defensive copy.

`ValidatedGraph` already applies the private-invariant approach: its fields
are inaccessible externally, its constructor validates once, and its immutable
borrow prevents mutation for the duration of use. Repeated canonicalization
then needs no padding scan.

An owned `DenseGraph` could provide the same guarantee without even the import
scan when constructed through checked edge insertion, if its storage and
dimensions became private and every constructor/mutator preserved the
invariants. Raw import would validate once. This would be an additional public
API change, and requires auditing all writers; simply hiding the search
functions or changing checks to debug assertions would not suffice. No such
API redesign or search-kernel optimization was made in this profiling pass.

## Artifacts and reproduction

All profiling artifacts are under
`/scr/kaldra/muhlir/canonaut-safety-fixes/targeted-profile/`:

- `*.callgrind` and `*.instructions.txt`: warmed-loop instruction profiles.
- `*.samples`, `*.summary`, `samples-summary.txt`: CPU samples and symbolized
  summaries. A few mixed C/Rust stdout lines were excluded during parsing.
- `refine-64-annotated.txt`, `rust-refine1.asm`, `c-singleton-scan.asm`: source
  annotation and assembly used for the hot-loop comparison.
- `summarize_samples.py`: symbolization script, run using `uv run python`.

The standalone harness is in `../bench/src/bin/targeted.rs`, with a
benchmark-only signal sampler in `../bench/sampler.c`. Original instruction
profiling and CPU sampling executables are retained in `../profile-native/`
and `../profile-sampler/`, respectively. LTO folds identical old/new helpers,
so some unchanged new-path routines appear under baseline symbol names.

Example instruction-profile command (all paths below relative to the
`canonaut-safety-fixes` artifact directory):

```sh
valgrind --tool=callgrind --collect-atstart=no --toggle-collect=profile_loop \
  --cache-sim=yes --branch-sim=yes \
  --callgrind-out-file=targeted-profile/validated-64.callgrind \
  profile-native/release/targeted validated 64 2000
```

Example CPU sampling command:

```sh
PROFILE_SAMPLES=1 taskset -c 2 profile-sampler/release/targeted \
  validated 64 5000000 > targeted-profile/validated-64.samples
```

Use modes `old`, `checked`, `validated`, `c`, or `validate` (constructor only).
Callgrind's simulated cache/branch results are not hardware measurements and
are not used to explain timing differences here. CPU sampling is statistical;
the percentages identify broad hotspots rather than precise per-line costs.
