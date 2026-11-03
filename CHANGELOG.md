# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 2026-11-03

### Input validation

- Support 32-bit targets such as wasm32: word size, bit and word indices come
  from `usize::BITS`, and the random generator keeps 64-bit state, so labelings
  match across targets (#1, by S. Joshua Swamidass).
- Reject nonzero adjacency padding before unchecked graph kernels. Add the
  allocation-free `ValidatedGraph` view and `canonize_validated` /
  `canonize_validated_with_colors` for validation once across repeated searches.
- **Breaking:** automorphism callbacks now take `&[u32]` for permutations and
  orbits, preventing safe callbacks from corrupting internal search state.
  Both slices have the current graph order as their length, including after reuse.
- Reject truncated/extra graph6 payloads, invalid order bytes, nonzero format
  padding, and unsupported orders before adjacency allocation. Validate each
  payload byte once outside the edge decoding loop. `Graph6Error` has new
  variants; exhaustive matches must handle them.
- Reject Python multigraph objects explicitly instead of silently dropping
  parallel-edge multiplicities.
- Add separate checked/validated input benchmarks; see
  `repro/input_validation_review.md` for before/after and C-nauty measurements.
- Reduce repeated manager setup checks using private buffer invariants and
  add a single-word input fast path. Keep the existing refinement scan after
  broader benchmarks rejected alternatives; see `repro/setup_optimization.md`.

### Correctness and performance

- Fix stale active-cell bits and unchecked out-of-bounds reads when reusing a
  manager on smaller graphs; verify group orders and orbit partitions as well
  as canonical adjacency.
- Support empty graphs through the manager and Python APIs; fix active-bitset
  copy lengths and first-call validation of canonical output buffers.
- Enable directed search in Python, reject mixed graph kinds in isomorphism
  tests, validate edge endpoints, and handle duplicate nodes/edges consistently.
- Make Python colored keys invariant to node order and include fixed color
  labels and directedness. **Breaking:** canonical byte keys use a new versioned
  format; stored keys must be regenerated. Color labels are stable built-in
  values (including tuples/frozensets), rather than arbitrary custom objects.
- **Breaking:** search/refinement kernels, algorithm options/contexts and
  Schreier machinery are private. Use `CanonautManager` and its
  builders; raw adjacency remains supported.
- Remove colored-call sorting allocations, redundant dimension inference and
  setup, and quadratic Python result-list copying. Size canonical output
  correctly at construction.
- Complete optional Schreier stabilizer pruning, including fixed-point base
  changes, rooted transversals, generator retention and reusable storage. Expose
  it through `CanonautManager::with_schreier()`; keep it disabled by default.
- Accelerate adjacency permutation and automorphism checks by clearing the
  lowest set bit instead of recomputing a shift to clear the highest one.
- Add optional C-nauty differential tests and configurable nauty include/library
  directories. See `repro/optimization_review.md` for measured performance.

### Added

- `test_data/smoke/`: a committed 1 965-graph cross-validation corpus with
  pynauty reference canonical forms, so a plain `cargo test` on a fresh clone
  verifies canonical forms against C nauty without running `./setup.sh`
  (`tools/gen_smoke_corpus.py`).
- `repro/`: scripts that regenerate every number and figure in the paper
  (`collect_criterion.py`, `make_tables.py`, `plot_figures.py`,
  `measure_memory.sh`) together with the collected measurement data.
- `src/bin/memcheck` reports heap traffic through a counting global allocator,
  separating set-up allocations from steady-state ones.
- Continuous integration (`.github/workflows/ci.yml`): tests on Linux/macOS/
  Windows, MSRV 1.85, clippy, Miri, the full pynauty cross-validation, and the
  Python bindings.
- `CITATION.cff`.
- Refinement-blind graph generators (`canonaut::generators`): Cai-Fürer-Immerman
  graphs over the 3-regular circulant (twisted and untwisted), Sylvester
  Hadamard bipartite graphs, and Paley graphs — plus a `hard` benchmark group
  and tests that use them as ground truth (twisted and untwisted CFI graphs
  must receive different canonical forms; |Aut| of the Paley graph on q
  vertices must be q(q-1)/2).
- `tests::test_word_boundary_paths`: fast canonicalizations covering
  m = 1, 2, 3, 4 words per vertex — the subset the Miri CI job runs.

### Changed

- Graph generators moved from `canonaut::tests` to `canonaut::generators`; the
  test suite and its fixtures are now `#[cfg(test)]` and no longer part of the
  published library.
- Corpus-driven tests skip with a message instead of failing when the
  (uncommitted) full corpora have not been regenerated.

### Fixed

- Short-prune condition in the automorphism-processing branches compared
  `gca_canonical_level` with itself, so `needs_short_prune` was never set. It
  now compares against `gca_first_level`, matching nauty's
  `if (gca_canon != gca_first) needshortprune = TRUE`.
- `fmperm_safe` allocated an `n`-byte visited array on every call for n > 64,
  i.e. once per automorphism found; it now uses a thread-local scratch buffer,
  making steady-state canonicalization allocation-free at every tested size.

### Removed

- Dead C-translation leftovers `itos`, `writeperm`, `fmperm` and `alloc_error`
  (all `pub unsafe fn` with no callers). No `unsafe fn` remains in the public
  API.
- Two no-op `std::mem::transmute` calls on `&[u32]` slices.

## [1.0.0]

- First release: dense canonical labeling, automorphism group generators,
  orbits via the Schreier–Sims chain, vertex-coloured canonicalization,
  digraphs, graph6 I/O, and the `canonaut` Python bindings.
