<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)"
            srcset="https://raw.githubusercontent.com/um-univie/canonaut/v1.0.0/assets/logo-dark.png">
    <img src="https://raw.githubusercontent.com/um-univie/canonaut/v1.0.0/assets/logo.png"
         alt="canonaut" width="380">
  </picture>
</p>

<p align="center">
  <a href="https://crates.io/crates/canonaut"><img alt="crates.io" src="https://img.shields.io/crates/v/canonaut.svg"></a>
  <a href="https://docs.rs/canonaut"><img alt="docs.rs" src="https://img.shields.io/docsrs/canonaut"></a>
  <a href="https://pypi.org/project/canonaut/"><img alt="PyPI" src="https://img.shields.io/pypi/v/canonaut"></a>
  <a href="https://github.com/um-univie/canonaut/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/um-univie/canonaut/actions/workflows/ci.yml/badge.svg"></a>
  <img alt="License" src="https://img.shields.io/badge/license-Apache--2.0-blue.svg">
</p>

A safe, native Rust implementation of [nauty](https://pallini.di.uniroma1.it/) — McKay's algorithm for computing graph automorphism groups and canonical labelings.

Given a graph, `canonaut` finds:
- **Canonical form** — a unique representative of the isomorphism class; two graphs are isomorphic iff their canonical forms are identical
- **Automorphism group** — generators of the symmetry group, plus orbit information and group size

## `DenseGraph` — building a graph

`DenseGraph` stores adjacency as a flat bit-matrix and tracks its own directedness and (optional) vertex coloring.

```rust
use canonaut::structs::DenseGraph;

let mut graph = DenseGraph::new(4);
// 4-cycle: 0-1-2-3-0
graph.add_edge(0, 1);
graph.add_edge(1, 2);
graph.add_edge(2, 3);
graph.add_edge(3, 0);
```

A graph's directedness is fixed at construction: `DenseGraph::new` / `add_edge` for
undirected graphs, `DenseGraph::new_directed` / `add_arc` for digraphs. Calling the
wrong one of `add_edge`/`add_arc` panics rather than silently producing a malformed
adjacency matrix.

## `CanonautManager` — high-level API

`CanonautManager` owns all working buffers and reuses them across calls. Repeated calls, including colored calls, allocate nothing once buffer capacities and search depth have been warmed up. Allocation counts are checked by `src/bin/memcheck` and the allocation regression test; see `repro/`.

### Canonicalization

```rust
use canonaut::structs::{CanonautManager, DenseGraph};

fn main() {
    let n = 4;
    let mut g = DenseGraph::new(n);
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    g.add_edge(2, 3);
    g.add_edge(3, 0);

    let mut manager = CanonautManager::new(n).with_canonization();
    manager.canonize_graph(&g);

    println!("Orbits: {}", manager.statistics().number_of_orbits);
    println!("Canon:  {:?}", manager.canonical_graph);
}
```

`canonize_graph` reads `graph.is_directed()` and configures digraph mode automatically —
no separate builder call needed for digraphs.

### Automorphism group info

After `canonize_graph`, `statistics()` reports the group:

```rust
let stats = manager.statistics();
// Group size = stats.group_size1 * 10^stats.group_size2
println!("group size: {}e{}", stats.group_size1, stats.group_size2);
println!("orbits:     {}", stats.number_of_orbits);
println!("generators: {}", stats.number_of_generators);
```

### Batch canonicalization (reuse manager)

Buffers are only resized when the graph order increases:

```rust
let mut manager = CanonautManager::new(max_n).with_canonization();
for g in graph_database.iter() {
    manager.canonize_graph(g); // g: &DenseGraph
    // use manager.canonical_graph, .labeling(), .orbits() ...
}
```

### Isomorphism testing

Two graphs are isomorphic iff their canonical forms coincide:

```rust
let canon_g = {
    let mut mgr = CanonautManager::new(n).with_canonization();
    mgr.canonize_graph(&g);
    mgr.canonical_graph.clone()
};
let canon_h = {
    let mut mgr = CanonautManager::new(n).with_canonization();
    mgr.canonize_graph(&h);
    mgr.canonical_graph.clone()
};
assert_eq!(canon_g, canon_h); // isomorphic iff equal
```

### Vertex-colored (partitioned) canonicalization

Color values are fixed labels across graphs. Coloring constrains the search: vertices of different colors are never identified by
an automorphism, and the canonical form respects the coloring. Attach a coloring to
the graph itself with `set_colors` (picked up automatically by `canonize_graph`), or
pass colors explicitly per call with `canonize_with_colors` / `canonize_graph_with_colors`
to canonicalize the same graph under several different colorings.

```rust
let mut path_p4 = DenseGraph::new(4); // 0-1-2-3
path_p4.add_edge(0, 1);
path_p4.add_edge(1, 2);
path_p4.add_edge(2, 3);

let mut manager = CanonautManager::new(4).with_canonization();
manager.canonize_graph(&path_p4);
// uncolored: 2 orbits, |Aut| = 2

path_p4.set_colors(vec![0, 0, 1, 1]);
manager.canonize_graph(&path_p4);
// colored: 4 orbits, |Aut| = 1 (picked up automatically)
```

For colored isomorphism testing, compare both canonical adjacency and the
canonical color sequence `colors[manager.labeling()[i] as usize]`. Adjacency alone
does not encode coloring.

### Directed graphs

```rust
let mut g = DenseGraph::new_directed(4);
g.add_arc(0, 1);
g.add_arc(1, 2);
g.add_arc(2, 3);
g.add_arc(3, 0);

let mut manager = CanonautManager::new(4).with_canonization();
manager.canonize_graph(&g); // digraph mode set automatically
println!("|Aut| = {}", manager.statistics().group_size1); // 4.0
```

### Automorphism generators

Register a read-only callback with nauty's `userautomproc` argument order to receive
every generator permutation found. Both slices have length `n`:

```rust
fn on_automorphism(generator_number: u32, perm: &[u32], orbits: &[u32], num_orbits: u32, stabvert: u32, n: u32) {
    println!("{:?}", perm);
}

let mut manager = CanonautManager::new(4)
    .with_canonization()
    .with_automorphism_callback(on_automorphism);
```

## Low-level API

Use `canonaut::structs::add_one_edge`/`add_one_arc` and `CanonautManager::canonize`
directly when working with raw adjacency slices instead of `DenseGraph` — e.g. when
integrating with `graph6` I/O (`canonaut::graph6`) or other external representations.

```rust
use canonaut::structs::{CanonautManager, add_one_edge};
use canonaut::utilities::words_needed;

let n = 6;
let m = words_needed(n) as u32;

// Build a cycle C6
let mut g = vec![0usize; n * m as usize];
for i in 0..n {
    add_one_edge(&mut g, i, (i + 1) % n, m);
}

let mut manager = CanonautManager::new(n).with_canonization();
manager.canonize(&g);
let canon = manager.take_canonical_graph();
println!("{:?}", canon);
```

Raw adjacency is checked for valid dimensions and zero padding before search.
For repeated canonicalizations of the same graph, borrow a `ValidatedGraph`
once to avoid repeating validation or inferring dimensions on each call:

```rust
use canonaut::structs::{CanonautManager, ValidatedGraph};

let adjacency = canonaut::generators::generate_cycle_graph(65);
let graph = ValidatedGraph::new(&adjacency); // no copy or allocation
let mut manager = CanonautManager::new(65).with_canonization();
manager.canonize_validated(&graph);
let colors: Vec<u32> = (0..65).map(|v| v % 3).collect();
manager.canonize_validated_with_colors(&graph, &colors);
```

The view holds an immutable borrow, so its adjacency cannot change while it is
in use. Configure directedness on the manager with `with_digraph()`, as with
raw slices; this is also required for graphs with loops. Coloring is supplied
per call. Search kernels do not perform additional validation.

## Key types

| Type | Location | Purpose |
|------|----------|---------|
| `CanonautManager` | `canonaut::structs` | High-level manager; owns all buffers |
| `DenseGraph` | `canonaut::structs` | Dense graph container (adjacency + directedness + coloring) |
| `ValidatedGraph` | `canonaut::structs` | Immutable raw adjacency view; validates once for repeated searches |
| `Statistics` | `canonaut::structs` | Output: group size, orbit count, search tree stats |
| `add_one_edge` | `canonaut::structs` | Add undirected edge to a raw flat adjacency bitmatrix |
| `add_one_arc` | `canonaut::structs` | Add directed edge (arc) to a raw flat adjacency bitmatrix |
| `words_needed` | `canonaut::utilities` | Compute `m = ⌈n/word_bits⌉` |

## Running tests

```sh
cargo test
```

runs the library and regression tests without external tooling, including the
committed `test_data/smoke/` corpus: 1 965 graphs (n = 9…129, sampled across
every corpus below) whose canonical forms are compared byte-for-byte against
pynauty reference output.

The 28 ignored library tests compare against the *full* corpora, which are not
committed (~125 MB) but regenerate deterministically — the `.g6` files are
exact `geng`/`networkx` output, the `.ref` files pynauty ground truth:

```sh
./setup.sh                    # installs nauty + a Python venv, regenerates
                              # test_data/, then runs the full test suite

./setup.sh --skip-testdata    # skip all of the above, just build + `cargo test`

cargo test -- --include-ignored   # once test_data/ exists
```

Tests that need a corpus that isn't present skip with a message instead of
failing. Regenerate the committed sample with
`.venv-testdata/bin/python3 tools/gen_smoke_corpus.py`.

Lints and undefined-behaviour checks, as run in CI:

```sh
cargo clippy --all-targets -- -D warnings
MIRIFLAGS=-Zmiri-disable-isolation cargo +nightly miri test --lib -- \
  --skip g6 --skip large --skip smoke --skip fuzz \
  --skip hard_instances --skip test_graphs --skip test_random --skip test_known
```

The Miri subset still enters every unsafe path (m = 1, 2, 3, 4 words per
vertex) via `tests::test_word_boundary_paths`; it takes about seven minutes,
whereas the corpus sweeps would take days.

## Hard instances

`canonaut::generators` also builds the standard refinement-blind families used
to stress the search rather than the refinement loops — `generate_cfi_graph`
(Cai–Fürer–Immerman, twisted or not), `generate_hadamard_graph` (Sylvester),
and `generate_paley_graph` (strongly regular). They are benchmarked by the
`hard` group and used as test ground truth.

## Reproducing the published measurements

`repro/` regenerates every table and figure in the paper — Criterion timings
with confidence intervals, peak RSS versus C nauty, and the allocation counts
behind the "no allocation in the steady state" claim. See `repro/README.md`.

## Benchmarks

The Criterion benchmarks compare against the C reference implementation. They
are gated behind the `c-nauty-bench` feature, which compiles a small C helper
and links the system libnauty (headers and library expected under
`/usr/local`):

```sh
cargo bench --features c-nauty-bench
```

A plain `cargo build` / `cargo test` has no C dependencies. Set
`NAUTY_INCLUDE_DIR` and `NAUTY_LIB_DIR` for installations outside `/usr/local`.
`cargo test --release --features c-nauty-bench --test c_reference` additionally
compares 68,600 cases against C nauty, including digraphs, colors, group orders,
orbits, and manager reuse across word boundaries.

For local performance builds, `RUSTFLAGS="-C target-cpu=native"` enables the
current CPU's instructions. Keep portable settings for binaries distributed to
other CPUs. See `repro/optimization_review.md` for measured kernel improvements.

The supported Rust API is `CanonautManager`, graph containers, graph6 I/O and
public construction utilities. Search/refinement kernels, dispatch tables and
algorithm contexts are internal. Optional stabilizer pruning is available through
`CanonautManager::with_schreier()`. It uses nauty-style randomized expansion:
undiscovered stabilizer elements reduce pruning, without changing correctness.
It remains off by default because its overhead can outweigh the saved search.
The implementation currently requires a 64-bit target.

## Python bindings

`bindings/python` (package `canonaut` on PyPI) exposes `canonize()` and
`are_isomorphic()` via [pyo3](https://pyo3.rs)/[maturin](https://www.maturin.rs),
with arbitrary hashable node labels and `networkx` interoperability on top of
the same `CanonautManager` core. It's a separate workspace member — not built
by a plain `cargo build`/`cargo test` at the repo root:

```sh
cd bindings/python
maturin develop --release   # or: maturin build --release
```

See `bindings/python/README.md` for the Python-side API.

## Logo

The mark is the Petersen graph, coloured by the vertex orbits `canonaut`
computes for it after individualising a single vertex — the picture is output
of the library it labels. Regenerate the SVGs (light, dark, square mark) with:

```sh
cargo run --example logo
```

## License

Licensed under the [Apache License, Version 2.0](LICENSE.txt).

This crate is a Rust port of [nauty and Traces](https://pallini.di.uniroma1.it/) (version 2.6+), which is also licensed under the Apache License, Version 2.0.

Original copyright holders:
- Brendan McKay, Australian National University (1984–)
- Adolfo Piperno, University of Rome "Sapienza" (2008–)
- Gunnar Brinkmann, University of Ghent (2009–)
- Magma project, University of Sydney
- Sampo Niskanen and Patric Östergård
