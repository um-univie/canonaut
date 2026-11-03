//! Graph canonicalization and automorphism group computation.
//!
//! `canonaut` is a safe Rust port of the [nauty](https://pallini.di.uniroma1.it/)
//! algorithm (McKay 1981; McKay & Piperno 2014) for computing canonical labelings
//! and automorphism groups of graphs.
//!
//! # Quick start
//!
//! ```rust
//! use canonaut::structs::{CanonautManager, DenseGraph};
//!
//! let vertex_count = 4;
//! let mut graph = DenseGraph::new(vertex_count);
//! graph.add_edge(0, 1);
//! graph.add_edge(1, 2);
//! graph.add_edge(2, 3);
//! graph.add_edge(3, 0);
//!
//! let mut manager = CanonautManager::new(vertex_count).with_canonization();
//! manager.canonize_graph(&graph);
//! assert_eq!(manager.statistics().number_of_orbits, 1); // 4-cycle is vertex-transitive
//! ```
//!
//! # Module overview
//!
//! | Module | Purpose |
//! |--------|---------|
//! [`structs`] | Public types: [`structs::CanonautManager`], [`structs::DenseGraph`], [`structs::SparseGraph`], [`structs::Statistics`] |
//! [`graph6`] | graph6 file format parser (sparse6/digraph6 are rejected) |
//! [`utilities`] | Bit-matrix dimensions: [`utilities::words_needed`] |
//! [`structs::bitset`] | [`structs::BitSetOpsMSB0`] trait — MSB-0 bitset operations on `[usize]` |

// The core search routines mirror nauty's C entry points (`nauty()`,
// `refine()`, `targetcell()`, …) argument-for-argument so that the port can be
// diffed against the reference implementation; splitting those signatures into
// parameter structs would obscure that correspondence.  The dispatch vector is
// likewise a direct translation of nauty's `dispatchvec` function pointers.
#![allow(clippy::too_many_arguments, clippy::type_complexity)]

// Internal kernels assume invariants established by the manager. Some unused
// nauty port helpers are retained for reference, but are not part of the API.
#[allow(dead_code)]
pub(crate) mod bit_ops;
pub mod generators;
pub mod graph6;
#[allow(dead_code)]
pub(crate) mod simd_ops;
pub(crate) mod consts;
pub mod io;
#[allow(dead_code)]
pub(crate) mod refine;
#[allow(dead_code)]
pub(crate) mod rng;
#[allow(dead_code)]
pub(crate) mod support;
#[allow(dead_code)]
pub(crate) mod search;
pub mod structs;
pub mod utilities;

#[cfg(test)]
mod test_vectors;
#[cfg(test)]
mod tests;
