//! Core types for configuring, running, and inspecting canonaut computations.
//!
//! The primary entry point for most users is [`CanonautManager`], which owns all
//! algorithm buffers and exposes a single [`CanonautManager::canonize`] method.
//! Search options and working contexts are internal implementation details.

use crate::consts::ALGORITHM_INFINITY;
use crate::io::OutputStream;
use crate::refine::DISPATCH_GRAPH;
use crate::refine::GraphContext;
use crate::search::canonicalize;
use crate::structs::{
    graph::{DenseGraph, ValidatedGraph},
    schreier_arena::SchreierSims,
    target_cell::TargetCellNodePool,
};
use crate::utilities::words_needed;

/// All mutable state threaded through a single canonicalization search.
///
/// Replaces the C global variables with a single owned struct. Buffer capacities grow
/// monotonically via [`AlgorithmContext::ensure_capacity`] and are reused across
/// repeated calls to avoid allocation pressure on the hot path.
#[derive(Default)]
#[allow(dead_code)]
pub(crate) struct AlgorithmContext {
    // ---- configuration ----
    /// Compute and store the canonical graph.
    pub get_canonical: bool,
    /// Treat the graph as a digraph (directed).
    pub digraph: bool,
    /// Write automorphism generators to the output stream.
    pub write_automorphisms: bool,
    /// Write partition markers to the output stream.
    pub do_markers: bool,
    /// Use the Cartesian product refinement.
    pub cartesian: bool,
    /// Enable the Schreier–Sims stabilizer chain.
    pub use_schreier_sims: bool,
    /// Maximum line length for automorphism output (`0` = unlimited).
    pub line_length: u32,
    /// Depth limit for target-cell selection (`100` = unlimited).
    pub target_cell_level: i32,
    /// Minimum search-tree depth at which the invariant procedure is called.
    pub min_invariant_level: i32,
    /// Maximum search-tree depth at which the invariant procedure is called.
    pub max_invariant_level: i32,
    /// Extra integer argument forwarded to the invariant procedure.
    pub invariant_argument: u32,
    /// Bump-allocator pool for target-cell node lists.
    pub target_cell_node_pool: TargetCellNodePool,

    // ---- user callbacks ----
    /// Called at each search-tree node before refinement.
    pub user_node_procedure:
        Option<fn(&[usize], &mut [u32], &mut [u32], u32, u32, i32, u32, u32, u32) -> ()>,
    /// Called whenever a new automorphism generator is found.
    pub user_automorphism_procedure: Option<fn(u32, &[u32], &[u32], u32, u32, u32) -> ()>,
    /// Called when the search descends to a new level.
    pub user_level_procedure: Option<
        fn(
            &mut [u32],
            &mut [u32],
            u32,
            &mut [u32],
            &mut Statistics,
            u32,
            u32,
            u32,
            u32,
            u32,
            u32,
        ) -> (),
    >,
    /// Called to override canonical-labeling comparison.
    pub user_canonical_procedure:
        Option<fn(&[usize], &mut [u32], &mut [usize], usize, u32, u32, u32) -> u32>,
    /// Vertex-invariant procedure for pruning the search tree.
    pub invariant_procedure: Option<
        fn(&[usize], &mut [u32], &mut [u32], i32, u32, u32, &mut [u32], u32, bool, u32, u32) -> (),
    >,
    /// Output stream for automorphism/marker output.
    pub outfile: OutputStream,
    /// Function dispatch table (dense vs. sparse graph operations).
    pub dispatch: DispatchVec,

    // ---- graph dimensions ----
    /// Words per adjacency row (`m = ⌈n/64⌉`).
    pub words_per_vertex: u32,
    /// Number of vertices.
    pub number_of_vertices: u32,

    // ---- search-tree state ----
    /// Number of times the invariant procedure was called.
    pub invariant_applications: usize,
    /// Number of times the invariant produced a useful pruning.
    pub invariant_successes: usize,
    /// Deepest level at which the invariant succeeded.
    pub invariant_success_level: i32,
    /// Greatest-common-ancestor level of the first path.
    pub gca_first_level: u32,
    /// Greatest-common-ancestor level of the canonical path.
    pub gca_canonical_level: i32,
    /// Shallowest level where the partition is non-trivial.
    pub first_nontrivial_level: u32,
    /// Depth at which the partition became uniform.
    pub uniform_partition_depth: u32,
    /// Level at which the first path last matched the canonical path.
    pub first_path_match_level: u32,
    /// Level at which the canonical path was last updated.
    pub canonical_match_level: i32,
    /// Result of the last canonical-path comparison (`-1`, `0`, or `1`).
    pub canonical_comparison: i32,
    /// Number of adjacency rows that matched in the last canonical comparison.
    pub canonical_matching_rows: u32,
    /// Search depth of the current canonical leaf.
    pub canonical_depth: u32,
    /// Representative vertex of the current orbit.
    pub orbit_representative: u32,
    /// Coset vertex used in the short-prune test.
    pub coset_vertex: u32,
    /// Whether a short prune can be applied at the current node.
    pub needs_short_prune: bool,

    // ---- reusable working buffers ----
    /// Main workspace; sized to `2 * 500 * m` words and reused across calls.
    pub workspace: Vec<usize>,
    /// Current high-water mark in `workspace`.
    pub workspace_top: usize,
    /// Logical free pointer within `workspace`.
    pub free_memory_pointer: usize,
    /// Schreier–Sims stabilizer chain.
    pub schreiersims: SchreierSims,

    /// Small secondary workspace (`2 * m` words).
    pub default_workspace: Vec<usize>,
    /// Scratch permutation buffer (`n` elements).
    pub work_permutation: Vec<u32>,
    /// Fixed-point bitset (`m` words).
    pub fixed_points: Vec<usize>,
    /// Active-cells bitset (`m` words).
    pub active: Vec<usize>,
    /// Labeling of the first leaf in the search tree (`n` elements).
    pub first_labeling: Vec<u32>,
    /// Labeling of the current canonical leaf (`n` elements).
    pub canonical_labeling: Vec<u32>,
    /// Invariant codes along the first path (`n + 2` entries).
    pub first_path_codes: Vec<i16>,
    /// Invariant codes along the canonical path (`n + 2` entries).
    pub canonical_path_codes: Vec<i16>,
    /// Target-cell indices along the first path (`n + 2` entries).
    pub first_target_cell: Vec<i32>,
}

#[allow(dead_code)]
impl AlgorithmContext {
    /// Creates a zeroed context pre-sized for `number_of_vertices` vertices.
    pub fn new(number_of_vertices: u32) -> Self {
        Self {
            get_canonical: false,
            digraph: false,
            write_automorphisms: false,
            do_markers: false,
            cartesian: false,
            use_schreier_sims: false,
            line_length: 0,
            target_cell_level: 0,
            min_invariant_level: 0,
            max_invariant_level: 0,
            invariant_argument: 0,
            user_node_procedure: None,
            user_automorphism_procedure: None,
            user_level_procedure: None,
            user_canonical_procedure: None,
            invariant_procedure: None,
            outfile: OutputStream::StandardOutput,
            dispatch: DispatchVec {
                is_automorphism: None,
                test_canonical_labeling: None,
                update_canonical_graph: None,
                refine: None,
                refine1: None,
                cheap_automorphism_check: None,
                target_cell: None,
                freedyn: None,
                init: None,
                cleanup: None,
            },
            words_per_vertex: 0,
            number_of_vertices: 0,
            schreiersims: SchreierSims::new(number_of_vertices as usize),
            ..Default::default()
        }
    }

    /// Ensure all internal algorithm buffers have sufficient capacity for the given
    /// graph dimensions. Allocated capacities are retained across calls; the
    /// active bitset also tracks its current logical length.
    pub fn ensure_capacity(&mut self, words_per_vertex: u32, number_of_vertices: u32) {
        let word_count = words_per_vertex as usize;
        let vertex_count = number_of_vertices as usize;
        // All vertex-sized buffers grow together here; first_path_codes also
        // includes the required two slack entries. The active set alone can
        // shrink, and must retain the current graph's exact word count.
        if self.first_path_codes.len() >= vertex_count + 2 && self.active.len() == word_count {
            debug_assert!(self.work_permutation.len() >= vertex_count);
            debug_assert!(self.first_labeling.len() >= vertex_count);
            debug_assert!(self.canonical_labeling.len() >= vertex_count);
            debug_assert!(self.canonical_path_codes.len() >= vertex_count + 2);
            debug_assert!(self.first_target_cell.len() >= vertex_count + 2);
            debug_assert!(self.fixed_points.len() >= word_count);
            debug_assert!(self.workspace.len() >= 2 * 500 * word_count);
            return;
        }
        if self.work_permutation.len() < vertex_count {
            self.work_permutation.resize(vertex_count, 0);
        }
        if self.default_workspace.len() < 2 * word_count {
            self.default_workspace.resize(2 * word_count, 0);
        }
        if self.fixed_points.len() < word_count {
            self.fixed_points.resize(word_count, 0);
        }
        // Refinement scans the whole active slice. Retain capacity, but never
        // expose stale words from a larger graph as active cells of this one.
        self.active.resize(word_count, 0);
        if self.first_labeling.len() < vertex_count {
            self.first_labeling.resize(vertex_count, 0);
        }
        if self.canonical_labeling.len() < vertex_count {
            self.canonical_labeling.resize(vertex_count, 0);
        }
        if self.first_path_codes.len() < vertex_count + 2 {
            self.first_path_codes.resize(vertex_count + 2, 0);
        }
        if self.canonical_path_codes.len() < vertex_count + 2 {
            self.canonical_path_codes.resize(vertex_count + 2, 0);
        }
        if self.first_target_cell.len() < vertex_count + 2 {
            self.first_target_cell.resize(vertex_count + 2, 0);
        }
        // Pre-size workspace so it can be reused across calls (default worksize = 500)
        let workspace_min = 2 * 500 * word_count;
        if self.workspace.len() < workspace_min {
            self.workspace.resize(workspace_min, 0);
        }
    }
}

/// Function dispatch table that selects dense- or sparse-graph implementations.
///
/// The default ([`crate::refine::DISPATCH_GRAPH`]) uses the
/// dense bit-matrix routines. Replacing individual slots enables custom
/// refinement, automorphism-check, or target-cell strategies.
#[derive(Copy, Clone, Default)]
#[repr(C)]
#[allow(dead_code)]
pub(crate) struct DispatchVec {
    pub is_automorphism: Option<fn(&[usize], &[u32], bool, u32, u32) -> bool>,
    pub test_canonical_labeling:
        Option<fn(&[usize], &[usize], &mut [u32], &mut u32, u32, u32, &mut GraphContext) -> i32>,
    pub update_canonical_graph:
        Option<fn(&[usize], &mut [usize], &[u32], u32, u32, u32, &mut GraphContext) -> ()>,
    pub refine: Option<
        fn(
            &[usize],
            &mut [u32],
            &mut [u32],
            i32,
            &mut u32,
            &mut [u32],
            &mut [usize],
            &mut u32,
            u32,
            u32,
            &mut GraphContext,
        ) -> (),
    >,
    pub refine1: Option<
        fn(
            &[usize],
            &mut [u32],
            &mut [u32],
            i32,
            &mut u32,
            &mut [u32],
            &mut [usize],
            &mut u32,
            u32,
            u32,
            &mut GraphContext,
        ) -> (),
    >,
    pub cheap_automorphism_check: Option<fn(&[u32], u32, bool, u32) -> bool>,
    pub target_cell: Option<
        fn(&[usize], &[u32], &[u32], i32, i32, bool, i32, u32, u32, &mut GraphContext) -> i32,
    >,
    pub freedyn: Option<fn() -> ()>,
    pub init: Option<
        unsafe fn(
            &mut [usize],
            &mut [usize],
            &mut [usize],
            &mut [usize],
            *mut u32,
            *mut u32,
            *mut usize,
            *mut CanonautOptions,
            *mut u32,
            u32,
            u32,
        ) -> (),
    >,
    pub cleanup: Option<
        unsafe fn(
            &mut [usize],
            &mut [usize],
            Option<&mut [usize]>,
            &mut [usize],
            *mut u32,
            *mut u32,
            *mut CanonautOptions,
            *mut Statistics,
            u32,
            u32,
        ) -> (),
    >,
}

/// Counters and measurements collected during a canonicalization run.
///
/// Returned by [`CanonautManager::statistics`] after [`CanonautManager::canonize`].
#[derive(Copy, Clone, Debug, Default)]
#[repr(C)]
pub struct Statistics {
    /// Floating-point mantissa of the automorphism group order
    /// (`group_size ≈ group_size1 × 10^group_size2`).
    pub group_size1: f64,
    /// Base-10 exponent of the automorphism group order.
    pub group_size2: u32,
    /// Number of vertex orbits under the automorphism group.
    pub number_of_orbits: u32,
    /// Number of generator permutations found.
    pub number_of_generators: u32,
    /// Non-zero on error; see nauty documentation for error codes.
    pub error_status_code: u32,
    /// Total search-tree nodes visited.
    pub number_of_nodes: usize,
    /// Leaf nodes pruned as non-canonical (bad leaves).
    pub number_of_bad_leaves: usize,
    /// Maximum depth reached in the search tree.
    pub maximum_search_depth: u32,
    /// Sum of target-cell sizes over all nodes visited.
    pub total_target_cell_members: usize,
    /// Number of times the canonical graph was updated.
    pub canonical_updates: usize,
    /// Total calls to the vertex-invariant procedure.
    pub invariant_procedure_calls: usize,
    /// Calls to the invariant that produced a useful pruning.
    pub successful_invariant_calls: usize,
    /// Shallowest search depth at which the invariant succeeded.
    pub shallowest_invariant_success_level: i32,
}

impl Statistics {
    /// Creates a zeroed statistics block.
    pub fn new() -> Self {
        Statistics {
            group_size1: 0.,
            group_size2: 0,
            number_of_orbits: 0,
            number_of_generators: 0,
            error_status_code: 0,
            number_of_nodes: 0,
            number_of_bad_leaves: 0,
            maximum_search_depth: 0,
            total_target_cell_members: 0,
            canonical_updates: 0,
            invariant_procedure_calls: 0,
            successful_invariant_calls: 0,
            shallowest_invariant_success_level: 0,
        }
    }
}

/// Full configuration block passed to the core `canonicalize()` function.
///
/// Configured through [`CanonautManager`] builders. Defaults to dense dispatch,
/// no callbacks and `line_length = 78`.
#[allow(dead_code)]
pub(crate) struct CanonautOptions {
    /// Compute and store the canonical graph after each call.
    pub get_canonical: bool,
    /// Treat the input as a directed graph.
    pub digraph: bool,
    /// Write automorphism generators to `outfile`.
    pub write_automorphisms: bool,
    /// Write partition markers to `outfile`.
    pub write_markers: bool,
    /// Start with the trivial (all-in-one-cell) partition.
    pub default_partition: bool,
    /// Use the Cartesian product refinement heuristic.
    pub cartesian: bool,
    /// Maximum output line length in characters (`0` = unlimited).
    pub line_length: u32,
    /// Output stream for automorphism/marker text (`None` = stdout).
    pub outfile: Option<OutputStream>,
    /// Replaces the built-in partition refinement with a custom procedure.
    pub user_refinement_procedure: Option<
        fn(
            &[usize],
            &mut [u32],
            &mut [u32],
            i32,
            &mut u32,
            &mut [u32],
            &mut [usize],
            &mut u32,
            u32,
            u32,
            &mut GraphContext,
        ) -> (),
    >,
    pub user_automorphism_procedure: Option<fn(u32, &[u32], &[u32], u32, u32, u32) -> ()>,
    pub user_level_procedure: Option<
        fn(
            &mut [u32],
            &mut [u32],
            u32,
            &mut [u32],
            &mut Statistics,
            u32,
            u32,
            u32,
            u32,
            u32,
            u32,
        ) -> (),
    >,
    pub user_node_procedure:
        Option<fn(&[usize], &mut [u32], &mut [u32], u32, u32, i32, u32, u32, u32) -> ()>,
    pub user_canonical_procedure:
        Option<fn(&[usize], &mut [u32], &mut [usize], usize, u32, u32, u32) -> u32>,
    pub invariant_procedure: Option<
        fn(&[usize], &mut [u32], &mut [u32], i32, u32, u32, &mut [u32], u32, bool, u32, u32) -> (),
    >,
    pub target_cell_level: i32,
    pub min_invariant_level: i32,
    pub max_invariant_level: i32,
    pub invariant_argument: u32,
    pub dispatch: DispatchVec,
    pub schreier: bool,
    pub extra_options: *mut std::ffi::c_void,
    pub algorithm_context: AlgorithmContext,
}

impl Default for CanonautOptions {
    fn default() -> Self {
        CanonautOptions {
            get_canonical: false,
            digraph: false,
            write_automorphisms: false,
            write_markers: false,
            default_partition: true,
            cartesian: false,
            line_length: 78,
            outfile: None,
            user_refinement_procedure: None,
            user_automorphism_procedure: None,
            user_level_procedure: None,
            user_node_procedure: None,
            user_canonical_procedure: None,
            invariant_procedure: None,
            target_cell_level: 100,
            min_invariant_level: 0,
            max_invariant_level: 1,
            invariant_argument: 0,
            dispatch: DISPATCH_GRAPH,
            schreier: false,
            extra_options: std::ptr::null() as *const std::os::raw::c_void
                as *mut std::os::raw::c_void,
            algorithm_context: AlgorithmContext::default(),
        }
    }
}

/// High-level owner of all algorithm buffers.
///
/// Pre-allocates output buffers at construction and reuses them across repeated
/// [`canonize`](CanonautManager::canonize) calls, so no allocation occurs on the
/// hot path once buffer capacity and search depth have been warmed up.
///
/// # Example
///
/// ```rust
/// use canonaut::structs::{CanonautManager, add_one_edge};
/// use canonaut::utilities::words_needed;
///
/// let vertex_count = 6;
/// let words_per_vertex = words_needed(vertex_count) as u32;
/// let mut adjacency = vec![0usize; vertex_count * words_per_vertex as usize];
/// // Build a 6-cycle
/// for vertex in 0..vertex_count { add_one_edge(&mut adjacency, vertex, (vertex + 1) % vertex_count, words_per_vertex); }
///
/// let mut manager = CanonautManager::new(vertex_count).with_canonization();
/// manager.canonize(&adjacency);
/// assert_eq!(manager.statistics().number_of_orbits, 1);
/// ```
pub struct CanonautManager {
    /// The canonical form of the last graph passed to [`canonize`](CanonautManager::canonize).
    ///
    /// Valid only when [`with_canonization`](CanonautManager::with_canonization) was called
    /// before the first run. Length is `n * m` words after the call. For colored
    /// isomorphism keys, also compare `colors[labeling[i]]` in canonical order;
    /// this adjacency matrix alone does not encode vertex colors.
    pub canonical_graph: Vec<usize>,
    labeling: Vec<u32>,
    partition: Vec<u32>,
    orbits: Vec<u32>,
    active: Option<Vec<usize>>,
    options: CanonautOptions,
    statistics: Statistics,
    graph_context: GraphContext,
}

impl CanonautManager {
    /// Creates a manager pre-sized for graphs on `number_of_vertices` vertices.
    ///
    /// Output and refinement buffers are allocated here. Search-tree nodes are
    /// allocated lazily and reused on subsequent calls.
    pub fn new(number_of_vertices: usize) -> Self {
        let mut options = CanonautOptions::default();
        // Pre-size all algorithm buffers for the expected graph dimensions so the
        // first call doesn't pay for allocation on the hot path.
        let word_count = words_needed(number_of_vertices);
        options
            .algorithm_context
            .ensure_capacity(word_count as u32, number_of_vertices as u32);

        let mut graph_context = GraphContext::new(number_of_vertices);
        // Ensure bucket has the +2 slack that ensure_bucket_size expects.
        graph_context.ensure_bucket_size(number_of_vertices as u32);

        CanonautManager {
            canonical_graph: vec![0; number_of_vertices * word_count],
            labeling: vec![0; number_of_vertices],
            partition: vec![0; number_of_vertices],
            orbits: vec![0; number_of_vertices],
            active: None,
            options,
            statistics: Statistics::new(),
            graph_context,
        }
    }

    /// Enables canonical-graph computation.
    ///
    /// After each [`canonize`](Self::canonize) call the result is stored in
    /// [`canonical_graph`](Self::canonical_graph).
    pub fn with_canonization(mut self) -> Self {
        self.options.get_canonical = true;
        self
    }

    /// Enable optional randomized stabilizer pruning for highly symmetric graphs.
    /// Missing random generators only reduces pruning; results remain exact.
    pub fn with_schreier(mut self) -> Self {
        self.options.schreier = true;
        self
    }

    /// Treats input graphs as digraphs: row `v`, bit `w` means the arc `v → w`.
    /// The adjacency matrix is not assumed to be symmetric.
    pub fn with_digraph(mut self) -> Self {
        self.options.digraph = true;
        self
    }

    /// Registers a callback invoked for every automorphism generator found.
    ///
    /// The arguments follow nauty's `userautomproc`:
    /// `(generator_number, perm, orbits, num_orbits, stabvert, n)`.
    /// The callback borrows read-only search state; both slices have length `n`.
    pub fn with_automorphism_callback(
        mut self,
        callback: fn(u32, &[u32], &[u32], u32, u32, u32),
    ) -> Self {
        self.options.user_automorphism_procedure = Some(callback);
        self
    }

    /// No-op; retained for API compatibility.
    pub fn preallocate(&mut self, _max_number_of_vertices: usize) {}

    /// Grows output buffers to fit a graph of `number_of_vertices` × `words_per_vertex`.
    ///
    /// Called automatically by [`canonize`](Self::canonize); callers may also
    /// use this to pre-size buffers before a run.
    pub fn prepare_for_run(&mut self, number_of_vertices: usize, words_per_vertex: usize) {
        // Logical lengths follow the current graph; allocated capacities are retained.
        let graph_words = number_of_vertices * words_per_vertex;
        Self::ensure_len(&mut self.canonical_graph, graph_words);
        // These private buffers always have the same logical length. Only the
        // public canonical_graph can be independently resized or taken out.
        if self.labeling.len() != number_of_vertices {
            Self::ensure_len(&mut self.labeling, number_of_vertices);
            Self::ensure_len(&mut self.partition, number_of_vertices);
            Self::ensure_len(&mut self.orbits, number_of_vertices);
        }
    }

    #[inline]
    fn ensure_len<T: Default>(vec: &mut Vec<T>, target_len: usize) {
        if vec.len() < target_len {
            vec.resize_with(target_len, T::default);
        }
        vec.truncate(target_len);
    }

    /// Clears and resizes `to_fill` to exactly `size` default-initialized elements.
    pub fn resize_and_fill<T>(to_fill: &mut Vec<T>, size: usize)
    where
        T: Default + Copy,
    {
        to_fill.clear();
        to_fill.resize(size, T::default());
    }

    /// Moves `canonical_graph` out of the manager, leaving an empty `Vec` in its place.
    ///
    /// Useful when the caller wants to own the result without cloning. A subsequent
    /// [`canonize`](Self::canonize) will reallocate the buffer as needed.
    pub fn take_canonical_graph(&mut self) -> Vec<usize> {
        std::mem::take(&mut self.canonical_graph)
    }

    /// Returns a reference to the statistics collected during the last [`canonize`](Self::canonize) call.
    pub fn statistics(&self) -> &Statistics {
        &self.statistics
    }

    /// Returns the vertex labeling from the last [`canonize`](Self::canonize) call.
    ///
    /// `labeling[i]` is the original vertex index of position `i` in the canonical order.
    pub fn labeling(&self) -> &[u32] {
        &self.labeling
    }

    /// Returns the orbit partition from the last [`canonize`](Self::canonize) call.
    ///
    /// `orbits[v]` is the canonical representative of vertex `v`'s orbit.
    pub fn orbits(&self) -> &[u32] {
        &self.orbits
    }

    /// Runs the canonicalization algorithm on the given adjacency slice.
    ///
    /// `graph` must be a flat bit-matrix with row stride `m = words_needed(n)` words,
    /// using MSB-0 bit ordering (bit `j` of row `i` set ↔ edge `i→j`).  The length
    /// `graph.len()` must equal `n * m` for some valid `(n, m)` pair, including
    /// an empty slice for the zero-vertex graph. Padding bits must be zero.
    /// Undirected input must be symmetric and loopless; use [`with_digraph`](Self::with_digraph)
    /// for directed graphs or graphs containing loops.
    ///
    /// After the call:
    /// - [`canonical_graph`](Self::canonical_graph) contains the canonical adjacency
    ///   matrix (if [`with_canonization`](Self::with_canonization) was set).
    /// - [`labeling`](Self::labeling) and [`orbits`](Self::orbits) are updated.
    /// - [`statistics`](Self::statistics) reflects this run.
    ///
    /// # Panics
    ///
    /// Panics if the dimensions are invalid or unsupported, or padding is nonzero.
    /// For repeated calls on the same adjacency, use [`ValidatedGraph`] and
    /// [`canonize_validated`](Self::canonize_validated) to validate only once.
    pub fn canonize(&mut self, graph: &[usize]) {
        self.canonize_validated(&ValidatedGraph::new(graph));
    }

    /// Like [`canonize`](Self::canonize), using a previously validated immutable
    /// adjacency view. No input scan or dimension inference is repeated.
    pub fn canonize_validated(&mut self, graph: &ValidatedGraph<'_>) {
        self.run(graph, None);
    }

    /// Like [`canonize_with_colors`](Self::canonize_with_colors), without
    /// revalidating the adjacency. `colors.len()` must equal the graph order.
    pub fn canonize_validated_with_colors(&mut self, graph: &ValidatedGraph<'_>, colors: &[u32]) {
        self.run(graph, Some(colors));
    }

    fn run(
        &mut self,
        graph: &ValidatedGraph<'_>,
        colors: Option<&[u32]>,
    ) {
        // Dimensions and padding were checked at the public input boundary.
        let number_of_vertices = graph.number_of_vertices;
        let words_per_vertex = graph.words_per_vertex;
        self.prepare_for_run(number_of_vertices, words_per_vertex);
        self.options.default_partition = colors.is_none();
        if let Some(colors) = colors {
            assert_eq!(
                colors.len(),
                number_of_vertices,
                "colors.len() must equal the vertex count"
            );
            for (vertex, label) in self.labeling.iter_mut().enumerate() {
                *label = vertex as u32;
            }
            // The vertex tie-breaker makes unstable sorting deterministic, and
            // using the labeling buffer avoids allocating on every colored call.
            self.labeling
                .sort_unstable_by_key(|&vertex| (colors[vertex as usize], vertex));
            for position in 0..number_of_vertices {
                let last_in_cell = position + 1 == number_of_vertices
                    || colors[self.labeling[position] as usize]
                        != colors[self.labeling[position + 1] as usize];
                self.partition[position] = if last_in_cell { 0 } else { ALGORITHM_INFINITY };
            }
        }

        let canonical_graph = if self.options.get_canonical {
            Some(self.canonical_graph.as_mut_slice())
        } else {
            None
        };
        let active = self.active.as_deref();

        // Pass an empty workspace slice — algorithm_context.workspace is pre-sized by
        // ensure_capacity (called inside canonicalize) and reused across calls.
        canonicalize(
            graph.adjacency,
            &mut self.labeling,
            &mut self.partition,
            active,
            &mut self.orbits,
            &mut self.options,
            &mut self.statistics,
            &mut [],
            words_per_vertex as u32,
            number_of_vertices as u32,
            canonical_graph,
            &mut self.graph_context,
        );
    }

    /// Like [`canonize`](Self::canonize), but takes a [`DenseGraph`] directly
    /// instead of a raw adjacency slice.  Digraph mode is set from
    /// [`graph.is_directed()`](DenseGraph::is_directed); no need to call
    /// [`with_digraph`](Self::with_digraph) yourself.  If the graph has a coloring
    /// attached via [`DenseGraph::set_colors`], it's used automatically (equivalent
    /// to calling [`canonize_with_colors`](Self::canonize_with_colors)); otherwise
    /// this is a plain uncolored [`canonize`](Self::canonize).
    pub fn canonize_graph(&mut self, graph: &DenseGraph) {
        let validated = ValidatedGraph::with_order(&graph.adjacency, graph.number_of_vertices());
        assert_eq!(graph.words_per_vertex as usize, validated.words_per_vertex,
            "DenseGraph row stride does not match the vertex count");
        self.options.digraph = graph.is_directed() || graph.has_loops();
        self.run(&validated, graph.colors());
    }

    /// Like [`canonize_graph`](Self::canonize_graph), but with `colors` passed
    /// explicitly — these are used instead of, and regardless of, any coloring
    /// attached to `graph` via [`DenseGraph::set_colors`].  Use this to canonicalize
    /// the same graph under several different colorings; see
    /// [`canonize_with_colors`](Self::canonize_with_colors).
    pub fn canonize_graph_with_colors(&mut self, graph: &DenseGraph, colors: &[u32]) {
        let validated = ValidatedGraph::with_order(&graph.adjacency, graph.number_of_vertices());
        assert_eq!(graph.words_per_vertex as usize, validated.words_per_vertex,
            "DenseGraph row stride does not match the vertex count");
        self.options.digraph = graph.is_directed() || graph.has_loops();
        self.run(&validated, Some(colors));
    }

    /// Like [`canonize`](Self::canonize), but with a vertex coloring (partition).
    ///
    /// `colors[v]` is an arbitrary `u32` color for vertex `v`; colors need not be
    /// contiguous or sorted. Vertices with different colors are never mapped onto
    /// each other by any automorphism, and the canonical form respects the coloring.
    /// `colors.len()` must equal the vertex count implied by `graph.len()`.
    ///
    /// Only this call is colored — the next plain [`canonize`](Self::canonize) call
    /// goes back to the trivial (uncolored) partition.
    ///
    /// # Panics
    ///
    /// Panics if `graph.len()` does not correspond to a valid `(n, m)` pair, or if
    /// `colors.len()` does not equal `n`, padding is nonzero, or the order is unsupported.
    pub fn canonize_with_colors(&mut self, graph: &[usize], colors: &[u32]) {
        assert_eq!(
            graph.len(),
            colors.len() * words_needed(colors.len()),
            "colors.len() must equal the vertex count implied by graph.len()"
        );
        let validated = ValidatedGraph::with_order(graph, colors.len());
        self.canonize_validated_with_colors(&validated, colors);
    }
}
