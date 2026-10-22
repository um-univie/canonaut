//! Graph representations and adjacency-matrix primitives.
//!
//! Two graph formats are provided:
//! - [`DenseGraph`] — flat bit-matrix, `n × m` `usize` words, MSB-0 bit ordering.
//! - [`SparseGraph`] — compressed adjacency-list (CSR-style).
//!
//! For most use-cases, [`DenseGraph`] is the right choice.  Build it with
//! [`DenseGraph::new`] and [`DenseGraph::add_edge`], then pass it directly to
//! [`CanonautManager::canonize_graph`](crate::structs::CanonautManager::canonize_graph).

use crate::structs::BitSetOpsMSB0;
use crate::utilities::words_needed;

/// An immutable adjacency slice with validated dimensions and zero padding.
///
/// Validate once when importing a raw graph, then reuse this view with
/// [`CanonautManager::canonize_validated`](crate::structs::CanonautManager::canonize_validated)
/// to avoid repeating the input scan. This allocates nothing. The borrow prevents
/// mutation of the underlying adjacency while the view is in use.
/// Directedness is configured on the manager, just as for raw adjacency. This
/// validates storage, not symmetry or the absence of loops: undirected input
/// must be symmetric and loopless; otherwise enable the manager's digraph mode.
///
/// ```
/// use canonaut::structs::{CanonautManager, ValidatedGraph};
/// let adjacency = canonaut::generators::generate_cycle_graph(65);
/// let graph = ValidatedGraph::new(&adjacency);
/// let mut manager = CanonautManager::new(65).with_canonization();
/// manager.canonize_validated(&graph);
/// ```
///
/// ```compile_fail
/// use canonaut::structs::{CanonautManager, ValidatedGraph};
/// let mut adjacency = vec![0usize; 2];
/// let graph = ValidatedGraph::new(&adjacency);
/// adjacency[0] = 1; // cannot invalidate a live view
/// CanonautManager::new(2).canonize_validated(&graph);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct ValidatedGraph<'a> {
    pub(crate) adjacency: &'a [usize],
    pub(crate) number_of_vertices: usize,
    pub(crate) words_per_vertex: usize,
}

impl<'a> ValidatedGraph<'a> {
    /// Borrows a flat MSB-0 adjacency matrix without copying it.
    ///
    /// # Panics
    /// Panics if the length is not `n * ceil(n / 64)` for a supported order,
    /// or if a padding bit outside `0..n` is set.
    #[inline]
    pub fn new(adjacency: &'a [usize]) -> Self {
        let len = adjacency.len();
        if len <= 64 {
            // Single-word rows have n == len, including the empty graph.
            // Avoid square-root inference and division on this common path.
            return Self::with_order(adjacency, len);
        }
        Self::infer_order(adjacency)
    }

    // Keep general dimension inference out of the single-word fast path.
    #[inline(never)]
    fn infer_order(adjacency: &'a [usize]) -> Self {
        let len = adjacency.len();
        // n*m is strictly increasing, with m = ceil(n/64).
        let words = ((len - 1) / 64).isqrt() + 1;
        let vertices = len / words;
        assert!(len % words == 0 && words_needed(vertices) == words,
            "graph slice length does not correspond to a valid (n, m) pair");
        Self::with_order(adjacency, vertices)
    }

    pub(crate) fn with_order(adjacency: &'a [usize], vertices: usize) -> Self {
        assert!(vertices <= (crate::consts::ALGORITHM_INFINITY - 2) as usize,
            "graph order exceeds the supported limit");
        let words = words_needed(vertices);
        assert_eq!(adjacency.len(), vertices * words,
            "graph slice length does not match the vertex count");
        let used = vertices % 64;
        if used != 0 {
            // Reduce first, then test once: the single-word case vectorizes.
            // Multiword graphs only need the last word of each adjacency row.
            let last_words = if words == 1 {
                adjacency.iter().fold(0, |bits, &word| bits | word)
            } else {
                adjacency.chunks_exact(words).fold(0, |bits, row| bits | row[words - 1])
            };
            assert_eq!(last_words & (usize::MAX >> used), 0,
                "graph padding bits must be zero");
        }
        Self { adjacency, number_of_vertices: vertices, words_per_vertex: words }
    }
}

/// Sparse graph in CSR (compressed-sparse-row) style.
///
/// Suitable for large, low-density graphs where the dense bit-matrix would
/// waste significant memory.  The sparse-graph code paths in the core algorithm
/// are currently incomplete; prefer [`DenseGraph`] for production use.
#[derive(Clone, Default)]
#[repr(C)]
pub struct SparseGraph {
    /// Number of directed edge slots (undirected loops counted once).
    pub number_of_directed_edges: usize,
    /// `edge_indices[v]` is the start of vertex `v`'s neighbor list in `adjacency_list`.
    pub edge_indices: Vec<usize>,
    /// Number of vertices.
    pub number_of_vertices: u32,
    /// Out-degree of each vertex.
    pub vertex_out_degrees: Vec<u32>,
    /// Flat neighbor list: neighbors of vertex `v` occupy
    /// `adjacency_list[edge_indices[v] .. edge_indices[v] + vertex_out_degrees[v]]`.
    pub adjacency_list: Vec<u32>,
    /// Optional per-edge integer weights (not yet used by the algorithm).
    pub weights: Option<Vec<i32>>,
}

impl SparseGraph {
    pub fn resize_adjacency_list(&mut self, capacity: usize) {
        if capacity > self.adjacency_list.len() {
            self.adjacency_list.resize(capacity, 0)
        }
    }
    pub fn resize_vertex_out_degrees(&mut self, capacity: usize) {
        if capacity > self.vertex_out_degrees.len() {
            self.vertex_out_degrees.resize(capacity, 0)
        }
    }
    pub fn resize_edge_indices(&mut self, capacity: usize) {
        if capacity > self.edge_indices.len() {
            self.edge_indices.resize(capacity, 0)
        }
    }
}

/// Dense graph represented as a flat bit-matrix with MSB-0 bit ordering.
///
/// The adjacency matrix is stored row-major: row `v` occupies
/// `adjacency[v * words_per_vertex .. (v + 1) * words_per_vertex]`.
/// Bit `j` in row `v` (MSB-0) indicates an edge from `v` to `j`.
///
/// Construct with [`DenseGraph::new`] and add edges with [`DenseGraph::add_edge`]:
///
/// ```rust
/// use canonaut::structs::DenseGraph;
///
/// let mut graph = DenseGraph::new(3);
/// graph.add_edge(0, 1);
/// graph.add_edge(1, 2);
/// ```
///
/// For digraphs, construct with [`DenseGraph::new_directed`] instead and add arcs
/// with [`DenseGraph::add_arc`].  A graph built one way rejects the other kind of
/// edge at the call site (a single `bool` check, not a hot-path cost) rather than
/// silently producing a matrix that disagrees with how it's canonicalized —
/// [`CanonautManager::canonize_graph`](crate::structs::CanonautManager::canonize_graph)
/// reads this flag and configures the manager's digraph mode to match automatically.
///
/// A vertex coloring can be attached with [`set_colors`](Self::set_colors); when
/// present, [`canonize_graph`](crate::structs::CanonautManager::canonize_graph) uses
/// it automatically.  To canonicalize the same graph under a different coloring
/// per call (rather than one fixed coloring per graph), pass colors explicitly to
/// [`canonize_graph_with_colors`](crate::structs::CanonautManager::canonize_graph_with_colors)
/// instead — it ignores whatever is stored on the graph.
#[derive(Clone, Debug, PartialEq)]
#[repr(C)]
pub struct DenseGraph {
    /// Flat bit-matrix: `n * m` words, where `m = words_per_vertex`.
    pub adjacency: Vec<usize>,
    /// Number of vertices `n`.
    pub number_of_vertices: u32,
    /// Words per adjacency row: `m = ⌈n / 64⌉` on 64-bit targets.
    pub words_per_vertex: u32,
    /// Whether this graph was constructed via [`new_directed`](Self::new_directed).
    /// Not `pub`: set at construction, read by [`is_directed`](Self::is_directed).
    pub(crate) directed: bool,
    /// Optional vertex coloring, one color per vertex. Not `pub`: length-validated
    /// on write by [`set_colors`](Self::set_colors), read by [`colors`](Self::colors).
    pub(crate) colors: Option<Vec<u32>>,
}

impl DenseGraph {
    /// Creates an edgeless undirected graph on `number_of_vertices` vertices.
    /// Edges are added with [`add_edge`](Self::add_edge).
    pub fn new(number_of_vertices: usize) -> Self {
        Self::new_with_directedness(number_of_vertices, false)
    }

    /// Creates an edgeless directed graph (digraph) on `number_of_vertices` vertices.
    /// Arcs are added with [`add_arc`](Self::add_arc).
    pub fn new_directed(number_of_vertices: usize) -> Self {
        Self::new_with_directedness(number_of_vertices, true)
    }

    fn new_with_directedness(number_of_vertices: usize, directed: bool) -> Self {
        let words_per_vertex = words_needed(number_of_vertices) as u32;
        DenseGraph {
            adjacency: vec![0usize; number_of_vertices * words_per_vertex as usize],
            number_of_vertices: number_of_vertices as u32,
            words_per_vertex,
            directed,
            colors: None,
        }
    }

    /// True if this graph was constructed via [`new_directed`](Self::new_directed).
    pub fn is_directed(&self) -> bool {
        self.directed
    }

    pub(crate) fn has_loops(&self) -> bool {
        (0..self.number_of_vertices()).any(|v| {
            self.adjacency[v * self.words_per_vertex as usize + v / 64]
                & (1usize << (63 - v % 64)) != 0
        })
    }

    /// Number of vertices.
    pub fn number_of_vertices(&self) -> usize {
        self.number_of_vertices as usize
    }

    /// True if the edge `vertex1 → vertex2` is present.
    ///
    /// For an undirected graph the adjacency matrix is symmetric, so argument
    /// order does not matter; for a digraph this asks specifically about the
    /// arc from `vertex1` to `vertex2`.
    ///
    /// # Panics
    ///
    /// Panics if either vertex is out of range.
    pub fn has_edge(&self, vertex1: usize, vertex2: usize) -> bool {
        assert!(
            vertex1 < self.number_of_vertices as usize && vertex2 < self.number_of_vertices as usize,
            "vertex out of range: has_edge({vertex1}, {vertex2}) on a graph with {} vertices",
            self.number_of_vertices
        );
        get_graph_row(&self.adjacency, vertex1, self.words_per_vertex).is_bit_set(vertex2)
    }

    /// Adds an undirected edge between `vertex1` and `vertex2`.
    ///
    /// # Panics
    ///
    /// Panics if this graph was constructed via [`new_directed`](Self::new_directed) —
    /// use [`add_arc`](Self::add_arc) instead.
    pub fn add_edge(&mut self, vertex1: usize, vertex2: usize) {
        assert!(
            !self.directed,
            "add_edge() called on a directed DenseGraph; use add_arc() instead"
        );
        assert!(vertex1 < self.number_of_vertices() && vertex2 < self.number_of_vertices(),
            "vertex out of range");
        add_one_edge(&mut self.adjacency, vertex1, vertex2, self.words_per_vertex);
    }

    /// Adds the directed edge (arc) `source_vertex → target_vertex`.
    ///
    /// # Panics
    ///
    /// Panics unless this graph was constructed via [`new_directed`](Self::new_directed) —
    /// use [`add_edge`](Self::add_edge) for an undirected graph.
    pub fn add_arc(&mut self, source_vertex: usize, target_vertex: usize) {
        assert!(
            self.directed,
            "add_arc() called on an undirected DenseGraph; use add_edge(), or construct with DenseGraph::new_directed()"
        );
        assert!(source_vertex < self.number_of_vertices() && target_vertex < self.number_of_vertices(),
            "vertex out of range");
        add_one_arc(&mut self.adjacency, source_vertex, target_vertex, self.words_per_vertex);
    }

    /// Attaches a vertex coloring: `colors[v]` is an arbitrary `u32` color for
    /// vertex `v`. Colors need not be contiguous or sorted; values are fixed
    /// labels across graphs. Call at construction or any time before
    /// [`canonize_graph`](crate::structs::CanonautManager::canonize_graph);
    /// replaces any coloring set previously.
    ///
    /// # Panics
    ///
    /// Panics if `colors.len()` does not equal `number_of_vertices`.
    pub fn set_colors(&mut self, colors: Vec<u32>) {
        assert_eq!(
            colors.len(),
            self.number_of_vertices as usize,
            "colors.len() ({}) must equal number_of_vertices ({})",
            colors.len(),
            self.number_of_vertices,
        );
        self.colors = Some(colors);
    }

    /// The vertex coloring attached via [`set_colors`](Self::set_colors), if any.
    pub fn colors(&self) -> Option<&[u32]> {
        self.colors.as_deref()
    }
}

/// Returns the adjacency row for vertex `row_number` as a slice of `words_per_vertex` words.
///
/// Use this when `words_per_vertex > 1` (i.e., `n > 64`).
pub fn get_graph_row(graph: &[usize], row_number: usize, words_per_vertex: u32) -> &[usize] {
    let stride = words_per_vertex as usize;
    let start = row_number * stride;
    let end = start + stride;
    &graph[start..end]
}

/// Returns the adjacency row for vertex `row_number` when `words_per_vertex == 1` (`n ≤ 64`).
pub fn get_graph_row_single(graph: &[usize], row_number: usize) -> &[usize] {
    let start = row_number;
    &graph[start..start + 1]
}

/// Returns a mutable adjacency row for vertex `row_number` when `words_per_vertex > 1`.
pub fn get_graph_row_mut(
    graph: &mut [usize],
    row_number: usize,
    words_per_vertex: u32,
) -> &mut [usize] {
    let stride = words_per_vertex as usize;
    let start = row_number * stride;
    let end = start + stride;
    &mut graph[start..end]
}

/// Returns a mutable adjacency row for vertex `row_number` when `words_per_vertex == 1`.
pub fn get_graph_row_single_mut(graph: &mut [usize], row_number: usize) -> &mut [usize] {
    let start = row_number;
    &mut graph[start..start + 1]
}

/// Adds an undirected edge between vertices `vertex1` and `vertex2` in-place.
///
/// Sets bit `vertex2` in row `vertex1` and bit `vertex1` in row `vertex2`, so the
/// adjacency matrix remains symmetric.  `words_per_vertex` is `= ⌈n/64⌉`.
pub fn add_one_edge(graph: &mut [usize], vertex1: usize, vertex2: usize, words_per_vertex: u32) {
    let row_of_vertex1 = get_graph_row_mut(graph, vertex1, words_per_vertex);
    row_of_vertex1.set_bit(vertex2);

    let row_of_vertex2 = get_graph_row_mut(graph, vertex2, words_per_vertex);
    row_of_vertex2.set_bit(vertex1);
}

/// Adds the directed edge (arc) `source_vertex → target_vertex` to a dense bit-matrix graph.
///
/// Sets bit `target_vertex` in row `source_vertex` only; use together with the `digraph` option.
/// `words_per_vertex` is `= ⌈n/64⌉`.
pub fn add_one_arc(graph: &mut [usize], source_vertex: usize, target_vertex: usize, words_per_vertex: u32) {
    let source_row = get_graph_row_mut(graph, source_vertex, words_per_vertex);
    source_row.set_bit(target_vertex);
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Graph Row Access Tests ---

    #[test]
    fn test_get_graph_row_stride_1() {
        // A graph with 3 vertices, 1 word per vertex.
        // Graph: [Row0, Row1, Row2]
        let graph = vec![0xAAA, 0xBBB, 0xCCC];
        let words_per_vertex = 1;

        // Check Row 0
        let row0 = get_graph_row(&graph, 0, words_per_vertex);
        assert_eq!(row0, &[0xAAA]);

        // Check Row 2
        let row2 = get_graph_row(&graph, 2, words_per_vertex);
        assert_eq!(row2, &[0xCCC]);
    }

    #[test]
    fn test_get_graph_row_stride_2() {
        // A graph with 2 vertices, 2 words per vertex (words_per_vertex=2).
        // Row 0: [Word0a, Word0b]
        // Row 1: [Word1a, Word1b]
        let graph = vec![10, 11, 20, 21];
        let words_per_vertex = 2;

        // Check Row 0
        let row0 = get_graph_row(&graph, 0, words_per_vertex);
        assert_eq!(row0, &[10, 11]);

        // Check Row 1
        let row1 = get_graph_row(&graph, 1, words_per_vertex);
        assert_eq!(row1, &[20, 21]);
    }

    // --- Bit Manipulation Tests (MSB-0 Logic) ---

    #[test]
    fn has_edge_reads_back_what_was_added() {
        let mut graph = DenseGraph::new(70); // two words per vertex
        graph.add_edge(0, 69);
        graph.add_edge(3, 64);

        assert!(graph.has_edge(0, 69));
        assert!(graph.has_edge(69, 0), "undirected adjacency is symmetric");
        assert!(graph.has_edge(3, 64));
        assert!(!graph.has_edge(0, 1));
        assert_eq!(graph.number_of_vertices(), 70);
    }

    #[test]
    fn has_edge_distinguishes_arc_direction() {
        let mut graph = DenseGraph::new_directed(4);
        graph.add_arc(0, 1);

        assert!(graph.has_edge(0, 1));
        assert!(!graph.has_edge(1, 0));
    }

    #[test]
    #[should_panic(expected = "vertex out of range")]
    fn has_edge_rejects_out_of_range_vertices() {
        DenseGraph::new(4).has_edge(0, 4);
    }

    #[test]
    fn test_msb0_ordering() {
        let mut data = [0usize];

        // In Nauty/MSB-0: Index 0 is the HIGHEST bit.
        data.set_bit(0);
        // 1000...000 should be 1 << 63 (on 64-bit)
        assert_eq!(
            data[0],
            1 << (usize::BITS - 1),
            "Index 0 should be the High Bit (MSB)"
        );

        // Index 63 (on 64-bit) is the LOWEST bit.
        let last_bit = usize::BITS as usize - 1;
        data.set_bit(last_bit);
        // Should now have MSB and LSB set.
        let expected = (1 << (usize::BITS - 1)) | 1;
        assert_eq!(data[0], expected, "Last index should be the Low Bit (LSB)");
    }

    #[test]
    fn test_contains_and_unset() {
        let mut data = [0usize];

        assert!(!data.is_bit_set(5));

        data.set_bit(5);
        assert!(data.is_bit_set(5));

        data.unset_bit(5);
        assert!(!data.is_bit_set(5));
        assert_eq!(data[0], 0);
    }

    #[test]
    fn test_clear_words() {
        let mut data = [usize::MAX, usize::MAX, usize::MAX];

        // Clear first 2 words (simulating m=2)
        data.clear_words(2);

        assert_eq!(data[0], 0);
        assert_eq!(data[1], 0);
        assert_eq!(data[2], usize::MAX); // Should remain untouched
    }

    #[test]
    fn test_next_set_bit() {
        // Setup: [0000...0000, 0010...0001]
        // We put no bits in the first word, and two bits in the second word.
        let mut data = [0usize, 0usize];

        let bit_a = usize::BITS as usize + 2; // 3rd bit of second word
        let bit_b = usize::BITS as usize + (usize::BITS as usize - 1); // Last bit of second word

        data.set_bit(bit_a);
        data.set_bit(bit_b);

        // 1. Scan from beginning (should skip empty word 0)
        assert_eq!(data.next_set_bit(0), Some(bit_a));

        // 2. Scan from exactly the first set bit
        assert_eq!(data.next_set_bit(bit_a), Some(bit_a));

        // 3. Scan from just after the first set bit
        assert_eq!(data.next_set_bit(bit_a + 1), Some(bit_b));

        // 4. Scan from after the last bit
        assert_eq!(data.next_set_bit(bit_b + 1), None);
    }

    #[test]
    fn test_next_set_bit_complex_masking() {
        // Test finding a bit in the SAME word, but after the start index
        // to ensure we are masking the "left" bits correctly.
        let mut data = [0usize];

        // Set bit 0 (MSB) and bit 10.
        data.set_bit(0);
        data.set_bit(10);

        // If we start searching at 1, we should skip bit 0 and find bit 10.
        assert_eq!(data.next_set_bit(1), Some(10));
    }

    #[test]
    fn test_intersect_inplace() {
        // Word 0: 1111 (binary for simplification)
        // Word 1: 1010
        let mut target = [0b1111, 0b1010];
        let source = vec![0b0110, 0b1100];

        target.intersect_inplace(&source);

        // Word 0: 1111 & 0110 = 0110
        assert_eq!(target[0], 0b0110);
        // Word 1: 1010 & 1100 = 1000
        assert_eq!(target[1], 0b1000);
    }

    #[test]
    fn test_intersect_inplace_size_mismatch() {
        let mut target = [usize::MAX, usize::MAX];
        let source = vec![0]; // source is smaller

        target.intersect_inplace(&source);

        // First word intersected
        assert_eq!(target[0], 0);
        // Second word cleared (because source didn't have it, treated as 0)
        assert_eq!(target[1], 0);
    }
}
