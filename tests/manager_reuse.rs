//! A `CanonautManager` is meant to be reused: that is the whole reason it owns
//! its buffers. Reuse has to survive a change in the *shape* of the graph, not
//! just its size -- in particular a drop in words-per-vertex, which happens the
//! moment a graph of more than 64 vertices is followed by a smaller one.

use canonaut::structs::CanonautManager;

#[test]
fn shrinking_active_cells_preserves_the_automorphism_group() {
    // The asymmetric graph can leave active bits above vertex 127. A reused
    // manager used to scan those stale words, reaching an unchecked read past
    // the smaller partition. In release builds it silently lost generators.
    let large = canonaut::generators::generate_random_graph(193, 0.1, 42);
    let n = 127;
    let small = vec![0; n * 2];
    let colors: Vec<u32> = (0..n).map(|v| if v < 3 { 0 } else { v as u32 }).collect();
    let mut manager = CanonautManager::new(193).with_canonization();
    for _ in 0..3 {
        manager.canonize(&large);
        manager.canonize_with_colors(&small, &colors);
        // Only the first three isolated vertices can be permuted: S_3.
        assert_eq!(manager.statistics().group_size1, 6.0);
        assert_eq!(manager.statistics().number_of_orbits, 125);
        assert_eq!(&manager.orbits()[..3], &[0, 0, 0]);
        assert_eq!(manager.canonical_graph, small);
    }
}

#[test]
fn empty_graph_resets_results_between_nonempty_runs() {
    use canonaut::structs::DenseGraph;
    let mut manager = CanonautManager::new(65).with_canonization();
    manager.canonize(&canonaut::generators::generate_cycle_graph(65));
    for colored in [false, true] {
        let mut empty = DenseGraph::new(0);
        if colored { empty.set_colors(vec![]); }
        manager.canonize_graph(&empty);
        assert!(manager.canonical_graph.is_empty());
        assert!(manager.labeling().is_empty());
        assert!(manager.orbits().is_empty());
        assert_eq!(manager.statistics().group_size1, 1.0);
        assert_eq!(manager.statistics().group_size2, 0);
        assert_eq!(manager.statistics().number_of_orbits, 0);
        assert_eq!(manager.statistics().number_of_generators, 0);
        manager.canonize_with_colors(&[0; 3], &[0; 3]);
        assert_eq!(manager.statistics().group_size1, 6.0);
    }
}

#[test]
fn raw_and_dense_entry_points_agree_across_word_boundaries() {
    use canonaut::structs::DenseGraph;
    let mut raw = CanonautManager::new(0).with_canonization();
    let mut dense = CanonautManager::new(0).with_canonization();
    for n in [1, 2, 63, 64, 65, 127, 128, 129, 257, 64, 0] {
        let mut graph = DenseGraph::new(n);
        for v in 1..n { graph.add_edge(v - 1, v); }
        raw.canonize(&graph.adjacency);
        dense.canonize_graph(&graph);
        assert_eq!(raw.canonical_graph, dense.canonical_graph, "n={n}");
        assert_eq!(raw.orbits(), dense.orbits(), "n={n}");
        assert_eq!(raw.statistics().group_size1, dense.statistics().group_size1);
    }
}

#[test]
#[should_panic(expected = "valid (n, m)")]
fn raw_adjacency_rejects_invalid_dimensions() {
    CanonautManager::new(0).canonize(&[0; 65]);
}

/// A deterministic pseudo-random graph, as a dense bitset of `n * m` words.
fn random_graph(n: usize, seed: &mut u64) -> Vec<usize> {
    let mut next = || {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        *seed
    };
    let words = n.div_ceil(64);
    let mut graph = vec![0usize; n * words];
    for a in 0..n {
        for b in (a + 1)..n {
            if next() % 4 == 0 {
                graph[a * words + b / 64] |= 1usize << (64 - 1 - (b % 64));
                graph[b * words + a / 64] |= 1usize << (64 - 1 - (a % 64));
            }
        }
    }
    graph
}

fn colors(n: usize) -> Vec<u32> {
    (0..n as u32).map(|vertex| vertex % 3).collect()
}


/// 91 vertices, 2 word(s) per vertex.
const EDGES_91: &[(u16, u16)] = &[
    (0,64), (1,65), (2,62), (3,59), (4,63), (5,34), (5,50), (6,84), (7,85), (8,86),
    (9,87), (10,88), (11,89), (12,41), (12,52), (12,57), (13,45), (13,55), (13,56), (14,38),
    (14,39), (14,40), (15,43), (15,47), (15,48), (16,37), (16,49), (16,61), (17,36), (17,44),
    (17,46), (18,66), (18,67), (19,68), (19,69), (20,70), (20,71), (21,72), (21,73), (22,68),
    (22,74), (23,75), (23,76), (24,77), (24,78), (25,79), (25,80), (26,66), (26,81), (27,70),
    (27,82), (28,74), (28,82), (29,72), (29,79), (30,77), (30,83), (31,73), (31,75), (32,58),
    (32,76), (32,80), (33,60), (33,81), (33,83), (34,67), (34,78), (35,60), (35,69), (35,71),
    (40,50), (41,43), (42,57), (42,61), (44,49), (45,52), (46,54), (47,59), (48,64), (51,53),
    (51,54), (53,58), (55,65), (56,63), (57,62), (58,90), (59,87), (60,90), (61,89), (62,88),
    (63,84), (64,86), (65,85),
];
const COLORS_91: &[u32] = &[
    0, 0, 0, 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 9, 9,
    9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11,
];

/// 36 vertices, 1 word(s) per vertex.
const EDGES_36: &[(u16, u16)] = &[
    (0,17), (0,20), (1,34), (2,35), (3,18), (3,24), (3,25), (4,5), (4,22), (5,23),
    (6,19), (6,23), (7,24), (7,26), (8,27), (8,28), (9,27), (9,29), (10,29), (10,30),
    (11,30), (11,31), (12,22), (12,26), (12,32), (13,28), (13,32), (13,33), (14,25), (14,31),
    (14,33), (15,21), (16,21), (18,21), (19,20), (22,34), (23,35),
];
const COLORS_36: &[u32] = &[
    0, 1, 1, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7,
    7, 8, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11,
];

/// Rebuilds a dense bitset graph from an edge list.
fn from_edges(n: usize, edges: &[(u16, u16)]) -> Vec<usize> {
    let words = n.div_ceil(64);
    let mut graph = vec![0usize; n * words];
    for &(a, b) in edges {
        let (a, b) = (a as usize, b as usize);
        graph[a * words + b / 64] |= 1usize << (64 - 1 - (b % 64));
        graph[b * words + a / 64] |= 1usize << (64 - 1 - (a % 64));
    }
    graph
}

/// One manager, a graph needing two words per vertex, then one needing a single
/// word. The second run used to index the first run's vertices: the target-cell
/// pool handed back a node whose bitset still held bits from the wider graph.
///
/// These two are real molecular graphs, in the order that first hit it. Random
/// graphs do not: they refine to discrete cells too quickly to ever ask the
/// pool for a reused node at a level where the stale bits are still set.
#[test]
fn manager_survives_a_real_drop_in_words_per_vertex() {
    let mut manager = CanonautManager::new(64).with_canonization();
    for (n, edges, colors) in [
        (91usize, EDGES_91, COLORS_91),
        (36usize, EDGES_36, COLORS_36),
    ] {
        let graph = from_edges(n, edges);
        manager.canonize_with_colors(&graph, colors);
        for &vertex in &manager.labeling()[..n] {
            assert!(
                (vertex as usize) < n,
                "n={n}: labelling names vertex {vertex}, which does not exist here"
            );
        }
    }
}

#[test]
fn manager_survives_a_drop_in_words_per_vertex() {
    let mut seed = 0x5EED_1234u64;
    let mut manager = CanonautManager::new(64).with_canonization();

    for round in 0..40 {
        for &n in &[91usize, 36] {
            let graph = random_graph(n, &mut seed);
            manager.canonize_with_colors(&graph, &colors(n));
            for &vertex in &manager.labeling()[..n] {
                assert!(
                    (vertex as usize) < n,
                    "round {round}, n={n}: labelling names vertex {vertex}, \
                     which does not exist in this graph"
                );
            }
        }
    }
}

/// The same graphs, each on a manager of its own, must give the same answer as
/// the reused one -- otherwise reuse is not merely crashing but silently
/// changing the canonical form.
#[test]
fn reuse_gives_the_same_labelling_as_a_fresh_manager() {
    let mut seed = 0xC0FF_EE00u64;
    let mut reused = CanonautManager::new(16).with_canonization();

    for &n in &[80usize, 20, 70, 33, 65, 12] {
        let graph = random_graph(n, &mut seed);
        let colors = colors(n);

        let mut fresh = CanonautManager::new(n).with_canonization();
        fresh.canonize_with_colors(&graph, &colors);
        let expected = fresh.labeling()[..n].to_vec();

        reused.canonize_with_colors(&graph, &colors);
        assert_eq!(&reused.labeling()[..n], &expected[..], "n={n}");
    }
}

#[test]
fn dense_undirected_loops_use_full_automorphism_checks() {
    use canonaut::structs::DenseGraph;
    let mut graph = DenseGraph::new(3);
    graph.add_edge(0, 0);
    let mut manager = CanonautManager::new(3).with_canonization();
    manager.canonize_graph(&graph);
    assert_eq!(manager.statistics().group_size1, 2.0);
    assert_eq!(manager.orbits(), &[0, 1, 1]);
    assert!(!graph.is_directed());
    // Switching back to a loopless graph must reset the effective mode.
    manager.canonize_graph(&DenseGraph::new(3));
    assert_eq!(manager.statistics().group_size1, 6.0);
}

#[test]
fn prepared_buffers_and_public_output_mutation_preserve_reuse() {
    use canonaut::generators::generate_cycle_graph;
    let mut reused = CanonautManager::new(8).with_canonization();
    for n in [8usize, 64, 3, 65, 0, 8] {
        let graph = generate_cycle_graph(n);
        let mut fresh = CanonautManager::new(n).with_canonization();
        fresh.canonize(&graph);
        // Public preparation can change the private buffers' logical lengths.
        reused.prepare_for_run(n + 1, (n + 1).div_ceil(64));
        for output_size in [0, 1, graph.len() + 1] {
            reused.canonical_graph.resize(output_size, usize::MAX);
            reused.canonize(&graph);
            assert_eq!(reused.canonical_graph, fresh.canonical_graph);
            assert_eq!(reused.labeling(), fresh.labeling());
            assert_eq!(reused.orbits(), fresh.orbits());
            let taken = reused.take_canonical_graph();
            assert_eq!(taken, fresh.canonical_graph);
        }
    }
}

#[test]
#[should_panic(expected = "vertex out of range")]
fn dense_graph_rejects_arcs_into_padding_bits() {
    canonaut::structs::DenseGraph::new_directed(3).add_arc(0, 3);
}
