//! A `CanonautManager` is meant to be reused: that is the whole reason it owns
//! its buffers. Reuse has to survive a change in the *shape* of the graph, not
//! just its size -- in particular a drop in words-per-vertex, which happens the
//! moment a graph of more than 64 vertices is followed by a smaller one.

use canonaut::structs::CanonautManager;

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
