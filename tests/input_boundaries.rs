use canonaut::{
    generators::generate_cycle_graph,
    structs::{CanonautManager, DenseGraph, ValidatedGraph},
};
use std::panic::{AssertUnwindSafe, catch_unwind};

#[test]
fn every_adjacency_entry_point_rejects_padding_before_search() {
    for n in [1usize, 2, 63, 65, 127, 129, 193] {
        for row in [0, n - 1] {
            let mut graph = DenseGraph::new(n);
            let words = n.div_ceil(64);
            graph.adjacency[row * words + words - 1] = 1;
            let colors: Vec<_> = (0..n as u32).collect();
            let mut manager = CanonautManager::new(n).with_canonization();
            // Even an already-discrete colored partition must reject padding:
            // otherwise it reaches the unchecked canonical permutation directly.
            for entry in 0..5 {
                assert!(catch_unwind(AssertUnwindSafe(|| match entry {
                    0 => manager.canonize(&graph.adjacency),
                    1 => manager.canonize_with_colors(&graph.adjacency, &colors),
                    2 => manager.canonize_graph(&graph),
                    3 => manager.canonize_graph_with_colors(&graph, &colors),
                    _ => { ValidatedGraph::new(&graph.adjacency); }
                })).is_err(), "n={n}, row={row}, entry={entry}");
            }
            // Rejected inputs leave the manager usable.
            graph.adjacency.fill(0);
            manager.canonize_graph_with_colors(&graph, &colors);
            assert_eq!(manager.statistics().group_size1, 1.0);
        }
    }
}

#[test]
fn validated_views_match_raw_search_across_word_boundaries() {
    let mut raw = CanonautManager::new(0).with_canonization();
    let mut checked = CanonautManager::new(0).with_canonization();
    for n in [0usize, 1, 2, 63, 64, 65, 127, 128, 129, 193, 3] {
        let graph = generate_cycle_graph(n);
        let view = ValidatedGraph::new(&graph);
        let colors: Vec<_> = (0..n).map(|v| (v % 3) as u32).collect();
        for colored in [false, true] {
            if colored {
                raw.canonize_with_colors(&graph, &colors);
                checked.canonize_validated_with_colors(&view, &colors);
            } else {
                raw.canonize(&graph);
                checked.canonize_validated(&view);
            }
            assert_eq!(raw.canonical_graph, checked.canonical_graph);
            assert_eq!(raw.orbits(), checked.orbits());
            assert_eq!(raw.statistics().group_size1, checked.statistics().group_size1);
        }
    }
}

#[test]
fn immutable_callback_observes_current_order_after_manager_shrinks() {
    fn observe(_: u32, permutation: &[u32], orbits: &[u32], _: u32, _: u32, n: u32) {
        assert_eq!(permutation.len(), n as usize);
        assert_eq!(orbits.len(), n as usize);
        let mut sorted = permutation.to_vec();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..n).collect::<Vec<_>>());
        assert!(orbits.iter().all(|&v| v < n));
    }
    for schreier in [false, true] {
        let mut manager = CanonautManager::new(193).with_automorphism_callback(observe);
        if schreier { manager = manager.with_schreier(); }
        for n in [193, 65, 4] {
            manager.canonize(&generate_cycle_graph(n));
            assert_eq!(manager.statistics().group_size1, (2 * n) as f64);
            assert!(manager.statistics().number_of_generators > 0);
        }
    }
}

#[test]
fn invalid_lengths_and_dense_metadata_are_rejected() {
    assert!(catch_unwind(|| ValidatedGraph::new(&[0; 65])).is_err());
    let mut graph = DenseGraph::new(3);
    graph.words_per_vertex = 0;
    assert!(catch_unwind(|| CanonautManager::new(3).canonize_graph(&graph)).is_err());
}
