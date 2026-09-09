//! Run with --features c-nauty-bench. Compares more than canonical adjacency:
//! orbit partitions and group orders catch lost generators during manager reuse.
#![cfg(feature = "c-nauty-bench")]
use canonaut::structs::CanonautManager;

// The reference library may be built without thread-local workspaces.
static C_NAUTY_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

unsafe extern "C" {
    fn c_nauty_reference(
        g: *const usize,
        lab: *mut i32,
        ptn: *mut i32,
        orbits: *mut i32,
        canon: *mut usize,
        m: i32,
        n: i32,
        directed: i32,
        colored: i32,
        size: *mut f64,
        exponent: *mut i32,
    );
}

fn compare(
    graph: &[usize],
    n: usize,
    directed: bool,
    colors: Option<&[u32]>,
    manager: &mut CanonautManager,
) {
    let _guard = C_NAUTY_LOCK.lock().unwrap();
    let words = n.div_ceil(64);
    let mut lab: Vec<i32> = (0..n as i32).collect();
    if let Some(colors) = colors {
        lab.sort_unstable_by_key(|&v| (colors[v as usize], v));
    }
    let mut partition: Vec<i32> = (0..n)
        .map(|i| {
            i32::from(
                i + 1 != n && colors.is_none_or(|c| c[lab[i] as usize] == c[lab[i + 1] as usize]),
            )
        })
        .collect();
    let mut orbits = vec![0; n];
    let mut canonical = vec![0; graph.len()];
    let (mut size, mut exponent) = (0.0, 0);
    // SAFETY: adjacency, permutation, partition and output buffers have the
    // dimensions required by densenauty; all labels and color cells are valid.
    unsafe {
        c_nauty_reference(
            graph.as_ptr(),
            lab.as_mut_ptr(),
            partition.as_mut_ptr(),
            orbits.as_mut_ptr(),
            canonical.as_mut_ptr(),
            words as i32,
            n as i32,
            directed as i32,
            colors.is_some() as i32,
            &mut size,
            &mut exponent,
        );
    }
    if let Some(colors) = colors {
        manager.canonize_with_colors(graph, colors);
    } else {
        manager.canonize(graph);
    }
    assert_eq!(
        manager.canonical_graph, canonical,
        "canonical form: n={n}, directed={directed}"
    );
    assert!(
        manager
            .orbits()
            .iter()
            .zip(&orbits)
            .all(|(&a, &b)| a == b as u32),
        "orbits: n={n}, directed={directed}"
    );
    let stats = manager.statistics();
    assert!(
        (stats.group_size1 / size - 1.0).abs() < 1e-9,
        "group size: n={n}, directed={directed}"
    );
    assert_eq!(stats.group_size2, exponent as u32);
}

#[test]
fn exhaustive_small_graphs_and_mixed_size_reuse_match_c_nauty() {
    for schreier in [false, true] {
        for directed in [false, true] {
            let mut manager = CanonautManager::new(4).with_canonization();
            if schreier {
                manager = manager.with_schreier();
            }
            if directed {
                manager = manager.with_digraph();
            }
            let slots: Vec<_> = (0..4)
                .flat_map(|i| {
                    (0..4)
                        .filter(move |&j| directed || j > i)
                        .map(move |j| (i, j))
                })
                .collect();
            for mask in 0..1usize << slots.len() {
                let mut graph = vec![0; 4];
                for (bit, &(u, v)) in slots.iter().enumerate() {
                    if mask & (1 << bit) != 0 {
                        graph[u] |= 1 << (63 - v);
                        if !directed {
                            graph[v] |= 1 << (63 - u);
                        }
                    }
                }
                compare(&graph, 4, directed, None, &mut manager);
            }
            let mut seed = 42u64;
            for trial in 0..1500 {
                let n: usize = [1, 2, 5, 9, 16, 33, 63, 64, 65, 70, 127, 128, 129, 193][trial % 14];
                let words = n.div_ceil(64);
                let mut graph = vec![0; n * words];
                for u in 0..n {
                    for v in 0..n {
                        seed ^= seed << 13;
                        seed ^= seed >> 7;
                        seed ^= seed << 17;
                        if (directed || v > u) && seed % 100 < [2, 10, 50, 90][trial % 4] {
                            graph[u * words + v / 64] |= 1 << (63 - v % 64);
                            if !directed {
                                graph[v * words + u / 64] |= 1 << (63 - u % 64);
                            }
                        }
                    }
                }
                let colors: Vec<u32> = (0..n).map(|v| ((v * 7 + trial) % 3) as u32).collect();
                compare(
                    &graph,
                    n,
                    directed,
                    if trial % 2 == 0 { Some(&colors) } else { None },
                    &mut manager,
                );
            }
        }
    }
}

#[test]
fn symmetric_families_with_schreier_match_c_nauty() {
    use canonaut::generators::*;
    let mut manager = CanonautManager::new(1).with_canonization().with_schreier();
    for n in [8, 64, 65, 127, 193, 9] {
        for graph in [
            generate_complete_graph(n),
            generate_cycle_graph(n),
            generate_star_graph(n),
        ] {
            compare(&graph, n, false, None, &mut manager);
            let colors: Vec<_> = (0..n).map(|v| (v % 3) as u32).collect();
            compare(&graph, n, false, Some(&colors), &mut manager);
        }
    }
    for n in [4, 10, 20] {
        compare(
            &generate_cfi_graph(n, false),
            n * 10,
            false,
            None,
            &mut manager,
        );
    }
    for n in [8, 16, 32] {
        compare(
            &generate_hadamard_graph(n),
            n * 2,
            false,
            None,
            &mut manager,
        );
    }
    for n in [13, 29, 61] {
        compare(&generate_paley_graph(n), n, false, None, &mut manager);
    }
}
