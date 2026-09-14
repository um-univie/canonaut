//! Repeat a workload for a sampling profiler or Callgrind.
//!
//! cargo run --profile flamegraph --bin profile -- random 512 1000
//! cargo run --profile flamegraph --bin profile -- colored 64 100000
//! cargo run --profile flamegraph --bin profile -- cfi 20 100
//!
//! For CFI, n is the base order (graph order 10n); for Hadamard it is
//! the matrix order (graph order 2n). All graph generation is outside the loop.
use canonaut::{
    generators::{
        generate_cfi_graph, generate_complete_graph, generate_cycle_graph, generate_hadamard_graph,
        generate_paley_graph, generate_path_graph, generate_random_graph, generate_star_graph,
    },
    structs::CanonautManager,
};
use std::hint::black_box;

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let family = args.get(1).map(String::as_str).unwrap_or("random");
    let n: usize = args
        .get(2)
        .map(|s| s.parse().expect("invalid order"))
        .unwrap_or(512);
    let iterations: usize = args
        .get(3)
        .map(|s| s.parse().expect("invalid iteration count"))
        .unwrap_or(1000);
    let (graph, vertices) = match family {
        "random" => (generate_random_graph(n, 0.5, 42), n),
        "random10" | "colored" => (generate_random_graph(n, 0.1, 42), n),
        "complete" => (generate_complete_graph(n), n),
        "cycle" => (generate_cycle_graph(n), n),
        "path" => (generate_path_graph(n), n),
        "star" => (generate_star_graph(n), n),
        "cfi" => (generate_cfi_graph(n, false), 10 * n),
        "hadamard" => (generate_hadamard_graph(n), 2 * n),
        "paley" => (generate_paley_graph(n), n),
        other => panic!(
            "unknown family {other}; use random|random10|colored|complete|cycle|path|star|cfi|hadamard|paley"
        ),
    };
    let mut manager = CanonautManager::new(vertices).with_canonization();
    if args.iter().any(|arg| arg == "--schreier") {
        manager = manager.with_schreier();
    }
    let start = std::time::Instant::now();
    let colors: Vec<u32> = (0..vertices).map(|v| (v % 3) as u32).collect();
    for _ in 0..iterations {
        if family == "colored" {
            manager.canonize_with_colors(black_box(&graph), black_box(&colors));
        } else {
            manager.canonize(black_box(&graph));
        }
        black_box(&manager.canonical_graph);
    }
    println!(
        "ns_per_graph={} nodes={} orbits={}",
        start.elapsed().as_nanos() / iterations.max(1) as u128,
        manager.statistics().number_of_nodes,
        manager.statistics().number_of_orbits
    );
}
