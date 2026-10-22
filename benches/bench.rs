use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main, SamplingMode};
use std::hint::black_box;
use canonaut::{
    generators::{
        generate_cfi_graph, generate_complete_graph, generate_cycle_graph,
        generate_hadamard_graph, generate_paley_graph, generate_path_graph,
        generate_random_graph, generate_star_graph,
    },
    graph6::parse_graph6_file,
    structs::{CanonautManager, ValidatedGraph},
    utilities::words_needed,
};

unsafe extern "C" {
    fn c_nauty_dense_canonize(
        graph: *mut u64,
        labeling: *mut i32,
        partition: *mut i32,
        orbits: *mut i32,
        canonical_graph: *mut u64,
        words_per_vertex: i32,
        vertex_count: i32,
    );
}

struct CNautyBuffers {
    labeling: Vec<i32>,
    partition: Vec<i32>,
    orbits: Vec<i32>,
    canonical_graph: Vec<u64>,
}

impl CNautyBuffers {
    fn new(vertex_count: usize, words_per_vertex: usize) -> Self {
        CNautyBuffers {
            labeling: vec![0; vertex_count],
            partition: vec![0; vertex_count],
            orbits: vec![0; vertex_count],
            canonical_graph: vec![0; vertex_count * words_per_vertex],
        }
    }
}

#[inline]
fn run_c_nauty(graph: &[usize], c_nauty_buffers: &mut CNautyBuffers, words_per_vertex: i32, vertex_count: i32) {
    unsafe {
        c_nauty_dense_canonize(
            graph.as_ptr() as *mut u64,
            c_nauty_buffers.labeling.as_mut_ptr(),
            c_nauty_buffers.partition.as_mut_ptr(),
            c_nauty_buffers.orbits.as_mut_ptr(),
            c_nauty_buffers.canonical_graph.as_mut_ptr(),
            words_per_vertex,
            vertex_count,
        );
    }
}

/// Sizes for O(n³)+ graphs (complete, star) — stop at 128 to avoid minute-long runs.
const SIZES_HARD: &[usize] = &[8, 16, 32, 64, 128];
/// Sizes for tractable graphs (cycle, path, random).
const SIZES_EASY: &[usize] = &[8, 16, 32, 64, 128, 256];

/// Keep validation-inclusive timings distinct from validation-once timings.
/// Existing canonaut-vs-C groups continue to use the checked raw-slice API.
fn bench_input_validation(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("input_validation");
    for vertex_count in [8, 36, 63, 64, 65, 129, 512] {
        let graph = generate_random_graph(vertex_count, 0.1, 42);
        let validated = ValidatedGraph::new(&graph);
        group.bench_function(BenchmarkId::new("checked", vertex_count), |bencher| {
            let mut manager = CanonautManager::new(vertex_count).with_canonization();
            bencher.iter(|| manager.canonize(black_box(&graph)));
        });
        group.bench_function(BenchmarkId::new("validated", vertex_count), |bencher| {
            let mut manager = CanonautManager::new(vertex_count).with_canonization();
            bencher.iter(|| manager.canonize_validated(black_box(&validated)));
        });
    }
    group.finish();
}

fn bench_complete(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("complete");
    for &vertex_count in SIZES_HARD {
        let words_per_vertex = words_needed(vertex_count);
        let graph = generate_complete_graph(vertex_count);

        group.bench_with_input(BenchmarkId::new("canonaut", vertex_count), &vertex_count, |bencher, _| {
            let mut manager = CanonautManager::new(vertex_count).with_canonization();
            bencher.iter(|| manager.canonize(black_box(&graph)));
        });

        group.bench_with_input(BenchmarkId::new("nauty", vertex_count), &vertex_count, |bencher, _| {
            let mut c_nauty_buffers = CNautyBuffers::new(vertex_count, words_per_vertex);
            bencher.iter(|| run_c_nauty(black_box(&graph), &mut c_nauty_buffers, words_per_vertex as i32, vertex_count as i32));
        });
    }
    group.finish();
}

fn bench_cycle(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("cycle");
    for &vertex_count in SIZES_EASY {
        let words_per_vertex = words_needed(vertex_count);
        let graph = generate_cycle_graph(vertex_count);

        group.bench_with_input(BenchmarkId::new("canonaut", vertex_count), &vertex_count, |bencher, _| {
            let mut manager = CanonautManager::new(vertex_count).with_canonization();
            bencher.iter(|| manager.canonize(black_box(&graph)));
        });

        group.bench_with_input(BenchmarkId::new("nauty", vertex_count), &vertex_count, |bencher, _| {
            let mut c_nauty_buffers = CNautyBuffers::new(vertex_count, words_per_vertex);
            bencher.iter(|| run_c_nauty(black_box(&graph), &mut c_nauty_buffers, words_per_vertex as i32, vertex_count as i32));
        });
    }
    group.finish();
}

fn bench_path(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("path");
    for &vertex_count in SIZES_EASY {
        let words_per_vertex = words_needed(vertex_count);
        let graph = generate_path_graph(vertex_count);

        group.bench_with_input(BenchmarkId::new("canonaut", vertex_count), &vertex_count, |bencher, _| {
            let mut manager = CanonautManager::new(vertex_count).with_canonization();
            bencher.iter(|| manager.canonize(black_box(&graph)));
        });

        group.bench_with_input(BenchmarkId::new("nauty", vertex_count), &vertex_count, |bencher, _| {
            let mut c_nauty_buffers = CNautyBuffers::new(vertex_count, words_per_vertex);
            bencher.iter(|| run_c_nauty(black_box(&graph), &mut c_nauty_buffers, words_per_vertex as i32, vertex_count as i32));
        });
    }
    group.finish();
}

fn bench_star(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("star");
    for &vertex_count in SIZES_HARD {
        let words_per_vertex = words_needed(vertex_count);
        let graph = generate_star_graph(vertex_count);

        group.bench_with_input(BenchmarkId::new("canonaut", vertex_count), &vertex_count, |bencher, _| {
            let mut manager = CanonautManager::new(vertex_count).with_canonization();
            bencher.iter(|| manager.canonize(black_box(&graph)));
        });

        group.bench_with_input(BenchmarkId::new("nauty", vertex_count), &vertex_count, |bencher, _| {
            let mut c_nauty_buffers = CNautyBuffers::new(vertex_count, words_per_vertex);
            bencher.iter(|| run_c_nauty(black_box(&graph), &mut c_nauty_buffers, words_per_vertex as i32, vertex_count as i32));
        });
    }
    group.finish();
}

fn bench_random(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("random_p50");
    for &vertex_count in SIZES_EASY {
        let words_per_vertex = words_needed(vertex_count);
        let graph = generate_random_graph(vertex_count, 0.5, 42);

        group.bench_with_input(BenchmarkId::new("canonaut", vertex_count), &vertex_count, |bencher, _| {
            let mut manager = CanonautManager::new(vertex_count).with_canonization();
            bencher.iter(|| manager.canonize(black_box(&graph)));
        });

        group.bench_with_input(BenchmarkId::new("nauty", vertex_count), &vertex_count, |bencher, _| {
            let mut c_nauty_buffers = CNautyBuffers::new(vertex_count, words_per_vertex);
            bencher.iter(|| run_c_nauty(black_box(&graph), &mut c_nauty_buffers, words_per_vertex as i32, vertex_count as i32));
        });
    }
    group.finish();
}

fn bench_g6_batch(criterion: &mut Criterion) {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));

    let files: &[(&str, &str)] = &[
        ("cubic_10", "test_data/cubic_10.g6"),
        ("all_9", "test_data/all_9.g6"),
        ("bipartite_10", "test_data/bipartite_10.g6"),
    ];

    let mut group = criterion.benchmark_group("g6_batch");
    group.sample_size(20);

    for (label, rel_path) in files {
        let data = std::fs::read(root.join(rel_path)).expect("read g6 file");
        let graphs: Vec<_> = parse_graph6_file(&data)
            .into_iter()
            .take(1000)
            .filter_map(Result::ok)
            .collect();

        let max_vertex_count = graphs.iter().map(|graph| graph.number_of_vertices as usize).max().unwrap_or(16);

        group.bench_function(format!("canonaut/{label}"), |bencher| {
            let mut manager = CanonautManager::new(max_vertex_count).with_canonization();
            bencher.iter(|| {
                for graph in black_box(&graphs) {
                    manager.canonize(&graph.adjacency);
                }
            });
        });

        group.bench_function(format!("nauty/{label}"), |bencher| {
            let max_words_per_vertex = words_needed(max_vertex_count);
            let mut c_nauty_buffers = CNautyBuffers::new(max_vertex_count, max_words_per_vertex);
            bencher.iter(|| {
                for graph in black_box(&graphs) {
                    let vertex_count = graph.number_of_vertices as usize;
                    let words_per_vertex = graph.words_per_vertex as usize;
                    if c_nauty_buffers.labeling.len() < vertex_count {
                        c_nauty_buffers.labeling.resize(vertex_count, 0);
                        c_nauty_buffers.partition.resize(vertex_count, 0);
                        c_nauty_buffers.orbits.resize(vertex_count, 0);
                        c_nauty_buffers.canonical_graph.resize(vertex_count * words_per_vertex, 0);
                    }
                    run_c_nauty(&graph.adjacency, &mut c_nauty_buffers, words_per_vertex as i32, vertex_count as i32);
                }
            });
        });
    }

    group.finish();
}

/// Large-n scaling benchmarks: 5 replicates, 1 s warm-up, flat sampling.
/// Covers sizes where the standard 100-sample run would take too long.
const SIZES_SCALE_EASY: &[usize] = &[256, 512, 1024, 2048];
const SIZES_SCALE_HARD: &[usize] = &[256, 512];

fn bench_scaling(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("scaling");
    group.sample_size(10);
    group.warm_up_time(std::time::Duration::from_secs(1));
    group.sampling_mode(SamplingMode::Flat);

    for &vertex_count in SIZES_SCALE_EASY {
        let words_per_vertex = words_needed(vertex_count);
        for (label, graph) in [
            ("cycle", generate_cycle_graph(vertex_count)),
            ("path",  generate_path_graph(vertex_count)),
            ("random_p50", generate_random_graph(vertex_count, 0.5, 42)),
        ] {
            group.bench_with_input(BenchmarkId::new(format!("canonaut/{label}"), vertex_count), &vertex_count, |bencher, _| {
                let mut manager = CanonautManager::new(vertex_count).with_canonization();
                bencher.iter(|| manager.canonize(black_box(&graph)));
            });
            group.bench_with_input(BenchmarkId::new(format!("nauty/{label}"), vertex_count), &vertex_count, |bencher, _| {
                let mut c_nauty_buffers = CNautyBuffers::new(vertex_count, words_per_vertex);
                bencher.iter(|| run_c_nauty(black_box(&graph), &mut c_nauty_buffers, words_per_vertex as i32, vertex_count as i32));
            });
        }
    }

    for &vertex_count in SIZES_SCALE_HARD {
        let words_per_vertex = words_needed(vertex_count);
        for (label, graph) in [
            ("complete", generate_complete_graph(vertex_count)),
            ("star",     generate_star_graph(vertex_count)),
        ] {
            group.bench_with_input(BenchmarkId::new(format!("canonaut/{label}"), vertex_count), &vertex_count, |bencher, _| {
                let mut manager = CanonautManager::new(vertex_count).with_canonization();
                bencher.iter(|| manager.canonize(black_box(&graph)));
            });
            group.bench_with_input(BenchmarkId::new(format!("nauty/{label}"), vertex_count), &vertex_count, |bencher, _| {
                let mut c_nauty_buffers = CNautyBuffers::new(vertex_count, words_per_vertex);
                bencher.iter(|| run_c_nauty(black_box(&graph), &mut c_nauty_buffers, words_per_vertex as i32, vertex_count as i32));
            });
        }
    }

    group.finish();
}

/// Refinement-blind ("hard") instances: colour refinement distinguishes
/// nothing on these, so runtime is dominated by the search tree instead of the
/// refinement inner loops.  Constructed in-process, so both implementations
/// canonicalize the identical buffer.
///   cfi/<n>      Cai-Furer-Immerman graph over the 3-regular circulant,
///                n = 10 x base vertex count
///   hadamard/<n> Sylvester Hadamard bipartite graph, n = 2 x matrix order
///   paley/<n>    Paley graph on n vertices (strongly regular)
fn bench_hard(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("hard");
    group.sample_size(10);
    group.warm_up_time(std::time::Duration::from_secs(1));
    group.sampling_mode(SamplingMode::Flat);

    let mut instances: Vec<(&str, usize, Vec<usize>)> = Vec::new();
    for base_vertex_count in [20usize, 40, 60, 100] {
        instances.push(("cfi", 10 * base_vertex_count, generate_cfi_graph(base_vertex_count, false)));
    }
    for order in [64usize, 128, 256] {
        instances.push(("hadamard", 2 * order, generate_hadamard_graph(order)));
    }
    for prime_order in [101usize, 257, 509] {
        instances.push(("paley", prime_order, generate_paley_graph(prime_order)));
    }

    for (label, vertex_count, graph) in instances {
        let words_per_vertex = words_needed(vertex_count);
        group.bench_with_input(BenchmarkId::new(format!("canonaut/{label}"), vertex_count), &vertex_count, |bencher, _| {
            let mut manager = CanonautManager::new(vertex_count).with_canonization();
            bencher.iter(|| manager.canonize(black_box(&graph)));
        });
        group.bench_with_input(BenchmarkId::new(format!("nauty/{label}"), vertex_count), &vertex_count, |bencher, _| {
            let mut c_nauty_buffers = CNautyBuffers::new(vertex_count, words_per_vertex);
            bencher.iter(|| run_c_nauty(black_box(&graph), &mut c_nauty_buffers, words_per_vertex as i32, vertex_count as i32));
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_complete,
    bench_cycle,
    bench_path,
    bench_star,
    bench_random,
    bench_g6_batch,
    bench_scaling,
    bench_hard,
    bench_input_validation,
);
criterion_main!(benches);
