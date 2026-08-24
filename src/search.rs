use crate::{
    consts::ALGORITHM_INFINITY,
    io::OutputStream,
    refine::GraphContext,
    support::{
        labelorg, kill_requested, individualize_vertex, apply_refinement, fmperm_safe, mark_partition_fixed_points, prune_by_generators, select_target_cell,
        check_build_consistency,
    },
    simd_ops,
    structs::{BitSetOpsMSB0, manager::*},
    utilities::safe_orbjoin,
};

// Approximation of real group size
pub fn multiply(number1: &mut f64, number2: &mut u32, factor: u32) {
    *number1 *= factor as f64;
    if *number1 > 1e10 {
        *number1 /= 1e10;
        *number2 += 10;
    }
}


/// Checks `words_per_vertex`/`num_vertices` against the algorithm's internal limits,
/// returning an error message if either is out of range (mirrors nauty's own limits).
///
/// # Examples
///
/// ```
/// use canonaut::search::validate_input;
/// use canonaut::consts::ALGORITHM_INFINITY;
/// let result = validate_input(1, 10);
/// assert!(result.is_none());
///
/// let result = validate_input(ALGORITHM_INFINITY / 64 + 2, 10);
/// assert!(result.is_some());
/// ```
pub fn validate_input(words_per_vertex: u32, num_vertices: u32) -> Option<String> {
    if words_per_vertex > ALGORITHM_INFINITY / 64 + 1 {
        return Some(format!(
            "canonaut: need words_per_vertex <= {}, but words_per_vertex={}",
            ALGORITHM_INFINITY / 64 + 1,
            words_per_vertex
        ));
    }

    if num_vertices > ALGORITHM_INFINITY - 2 || num_vertices > 64 * words_per_vertex {
        return Some(format!(
            "canonaut: need num_vertices <= min({},{}*words_per_vertex), but num_vertices={}",
            ALGORITHM_INFINITY - 2,
            64,
            num_vertices
        ));
    }

    None
}

pub fn init_empty_stats(stats_block: &mut Statistics) {
    stats_block.group_size1 = 1.0f64;
    stats_block.group_size2 = 0;
    stats_block.number_of_orbits = 0;
    stats_block.number_of_generators = 0;
    stats_block.error_status_code = 0;
    stats_block.number_of_nodes = 1;
    stats_block.number_of_bad_leaves = 0;
    stats_block.maximum_search_depth = 1;
    stats_block.total_target_cell_members = 0;
    stats_block.canonical_updates = 0;
    stats_block.invariant_procedure_calls = 0;
    stats_block.successful_invariant_calls = 0;
    stats_block.shallowest_invariant_success_level = 0;
}

/// True if the dispatch vector has every function pointer `canonicalize()` requires.
pub fn validate_dispatch_vector(dispatch_vec: &DispatchVec) -> bool {
    dispatch_vec.refine.is_some()
        && dispatch_vec.update_canonical_graph.is_some()
        && dispatch_vec.target_cell.is_some()
        && dispatch_vec.cheap_automorphism_check.is_some()
}


/// Finds generators for the automorphism group of a vertex-colored graph and
/// optionally produces a canonically labeled isomorph.
///
/// # Arguments
///
/// * `graph` - The input graph represented as a flattened adjacency list or matrix.
/// * `labeling` - The partition nest defining the coloring of the graph. On return,
///   if `options.getcanon` is true, this contains the canonical labeling.
/// * `partition` - Array defining the initial partition of vertices. If `options.defaultpartition`
///   is true, the program sets this. Otherwise, it must be provided.
/// * `active` - An optional initial set of active colors. Used if `options.defaultpartition` is false.
/// * `orbits` - On return, `orbits[i]` contains the index of the lowest-numbered vertex
///   in the same orbit as vertex `i`.
/// * `options` - Configuration options for the algorithm (e.g., getting canonical labels, digraphs).
/// * `stats` - Output structure for algorithm statistics (e.g., group size, number of nodes).
/// * `workspace` - A mutable chunk of memory used for working storage.
/// * `words_per_vertex` - The number of `usize` words required to represent sets for a single vertex.
/// * `number_of_vertices` - The total number of vertices (`n`) in the graph.
/// * `canonical_graph` - The canonically labeled isomorph of the graph. This is only populated
///   if `options.getcanon` is true.
///
/// # Panics
///
/// * Panics if the dispatch vector (`options.dispatch`) is improperly configured.
/// * Panics if the internal capacity checks for words and vertices exceed physical limits.
///
/// # Side Effects
///
/// * Modifies `stats.error_status_code` if an error occurs during initialization or execution.
pub fn canonicalize(
    graph: &[usize],
    labeling: &mut [u32],
    partition: &mut [u32],
    active: Option<&[usize]>,
    orbits: &mut [u32],
    options: &mut CanonautOptions,
    stats: &mut Statistics,
    workspace: &mut [usize],
    words_per_vertex: u32,
    number_of_vertices: u32,
    canonical_graph: Option<&mut [usize]>,
    graph_context: &mut GraphContext,
) {
    let mut vertex_index: u32;
    let mut number_of_cells: u32;

    if options.get_canonical && canonical_graph.is_none() {
        stats.error_status_code = 3;
        eprintln!("Canonical graph requested but no graph provided\n",);
        return;
    }

    let canonical_graph = canonical_graph.unwrap_or(&mut []);

    let initstatus: u32;

    options.algorithm_context.dispatch = options.dispatch;

    // Set up refine function
    if (options.user_refinement_procedure).is_some() {
        options.algorithm_context.dispatch.refine = options.user_refinement_procedure;
    } else if (options.algorithm_context.dispatch.refine1).is_some() && words_per_vertex == 1 {
        options.algorithm_context.dispatch.refine = options.algorithm_context.dispatch.refine1;
    }

    if !validate_dispatch_vector(&options.algorithm_context.dispatch) {
        eprintln!(">E bad dispatch vector");
        panic!();
    }

    if let Some(error_msg) = validate_input(words_per_vertex, number_of_vertices) {
        stats.error_status_code = if words_per_vertex > ALGORITHM_INFINITY / 64 + 1 {
            2
        } else {
            1
        };
        eprintln!("{error_msg}");
        return;
    }

    // Handle empty graph case
    if number_of_vertices == 0 {
        init_empty_stats(&mut *stats);
        initstatus = 0;
        if initstatus != 0 {
            stats.error_status_code = initstatus;
        }
        return;
    }

    // Initialize options.algorithm_context from parameters and options
    options.algorithm_context.words_per_vertex = words_per_vertex;
    options.algorithm_context.number_of_vertices = number_of_vertices;
    // options.algorithm_context.orbits = orbits;

    // Configuration from options
    options.algorithm_context.get_canonical = options.get_canonical;
    options.algorithm_context.digraph = options.digraph;
    options.algorithm_context.write_automorphisms = options.write_automorphisms;
    options.algorithm_context.do_markers = options.write_markers;
    options.algorithm_context.cartesian = options.cartesian;
    options.algorithm_context.use_schreier_sims = options.schreier;
    options.algorithm_context.line_length = options.line_length;
    options.algorithm_context.target_cell_level = if options.algorithm_context.digraph {
        0
    } else {
        options.target_cell_level
    };
    options.algorithm_context.outfile = match options.outfile.as_mut() {
        Some(stream) => std::mem::take(stream),
        None => OutputStream::StandardOutput,
    };

    // Function pointers
    options.algorithm_context.user_node_procedure = options.user_node_procedure;
    options.algorithm_context.user_automorphism_procedure = options.user_automorphism_procedure;
    options.algorithm_context.user_level_procedure = options.user_level_procedure;
    options.algorithm_context.user_canonical_procedure = options.user_canonical_procedure;
    options.algorithm_context.invariant_procedure = options.invariant_procedure;

    // Invariant levels
    if options.min_invariant_level < 0 && options.get_canonical {
        options.algorithm_context.min_invariant_level = -options.min_invariant_level;
    } else {
        options.algorithm_context.min_invariant_level = options.min_invariant_level;
    }
    if options.max_invariant_level < 0 && options.get_canonical {
        options.algorithm_context.max_invariant_level = -options.max_invariant_level;
    } else {
        options.algorithm_context.max_invariant_level = options.max_invariant_level;
    }
    options.algorithm_context.invariant_argument = options.invariant_argument;

    check_build_consistency(64, 28090);

    // Ensure all internal buffers have sufficient capacity.  Buffers only grow
    // so they are safely reused across repeated calls (no allocation on the hot path).
    options.algorithm_context.ensure_capacity(
        options.algorithm_context.words_per_vertex,
        options.algorithm_context.number_of_vertices,
    );
    graph_context.ensure_workperm_size(options.algorithm_context.number_of_vertices);
    graph_context.ensure_workset_size(options.algorithm_context.words_per_vertex);
    graph_context.ensure_bucket_size(options.algorithm_context.number_of_vertices);

    // Reset the bump allocator so this run reuses the same O(depth) nodes.
    options.algorithm_context.target_cell_node_pool.reset();

    // Initialize partition
    if options.default_partition {
        let vertex_count = options.algorithm_context.number_of_vertices as usize;
        // Separate fills so the compiler can SIMD-vectorize each independently.
        partition[..vertex_count].fill(ALGORITHM_INFINITY);
        partition[vertex_count - 1] = 0;
        for (position, label) in labeling[..vertex_count].iter_mut().enumerate() {
            *label = position as u32;
        }
        options
            .algorithm_context
            .active
            .clear_words(words_per_vertex as usize);
        options.algorithm_context.active.set_bit(0);
        number_of_cells = 1;
    } else {
        partition[(options.algorithm_context.number_of_vertices - 1) as usize] = 0;

        number_of_cells = 0;
        for cell_end_marker in partition
            .iter_mut()
            .take(options.algorithm_context.number_of_vertices as usize)
        {
            if *cell_end_marker != 0 {
                *cell_end_marker = ALGORITHM_INFINITY;
            } else {
                number_of_cells += 1;
            }
        }
        if let Some(active) = active {
            let word_count = words_per_vertex as usize;
            options.algorithm_context.active[..word_count].copy_from_slice(&active[..word_count]);
        } else {
            options
                .algorithm_context
                .active
                .clear_words(words_per_vertex as usize);
            vertex_index = 0;
            while vertex_index < options.algorithm_context.number_of_vertices {
                options.algorithm_context.active.set_bit(vertex_index as usize);
                while partition[vertex_index as usize] != 0 {
                    vertex_index += 1;
                }
                vertex_index += 1;
            }
        }
    }

    // options.algorithm_graph = &mut [];
    initstatus = 0;

    // if (context.dispatch.init).is_some() {
    //     (context.dispatch.init).expect("non-null function pointer")(
    //         graph,
    //         &mut graph,
    //         canonical_graph,
    //         &mut context.canonical_graph,
    //         labeling.as_mut_ptr(),
    //         partition.as_mut_ptr(),
    //         active.unwrap().as_mut_ptr(),
    //         options,
    //         &mut initstatus,
    //         words_per_vertex,
    //         number_of_vertices,
    //     );
    // }

    if initstatus != 0 {
        stats.error_status_code = initstatus;
        return;
    }

    if options.algorithm_context.use_schreier_sims {
        options.algorithm_context.schreiersims.reset(options.algorithm_context.number_of_vertices as usize);
    }

    // Initialize orbits and stats
    for (vertex_index, orbit) in orbits.iter_mut().enumerate() {
        *orbit = vertex_index as u32;
    }

    stats.group_size1 = 1.0f64;
    stats.group_size2 = 0;
    stats.number_of_generators = 0;
    stats.number_of_nodes = 0;
    stats.number_of_bad_leaves = 0;
    stats.total_target_cell_members = 0;
    stats.canonical_updates = 0;
    stats.number_of_orbits = options.algorithm_context.number_of_vertices;

    options
        .algorithm_context
        .fixed_points
        .clear_words(words_per_vertex as usize);

    // Initialize algorithm state
    options.algorithm_context.first_nontrivial_level = 1;
    options.algorithm_context.canonical_match_level = -(1);
    options.algorithm_context.invariant_success_level = ALGORITHM_INFINITY as i32;
    options.algorithm_context.invariant_successes = 0;
    options.algorithm_context.invariant_applications = options.algorithm_context.invariant_successes;
    options.algorithm_context.needs_short_prune = false;

    // Set up workspace.  Prefer the caller-supplied workspace if it is larger than
    // the one already owned by algorithm_context (e.g. standalone canonicalize() call with
    // a big workspace).  Otherwise reuse algorithm_context.workspace as-is — it was
    // pre-sized by ensure_capacity so no allocation is needed.
    if workspace.len() > options.algorithm_context.workspace.len() {
        options.algorithm_context.workspace = workspace.to_vec();
    } else if options.algorithm_context.workspace.len()
        < (2 * options.algorithm_context.words_per_vertex) as usize
    {
        // Fallback: if somehow workspace is still too small, use default_workspace
        options.algorithm_context.workspace =
            std::mem::take(&mut options.algorithm_context.default_workspace);
    }
    options.algorithm_context.workspace_top = options.algorithm_context.workspace.len()
        - options.algorithm_context.workspace.len()
            % (2 * options.algorithm_context.words_per_vertex) as usize;
    options.algorithm_context.free_memory_pointer = 0;

    stats.error_status_code = 0;

    // Start main algorithm - using root pool index instead of TargetCellNode pointer
    let backtrack_level: i32 = first_path_node_0(
        graph,
        orbits,
        canonical_graph,
        labeling,
        partition,
        1,
        number_of_cells,
        None,
        &mut options.algorithm_context,
        stats,
        graph_context,
    );

    if backtrack_level == -(11) {
        stats.error_status_code = 4;
    } else if backtrack_level == -(12) {
        stats.error_status_code = 5;
    } else {
        if options.algorithm_context.get_canonical {
            (options.algorithm_context.dispatch.update_canonical_graph).expect("non-null function pointer")(
                graph,
                canonical_graph,
                &options.algorithm_context.canonical_labeling,
                options.algorithm_context.canonical_matching_rows,
                words_per_vertex,
                number_of_vertices,
                graph_context,
            );
            vertex_index = 0;
            while vertex_index < options.algorithm_context.number_of_vertices {
                labeling[vertex_index as usize] =
                    options.algorithm_context.canonical_labeling[vertex_index as usize];
                vertex_index += 1;
            }
        }
        stats.shallowest_invariant_success_level =
            if options.algorithm_context.invariant_success_level == ALGORITHM_INFINITY as i32 {
                0
            } else {
                options.algorithm_context.invariant_success_level
            };
        stats.invariant_procedure_calls = options.algorithm_context.invariant_applications;
        stats.successful_invariant_calls = options.algorithm_context.invariant_successes;
    }

    // Cleanup
    if options.algorithm_context.number_of_vertices >= 320 {
        if let Some(procedure) = options.algorithm_context.dispatch.freedyn {
            procedure();
        }
        search_freedyn();
    }

}

/// Processes one node along the first (leftmost) path of the search tree.
/// `parent_node_index` indexes into `context.target_cell_node_pool`; returns
/// the backtrack level.
fn first_path_node_0(
    graph: &[usize],
    orbits: &mut [u32],
    canonical_graph: &mut [usize],
    labeling: &mut [u32],
    partition: &mut [u32],
    level: i32,
    mut number_of_cells: u32,
    parent_node_index: Option<usize>,
    context: &mut AlgorithmContext,
    stats_ptr: &mut Statistics,
    graph_context: &mut GraphContext,
) -> i32 {
    let mut target_vertex: Option<usize>;

    let mut same_orbit_count: u32;
    let mut backtrack_level: i32;
    let mut target_cell_size: u32 = 0;
    let mut target_cell_index: i32;
    let mut child_count: u32 = 0;
    let mut quick_invariant: u32 = 0;
    let mut hash_code: u32 = 0;

    // Get or allocate node using pure Rust TargetCellNodePool
    let current_node_index = {
        // Check if we have a next node from parent
        let existing_index = if let Some(parent_idx) = parent_node_index {
            context.target_cell_node_pool.get_next_index(parent_idx)
        } else {
            None
        };

        if let Some(existing_idx) = existing_index {
            existing_idx
        } else {
            // Allocate new node
            let new_index = context
                .target_cell_node_pool
                .allocate_node(context.words_per_vertex as usize);

            // Link to parent if we have one
            if let Some(parent_idx) = parent_node_index {
                context
                    .target_cell_node_pool
                    .link_nodes(parent_idx, new_index);
            }
            new_index
        }
    };

    stats_ptr.number_of_nodes += 1;

    apply_refinement(
        graph,
        labeling,
        partition,
        level,
        &mut number_of_cells,
        &mut quick_invariant,
        &mut context.work_permutation,
        &mut context.active,
        &mut hash_code,
        context.dispatch.refine,
        context.invariant_procedure,
        context.min_invariant_level,
        context.max_invariant_level,
        context.invariant_argument,
        context.digraph,
        context.words_per_vertex,
        context.number_of_vertices,
        graph_context,
    );

    context.first_path_codes[level as usize] = hash_code as i16;
    if quick_invariant > 0 {
        context.invariant_applications += 1;
        if quick_invariant == 2 {
            context.invariant_successes += 1;
            if context.min_invariant_level < 0 {
                context.min_invariant_level = level;
            }
            if context.max_invariant_level < 0 {
                context.max_invariant_level = level;
            }
            if level < context.invariant_success_level {
                context.invariant_success_level = level;
            }
        }
    }

    target_cell_index = -(1);
    if number_of_cells != context.number_of_vertices {
        // SAFETY: current_node_index was just allocated/retrieved from the pool
        let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };

        select_target_cell(
            graph,
            labeling,
            partition,
            level,
            &mut target_cell.cell_data,
            &mut target_cell_size,
            &mut target_cell_index,
            context.target_cell_level,
            context.digraph,
            -(1),
            context.dispatch.target_cell,
            context.words_per_vertex,
            context.number_of_vertices,
            graph_context,
        );
        stats_ptr.total_target_cell_members += target_cell_size as usize;
    }

    context.first_target_cell[level as usize] = target_cell_index;
    if let Some(procedure) = context.user_node_procedure {
        procedure(
            graph,
            labeling,
            partition,
            level as u32,
            number_of_cells,
            target_cell_index,
            context.first_path_codes[level as usize] as u32,
            context.words_per_vertex,
            context.number_of_vertices,
        );
    }

    if number_of_cells == context.number_of_vertices {
        firstterminal(labeling, level as u32, context, stats_ptr);
        if let Some(procedure) = context.user_level_procedure {
            procedure(
                labeling,
                partition,
                level as u32,
                orbits,
                stats_ptr,
                0,
                1,
                1,
                context.number_of_vertices,
                0,
                context.number_of_vertices,
            );
        }
        if let Some(canonical_procedure) =
            context.user_canonical_procedure.filter(|_| context.get_canonical)
        {
            (context.dispatch.update_canonical_graph).expect("non-null function pointer")(
                graph,
                canonical_graph,
                &context.canonical_labeling,
                context.canonical_matching_rows,
                context.words_per_vertex,
                context.number_of_vertices,
                graph_context,
            );
            context.canonical_matching_rows = context.number_of_vertices;
            if canonical_procedure(
                graph,
                &mut context.canonical_labeling,
                canonical_graph,
                stats_ptr.canonical_updates,
                context.canonical_path_codes[level as usize] as u32,
                context.words_per_vertex,
                context.number_of_vertices,
            ) != 0
            {
                return -(11);
            }
        }
        return level - 1;
    }

    if kill_requested() {
        return -(12);
    }

    if context.first_nontrivial_level >= level as u32
        && !(context.dispatch.cheap_automorphism_check).expect("non-null function pointer")(
            partition,
            level as u32,
            context.digraph,
            context.number_of_vertices,
        )
    {
        context.first_nontrivial_level = level as u32 + 1;
    }

    same_orbit_count = 0;
    // SAFETY: current_node_index was just allocated/retrieved from the pool
    let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
    target_vertex = target_cell.cell_data.next_set_bit(0);
    let first_target_vertex: Option<usize> = target_vertex;
    while let Some(target_vertex_value) = target_vertex {
        if orbits[target_vertex_value] == target_vertex_value as u32 {
            individualize_vertex(
                labeling,
                partition,
                level + 1,
                target_cell_index,
                target_vertex_value as u32,
                &mut context.active,
                context.words_per_vertex,
            );
            context.fixed_points.set_bit(target_vertex_value);
            context.coset_vertex = target_vertex_value as u32;
            if target_vertex == first_target_vertex {
                backtrack_level = first_path_node_0(
                    graph,
                    orbits,
                    canonical_graph,
                    labeling,
                    partition,
                    level + 1,
                    number_of_cells + 1,
                    Some(current_node_index),
                    context,
                    stats_ptr,
                    graph_context,
                );
                child_count = 1;
                context.gca_first_level = level as u32;
                context.orbit_representative = first_target_vertex.unwrap() as u32;
            } else {
                backtrack_level = othernode0_rust(
                    graph,
                    orbits,
                    canonical_graph,
                    labeling,
                    partition,
                    level as u32 + 1,
                    number_of_cells + 1,
                    Some(current_node_index),
                    context,
                    stats_ptr,
                    graph_context,
                );
                child_count += 1;
            }
            context.fixed_points.unset_bit(target_vertex_value);
            if backtrack_level < level {
                return backtrack_level;
            }
            if context.needs_short_prune {
                context.needs_short_prune = false;
                let words = context.words_per_vertex as usize;
                let src_start = context.free_memory_pointer - words;
                // SAFETY: current_node_index was just allocated/retrieved from the pool
                let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
                simd_ops::intersect_inplace(
                    &mut target_cell.cell_data,
                    &context.workspace[src_start..src_start + words],
                    words,
                );
            }
            recover(partition, level as u32, context);
        }
        if orbits[target_vertex_value] == first_target_vertex.unwrap() as u32 {
            same_orbit_count += 1;
        }
        // SAFETY: current_node_index was just allocated/retrieved from the pool
        let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
        target_vertex = target_cell.cell_data.next_set_bit(target_vertex_value + 1);
    }

    multiply(
        &mut stats_ptr.group_size1,
        &mut stats_ptr.group_size2,
        same_orbit_count,
    );

    if target_cell_size == same_orbit_count && context.uniform_partition_depth == (level + 1) as u32 {
        context.uniform_partition_depth -= 1;
    }
    if context.do_markers {
        writemarker(
            level as u32,
            match first_target_vertex {
                Some(vertex) => vertex as i32,
                None => -1,
            },
            same_orbit_count,
            target_cell_size,
            stats_ptr.number_of_orbits,
            number_of_cells,
            context.number_of_vertices,
            &mut context.outfile,
        );
    }

    if let Some(procedure) = context.user_level_procedure {
        procedure(
            labeling,
            partition,
            level as u32,
            orbits,
            stats_ptr,
            first_target_vertex.unwrap() as u32,
            same_orbit_count,
            target_cell_size,
            number_of_cells,
            child_count,
            context.number_of_vertices,
        );
    }
    level - 1
}

/// Processes a node off the first path — same shape as `first_path_node_0`, but
/// this branch of the tree isn't guaranteed canonical yet, so it runs the
/// automorphism/pruning checks that the first path skips.
fn othernode0_rust(
    graph: &[usize],
    orbits: &mut [u32],
    canonical_graph: &mut [usize],
    labeling: &mut [u32],
    partition: &mut [u32],
    level: u32,
    mut number_of_cells: u32,
    parent_node_index: Option<usize>,
    context: &mut AlgorithmContext,
    stats_ptr: &mut Statistics,
    graph_context: &mut GraphContext,
) -> i32 {
    let mut target_vertex: Option<usize>;

    let mut hash_code: u32 = 0;
    let mut backtrack_level: i32;
    let mut target_cell_size: u32 = 0;
    let mut target_cell_index: i32;
    let mut quick_invariant: u32 = 0;

    // Get or allocate node using pure Rust TargetCellNodePool
    let current_node_index = {
        // Check if we have a next node from parent
        let existing_index = if let Some(parent_idx) = parent_node_index {
            context.target_cell_node_pool.get_next_index(parent_idx)
        } else {
            None
        };

        if let Some(existing_idx) = existing_index {
            existing_idx
        } else {
            // Allocate new node
            let new_index = context
                .target_cell_node_pool
                .allocate_node(context.words_per_vertex as usize);

            // Link to parent if we have one
            if let Some(parent_idx) = parent_node_index {
                context
                    .target_cell_node_pool
                    .link_nodes(parent_idx, new_index);
            }

            new_index
        }
    };

    if kill_requested() {
        return -(12);
    }

    stats_ptr.number_of_nodes += 1;
    apply_refinement(
        graph,
        labeling,
        partition,
        level as i32,
        &mut number_of_cells,
        &mut quick_invariant,
        &mut context.work_permutation,
        &mut context.active,
        &mut hash_code,
        context.dispatch.refine,
        context.invariant_procedure,
        context.min_invariant_level,
        context.max_invariant_level,
        context.invariant_argument,
        context.digraph,
        context.words_per_vertex,
        context.number_of_vertices,
        graph_context,
    );

    let path_code: i16 = hash_code as i16;
    if quick_invariant > 0 {
        context.invariant_applications = context.invariant_applications.wrapping_add(1);
        if quick_invariant == 2 {
            context.invariant_successes = context.invariant_successes.wrapping_add(1);
            if (level as i32) < context.invariant_success_level {
                context.invariant_success_level = level as i32;
            }
        }
    }

    if context.first_path_match_level == level - 1 && path_code as u32 == context.first_path_codes[level as usize] as u32 {
        context.first_path_match_level = level;
    }

    if context.get_canonical {
        if context.canonical_match_level == (level - 1) as i32 {
            if (path_code as u32) < context.canonical_path_codes[level as usize] as u32 {
                context.canonical_comparison = -(1);
            } else if path_code as u32 > context.canonical_path_codes[level as usize] as u32 {
                context.canonical_comparison = 1;
            } else {
                context.canonical_comparison = 0;
                context.canonical_match_level = level as i32;
            }
        }
        if context.canonical_comparison > 0 {
            context.canonical_path_codes[level as usize] = path_code;
        }
    }

    target_cell_index = -(1);
    if number_of_cells < context.number_of_vertices
        && (context.first_path_match_level == level || context.get_canonical && context.canonical_comparison >= 0)
    {
        if !context.get_canonical || context.canonical_comparison < 0 {
            // SAFETY: current_node_index was just allocated/retrieved from the pool
            let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
            select_target_cell(
                graph,
                labeling,
                partition,
                level as i32,
                &mut target_cell.cell_data,
                &mut target_cell_size,
                &mut target_cell_index,
                context.target_cell_level,
                context.digraph,
                context.first_target_cell[level as usize],
                context.dispatch.target_cell,
                context.words_per_vertex,
                context.number_of_vertices,
                graph_context,
            );
            if target_cell_index != context.first_target_cell[level as usize] {
                context.first_path_match_level = level - 1;
            }
        } else {
            // SAFETY: current_node_index was just allocated/retrieved from the pool
            let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
            select_target_cell(
                graph,
                labeling,
                partition,
                level as i32,
                &mut target_cell.cell_data,
                &mut target_cell_size,
                &mut target_cell_index,
                context.target_cell_level,
                context.digraph,
                -(1),
                context.dispatch.target_cell,
                context.words_per_vertex,
                context.number_of_vertices,
                graph_context,
            );
        }
        stats_ptr.total_target_cell_members =
            (stats_ptr.total_target_cell_members).wrapping_add(target_cell_size as usize);
    }

    if let Some(procedure) = context.user_node_procedure {
        procedure(
            graph,
            labeling,
            partition,
            level,
            number_of_cells,
            target_cell_index,
            path_code as u32,
            context.words_per_vertex,
            context.number_of_vertices,
        );
    }

    backtrack_level = processnode(
        graph,
        orbits,
        canonical_graph,
        labeling,
        partition,
        level,
        number_of_cells,
        context,
        stats_ptr,
        graph_context,
    );
    if backtrack_level < level as i32 {
        return backtrack_level;
    }

    if context.needs_short_prune {
        context.needs_short_prune = false;
        let words = context.words_per_vertex as usize;
        let src_start = context.free_memory_pointer - words;
        // SAFETY: current_node_index was just allocated/retrieved from the pool
        let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
        simd_ops::intersect_inplace(
            &mut target_cell.cell_data,
            &context.workspace[src_start..src_start + words],
            words,
        );
    }

    if !(context.dispatch.cheap_automorphism_check).expect("non-null function pointer")(
        partition,
        level,
        context.digraph,
        context.number_of_vertices,
    ) {
        context.first_nontrivial_level = level + 1;
    }

    // SAFETY: current_node_index was just allocated/retrieved from the pool
    let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
    target_vertex = target_cell.cell_data.next_set_bit(0);
    let first_target_vertex: Option<usize> = target_vertex;
    while let Some(target_vertex_value) = target_vertex {
        individualize_vertex(
            labeling,
            partition,
            level as i32 + 1,
            target_cell_index,
            target_vertex_value as u32,
            &mut context.active,
            context.words_per_vertex,
        );
        context.fixed_points.set_bit(target_vertex_value);
        backtrack_level = othernode0_rust(
            graph,
            orbits,
            canonical_graph,
            labeling,
            partition,
            level + 1,
            number_of_cells + 1,
            Some(current_node_index),
            context,
            stats_ptr,
            graph_context,
        );
        context.fixed_points.unset_bit(target_vertex_value);
        if backtrack_level < level as i32 {
            return backtrack_level;
        }
        if context.needs_short_prune {
            context.needs_short_prune = false;
            let words = context.words_per_vertex as usize;
            let src_start = context.free_memory_pointer - words;
            // SAFETY: current_node_index was just allocated/retrieved from the pool
            let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
            simd_ops::intersect_inplace(
                &mut target_cell.cell_data,
                &context.workspace[src_start..src_start + words],
                words,
            );
        }
        if target_vertex == first_target_vertex {
            // SAFETY: current_node_index was just allocated/retrieved from the pool
            let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
            prune_by_generators(
                &mut target_cell.cell_data,
                &context.fixed_points,
                &context.workspace[..context.free_memory_pointer],
                context.words_per_vertex,
            );
            if context.use_schreier_sims {
                // SAFETY: current_node_index was just allocated/retrieved from the pool
                let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
                context
                    .schreiersims
                    .prune_set(&context.fixed_points, &mut target_cell.cell_data);
            }
        }
        recover(partition, level, context);
        // SAFETY: current_node_index was just allocated/retrieved from the pool
        let target_cell = unsafe { context.target_cell_node_pool.get_node_unchecked_mut(current_node_index) };
        target_vertex = target_cell.cell_data.next_set_bit(target_vertex_value + 1);
    }
    (level - 1) as i32
}

/// Records the labeling at the first terminal node reached — becomes the
/// reference (`first_labeling`/`first_path_codes`) that later branches compare against.
fn firstterminal(
    labeling: &mut [u32],
    level: u32,
    context: &mut AlgorithmContext,
    stats_ptr: &mut Statistics,
) {
    let mut vertex_index: u32;

    stats_ptr.maximum_search_depth = level;
    context.first_path_match_level = level;
    context.uniform_partition_depth = context.first_path_match_level;
    context.gca_first_level = context.uniform_partition_depth;
    context.first_path_codes[(level + 1) as usize] = 0o77777_u32 as i16;
    context.first_target_cell[(level + 1) as usize] = -(1);

    vertex_index = 0;
    while vertex_index < context.number_of_vertices {
        context.first_labeling[vertex_index as usize] = labeling[vertex_index as usize];
        vertex_index += 1;
    }

    if context.get_canonical {
        context.gca_canonical_level = level as i32;
        context.canonical_match_level = context.gca_canonical_level;
        context.canonical_depth = context.canonical_match_level as u32;
        context.canonical_comparison = 0;
        context.canonical_matching_rows = 0;

        vertex_index = 0;
        while vertex_index < context.number_of_vertices {
            context.canonical_labeling[vertex_index as usize] = labeling[vertex_index as usize];
            vertex_index += 1;
        }

        vertex_index = 0;
        while vertex_index <= level {
            context.canonical_path_codes[vertex_index as usize] = context.first_path_codes[vertex_index as usize];
            vertex_index += 1;
        }
        context.canonical_path_codes[(level + 1) as usize] = 0o77777_u32 as i16;
        stats_ptr.canonical_updates = 1;
    }
}

/// Classifies a terminal node (found generator / canonical improved / pruned /
/// no-op — the `node_outcome` values below) and dispatches the matching bookkeeping.
/// Returns the backtrack level.
#[inline]
fn processnode(
    graph: &[usize],
    mut orbits: &mut [u32],
    canonical_graph: &mut [usize],
    labeling: &mut [u32],
    partition: &mut [u32],
    level: u32,
    number_of_cells: u32,
    context: &mut AlgorithmContext,
    stats_ptr: &mut Statistics,
    graph_context: &mut GraphContext,
) -> i32 {
    let mut vertex_index: u32;
    let mut node_outcome: u32;

    let short_prune_applicable: bool;
    let mut matching_rows: u32 = 0;

    node_outcome = 0;
    if context.first_path_match_level != level && (!context.get_canonical || context.canonical_comparison < 0) {
        node_outcome = 4;
    } else if number_of_cells == context.number_of_vertices {
        if context.first_path_match_level == level {
            vertex_index = 0;
            while vertex_index < context.number_of_vertices {
                context.work_permutation[context.first_labeling[vertex_index as usize] as usize] = labeling[vertex_index as usize];
                vertex_index += 1;
            }
            if context.gca_first_level >= context.first_nontrivial_level
                || (context.dispatch.is_automorphism).expect("non-null function pointer")(
                    graph,
                    &context.work_permutation,
                    context.digraph,
                    context.words_per_vertex,
                    context.number_of_vertices,
                )
            {
                node_outcome = 1;
            }
        }
        if node_outcome == 0 {
            if context.get_canonical {
                matching_rows = 0;
                if context.canonical_comparison == 0 {
                    if level < context.canonical_depth {
                        context.canonical_comparison = 1;
                    } else {
                        (context.dispatch.update_canonical_graph).expect("non-null function pointer")(
                            graph,
                            canonical_graph,
                            &context.canonical_labeling,
                            context.canonical_matching_rows,
                            context.words_per_vertex,
                            context.number_of_vertices,
                            graph_context,
                        );
                        context.canonical_matching_rows = context.number_of_vertices;
                        context.canonical_comparison = (context.dispatch.test_canonical_labeling)
                            .expect("non-null function pointer")(
                            graph,
                            canonical_graph,
                            labeling,
                            &mut matching_rows,
                            context.words_per_vertex,
                            context.number_of_vertices,
                            graph_context,
                        );
                    }
                }
                if context.canonical_comparison == 0 {
                    vertex_index = 0;
                    while vertex_index < context.number_of_vertices {
                        context.work_permutation[context.canonical_labeling[vertex_index as usize] as usize] =
                            labeling[vertex_index as usize];
                        vertex_index += 1;
                    }
                    node_outcome = 2;
                } else if context.canonical_comparison > 0 {
                    node_outcome = 3;
                } else {
                    node_outcome = 4;
                }
            } else {
                node_outcome = 4;
            }
        }
    }

    if node_outcome != 0 && level > stats_ptr.maximum_search_depth {
        stats_ptr.maximum_search_depth = level;
    }

    match node_outcome {
        0 => return level as i32,
        1 => {
            if context.free_memory_pointer == context.workspace_top {
                context.free_memory_pointer -= (2 * context.words_per_vertex) as usize;
            }
            {
                let workspace_offset = context.free_memory_pointer;
                let word_count = context.words_per_vertex as usize;
                let vertex_count = context.number_of_vertices;
                let (fixed_points_slice, min_cycle_representatives_slice) =
                    context.workspace[workspace_offset..workspace_offset + 2 * word_count].split_at_mut(word_count);
                fmperm_safe(&context.work_permutation[..vertex_count as usize], fixed_points_slice, min_cycle_representatives_slice, context.words_per_vertex, vertex_count);
            }
            context.free_memory_pointer += (2 * context.words_per_vertex) as usize;
            if context.write_automorphisms {
                write_permutation(
                    &mut context.outfile,
                    &context.work_permutation,
                    context.cartesian,
                    context.line_length,
                    context.number_of_vertices,
                );
            }
            let orbits_slice = &mut orbits;
            let workperm_slice = &mut context.work_permutation[..context.number_of_vertices as usize];
            stats_ptr.number_of_orbits =
                safe_orbjoin(orbits_slice, workperm_slice, context.number_of_vertices);
            stats_ptr.number_of_generators += 1;
            if let Some(procedure) = context.user_automorphism_procedure {
                procedure(
                    stats_ptr.number_of_generators,
                    &mut context.work_permutation,
                    orbits,
                    stats_ptr.number_of_orbits,
                    context.orbit_representative,
                    context.number_of_vertices,
                );
            }
            if context.use_schreier_sims {
                context
                    .schreiersims
                    .add_generator(&context.work_permutation[..context.number_of_vertices as usize]);
            }
            // This automorphism reaches the first leaf. Backtracking relative
            // to the canonical leaf can skip branches and lose generators.
            return context.gca_first_level as i32;
        }
        2 => {
            if context.free_memory_pointer == context.workspace_top {
                context.free_memory_pointer -= (2 * context.words_per_vertex) as usize;
            }
            {
                let workspace_offset = context.free_memory_pointer;
                let word_count = context.words_per_vertex as usize;
                let vertex_count = context.number_of_vertices;
                let (fixed_points_slice, min_cycle_representatives_slice) =
                    context.workspace[workspace_offset..workspace_offset + 2 * word_count].split_at_mut(word_count);
                fmperm_safe(&context.work_permutation[..vertex_count as usize], fixed_points_slice, min_cycle_representatives_slice, context.words_per_vertex, vertex_count);
            }
            context.free_memory_pointer += (2 * context.words_per_vertex) as usize;
            let prior_orbit_count = stats_ptr.number_of_orbits as i32;
            let orbits_slice = &mut orbits;
            let workperm_slice = &mut context.work_permutation[..context.number_of_vertices as usize];
            stats_ptr.number_of_orbits =
                safe_orbjoin(orbits_slice, workperm_slice, context.number_of_vertices);
            if stats_ptr.number_of_orbits as i32 == prior_orbit_count {
                if context.gca_canonical_level != context.gca_first_level as i32 {
                    context.needs_short_prune = true;
                }
                return context.gca_canonical_level;
            }
            if context.write_automorphisms {
                write_permutation(
                    &mut context.outfile,
                    &context.work_permutation,
                    context.cartesian,
                    context.line_length,
                    context.number_of_vertices,
                );
            }
            stats_ptr.number_of_generators += 1;
            if let Some(procedure) = context.user_automorphism_procedure {
                procedure(
                    stats_ptr.number_of_generators,
                    &mut context.work_permutation,
                    orbits,
                    stats_ptr.number_of_orbits,
                    context.orbit_representative,
                    context.number_of_vertices,
                );
            }
            if context.use_schreier_sims {
                context
                    .schreiersims
                    .add_generator(&context.work_permutation[..context.number_of_vertices as usize]);
            }
            if orbits[context.coset_vertex as usize] < context.coset_vertex {
                return context.gca_first_level as i32;
            }
            if context.gca_canonical_level != context.gca_first_level as i32 {
                context.needs_short_prune = true;
            }
            return context.gca_canonical_level;
        }
        3 => {
            stats_ptr.canonical_updates = (stats_ptr.canonical_updates).wrapping_add(1);
            vertex_index = 0;
            while vertex_index < context.number_of_vertices {
                context.canonical_labeling[vertex_index as usize] = labeling[vertex_index as usize];
                vertex_index += 1;
            }
            context.gca_canonical_level = level as i32;
            context.canonical_match_level = context.gca_canonical_level;
            context.canonical_depth = context.canonical_match_level as u32;
            context.canonical_comparison = 0;
            context.canonical_path_codes[(level + 1) as usize] = 0o77777_u32 as i16;
            context.canonical_matching_rows = matching_rows;
            if let Some(canonical_procedure) =
            context.user_canonical_procedure.filter(|_| context.get_canonical)
        {
                (context.dispatch.update_canonical_graph).expect("non-null function pointer")(
                    graph,
                    canonical_graph,
                    &context.canonical_labeling,
                    context.canonical_matching_rows,
                    context.words_per_vertex,
                    context.number_of_vertices,
                    graph_context,
                );
                context.canonical_matching_rows = context.number_of_vertices;
                if canonical_procedure(
                    graph,
                    &mut context.canonical_labeling,
                    canonical_graph,
                    stats_ptr.canonical_updates,
                    context.canonical_path_codes[level as usize] as u32,
                    context.words_per_vertex,
                    context.number_of_vertices,
                ) != 0
                {
                    return -(11);
                }
            }
        }
        4 => {
            stats_ptr.number_of_bad_leaves = (stats_ptr.number_of_bad_leaves).wrapping_add(1);
        }
        _ => {}
    }

    if level != context.first_nontrivial_level {
        short_prune_applicable = true;
        if context.free_memory_pointer == context.workspace_top {
            context.free_memory_pointer -= (2 * context.words_per_vertex) as usize;
        }
        // Mutably split at index to make this work safely
        let (first, second) = context
            .workspace
            .split_at_mut(context.free_memory_pointer + context.words_per_vertex as usize);

        mark_partition_fixed_points(
            labeling,
            partition,
            context.first_nontrivial_level,
            &mut first[context.free_memory_pointer..],
            second,
            context.words_per_vertex,
            context.number_of_vertices,
        );
        context.free_memory_pointer += (2 * context.words_per_vertex) as usize;
    } else {
        short_prune_applicable = false;
    }

    let candidate_level = if context.uniform_partition_depth as i32 > context.canonical_match_level {
        context.uniform_partition_depth as i32 - 1
    } else {
        context.canonical_match_level
    };
    let new_backtrack_level: i32 = if context.first_nontrivial_level as i32 <= candidate_level {
        context.first_nontrivial_level as i32 - 1
    } else {
        candidate_level
    };
    if short_prune_applicable && new_backtrack_level != context.gca_first_level as i32 {
        context.needs_short_prune = true;
    }
    new_backtrack_level
}

/// Resets partition entries finer than `level` after backtracking, and pulls
/// back any of `context`'s per-level tracking fields the backtrack invalidated.
#[inline]
fn recover(partition: &mut [u32], level: u32, context: &mut AlgorithmContext) {
    // SAFETY: partition has length number_of_vertices, vertex_index < number_of_vertices.
    let vertex_count = context.number_of_vertices as usize;
    for vertex_index in 0..vertex_count {
        let partition_entry = unsafe { partition.get_unchecked_mut(vertex_index) };
        if *partition_entry > level {
            *partition_entry = ALGORITHM_INFINITY;
        }
    }

    if level < context.first_nontrivial_level {
        context.first_nontrivial_level = level + 1;
    }
    if level < context.first_path_match_level {
        context.first_path_match_level = level;
    }
    if context.get_canonical {
        if (level as i32) < context.gca_canonical_level {
            context.gca_canonical_level = level as i32;
        }
        if level as i32 <= context.canonical_match_level {
            context.canonical_match_level = level as i32;
            context.canonical_comparison = 0;
        }
    }
}

/// Writes one `-m` progress line (level, cell/orbit counts, fixed vertex) to `output`.
fn writemarker(
    level: u32,
    target_vertex: i32,
    same_orbit_count: u32,
    target_cell_size: u32,
    number_of_orbits: u32,
    number_of_cells: u32,
    mut _number_of_vertices: u32,
    output: &mut OutputStream,
) {
    let _ = output.write(b"level ");
    let number_string = level.to_string();
    let _ = output.write(number_string.as_bytes());
    let _ = output.write(b":  ");
    if number_of_cells != number_of_orbits {
        let number_string = number_of_cells.to_string();
        let _ = output.write(number_string.as_bytes());
        let _ = output.write(b" cell\0");
        if number_of_cells == 1 {
            let _ = output.write(b"; \0");
        } else {
            let _ = output.write(b"s; \0");
        }
    }
    let number_of_orbits_string = number_of_orbits.to_string();
    let _ = output.write(number_of_orbits_string.as_bytes());
    let _ = output.write(b" orbit");
    if number_of_orbits == 1 {
        let _ = output.write(b"; ");
    } else {
        let _ = output.write(b"s; ");
    }

    let target_vertex_string = (target_vertex + labelorg() as i32).to_string();
    let _ = output.write(target_vertex_string.as_bytes());
    let _ = output.write(b" fixed; index ");
    let index_string = same_orbit_count.to_string();
    let _ = output.write(index_string.as_bytes());
    if target_cell_size != same_orbit_count {
        let _ = output.write(b"/");
        let cell_size_string = target_cell_size.to_string();
        let _ = output.write(cell_size_string.as_bytes());
    }
    let _ = output.write(b"\n");
}

/// No-op: kept for API parity with nauty's C `freedyn`. Buffers live in
/// `AlgorithmContext` and are grown/reused by `ensure_capacity`, not freed here.
pub fn search_freedyn() {}

/// Helper function to handle line wrapping (logic from C macro CONDNL)
/// "writes end-of-line and 3 spaces if x characters won't fit"
fn conditional_newline(
    len_needed: u32,
    output: &mut OutputStream,
    current_length: &mut u32,
    linelength: u32,
) {
    if linelength > 0 && *current_length + len_needed > linelength {
        // C: putstring(f,"\n   "); curlen=3;
        let _ = output.write(b"\n   ");
        *current_length = 3;
    }
}

fn write_permutation(
    output: &mut OutputStream,
    permutation: &[u32],
    cartesian: bool,
    linelength: u32,
    vertex_count: u32,
) {
    let mut current_length = 0;

    if cartesian {
        for vertex_index in 0..vertex_count {
            let labeled_value = permutation[vertex_index as usize] + labelorg();
            let value_string = labeled_value.to_string();
            let value_string_len = value_string.len() as u32;

            conditional_newline(value_string_len + 1, output, &mut current_length, linelength);

            let _ = output.write(b" ");
            let _ = output.write(value_string.as_bytes());

            current_length += value_string_len + 1;
        }
        let _ = output.write(b"\n");
    } else {
        // C: DYNALLOC1 ... workperm ...
        // We use a vector to track visited elements for cycle decomposition
        let mut visited = vec![0; vertex_count as usize];

        for vertex_index in 0..vertex_count {
            let visited_index = vertex_index as usize;

            // If not visited yet AND not a fixed point (permutation[vertex_index] != vertex_index)
            if visited[visited_index] == 0 && permutation[visited_index] != vertex_index {
                let mut cycle_current_vertex = vertex_index;
                let cycle_start_label = cycle_current_vertex + labelorg();
                let cycle_start_string = cycle_start_label.to_string();
                let cycle_start_len = cycle_start_string.len() as u32;

                // C: if (curlen > 3) CONDNL(2*intlen+4);
                if current_length > 3 {
                    conditional_newline(
                        2 * cycle_start_len + 4,
                        output,
                        &mut current_length,
                        linelength,
                    );
                }

                let _ = output.write(b"(");

                // Do-while loop simulation
                loop {
                    let cycle_step_string = (cycle_current_vertex + labelorg()).to_string();
                    let cycle_step_len = cycle_step_string.len() as u32;

                    let _ = output.write(cycle_step_string.as_bytes());
                    current_length += cycle_step_len + 1;

                    let cycle_previous_vertex = cycle_current_vertex;
                    cycle_current_vertex = permutation[cycle_current_vertex as usize];
                    visited[cycle_previous_vertex as usize] = 1; // Mark visited

                    if cycle_current_vertex != vertex_index {
                        let cycle_next_string = (cycle_current_vertex + labelorg()).to_string();
                        let cycle_next_len = cycle_next_string.len() as u32;

                        // C: CONDNL(intlen+2);
                        conditional_newline(
                            cycle_next_len + 2,
                            output,
                            &mut current_length,
                            linelength,
                        );
                        let _ = output.write(b" ");
                    } else {
                        break;
                    }
                }

                let _ = output.write(b")");
                current_length += 1;
            }
        }

        // C: if (curlen == 0) putstring(f,"(1)\n");
        if current_length == 0 {
            let _ = output.write(b"(1)\n");
        } else {
            let _ = output.write(b"\n");
        }
    }
}

#[cfg(test)]
mod entry_tests {
    use super::*;
    use crate::refine::GraphContext;

    #[test]
    fn explicit_active_cells_use_words_not_vertices() {
        let mut options = CanonautOptions { default_partition: false, ..Default::default() };
        let mut stats = Statistics::new();
        let mut orbits = [0; 4];
        canonicalize(&[0; 4], &mut [0, 1, 2, 3], &mut [1, 1, 1, 0],
            Some(&[1usize << 63]), &mut orbits, &mut options, &mut stats,
            &mut [], 1, 4, None, &mut GraphContext::new(4));
        assert_eq!(stats.group_size1, 24.0);
        assert_eq!(orbits, [0; 4]);
    }

    #[test]
    fn missing_canonical_output_is_reported_on_the_first_call() {
        let mut options = CanonautOptions { get_canonical: true, ..Default::default() };
        let mut stats = Statistics::new();
        canonicalize(&[0], &mut [0], &mut [0], None, &mut [0],
            &mut options, &mut stats, &mut [], 1, 1, None, &mut GraphContext::new(1));
        assert_eq!(stats.error_status_code, 3);
    }
}
