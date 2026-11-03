use crate::bit_ops::{BIT_USIZE, W, WMASK, WSHIFT};
use crate::bit_ops::WORDSIZE;
use crate::io::OutputStream;
use crate::refine::GraphContext;
use crate::simd_ops::{andn_any, intersect_inplace};
use crate::structs::BitSetOpsMSB0;

use std::sync::atomic::{AtomicU32, Ordering};

/// Origin added to vertex labels in diagnostic output (nauty's `labelorg`).
pub static LABELORG: AtomicU32 = AtomicU32::new(0);
/// Set to a non-zero value to request early termination of a running search
/// (nauty's `nauty_kill_request`).
pub static KILL_REQUEST: AtomicU32 = AtomicU32::new(0);

/// Current label origin for diagnostic output.
#[inline]
pub fn labelorg() -> u32 {
    LABELORG.load(Ordering::Relaxed)
}

/// Whether early termination of the running search has been requested.
#[inline]
pub fn kill_requested() -> bool {
    KILL_REQUEST.load(Ordering::Relaxed) != 0
}

pub fn permute_set(set1: &[usize], set2: &mut [usize], words_per_vertex: u32, permutation: &[u32]) {
    // Clearing the lowest bit keeps iteration independent of the index lookup.
    match words_per_vertex {
        1 => permute_set_single_word(set1, set2, permutation),
        2 => permute_set_two_words(set1, set2, permutation),
        _ => permute_set_multiword(set1, set2, words_per_vertex, permutation),
    }
}

/// permute_set for words_per_vertex == 1 (n ≤ 64): no bounds checks, no trait dispatch.
#[inline(always)]
fn permute_set_single_word(set1: &[usize], set2: &mut [usize], permutation: &[u32]) {
    let mut permuted_word = 0usize;
    // SAFETY: caller guarantees set1.len() >= 1.
    let mut remaining_bits = unsafe { *set1.get_unchecked(0) };
    while remaining_bits != 0 {
        let source_bit_position = usize::BITS - 1 - remaining_bits.trailing_zeros();
        remaining_bits &= remaining_bits - 1;
        // SAFETY: padding bits are zero, so source_bit_position < n and the
        // permutation maps that vertex to another vertex below n <= 64.
        let mapped_bit_position = unsafe { *permutation.get_unchecked(source_bit_position as usize) };
        permuted_word |= 1usize << (usize::BITS - 1 - mapped_bit_position);
    }
    // SAFETY: caller guarantees set2.len() >= 1.
    unsafe { *set2.get_unchecked_mut(0) = permuted_word };
}

/// permute_set for words_per_vertex == 2 (64 < n ≤ 128): loads both words into
/// registers, writes both output words directly — no trait dispatch, no
/// bounds-checked `get_mut` per bit.
#[inline(always)]
fn permute_set_two_words(set1: &[usize], set2: &mut [usize], permutation: &[u32]) {
    let mut permuted_word0 = 0usize;
    let mut permuted_word1 = 0usize;

    // SAFETY: caller guarantees set1.len() >= 2.
    let mut source_word0 = unsafe { *set1.get_unchecked(0) };
    while source_word0 != 0 {
        let bit_index = (usize::BITS - 1 - source_word0.trailing_zeros()) as usize;
        source_word0 &= source_word0 - 1;
        // SAFETY: bit_index < 64 < n = permutation.len().
        let mapped_bit_position = unsafe { *permutation.get_unchecked(bit_index) } as usize;
        if mapped_bit_position < W {
            permuted_word0 |= 1usize << (WMASK - mapped_bit_position);
        } else {
            permuted_word1 |= 1usize << (2 * W - 1 - mapped_bit_position);
        }
    }
    let mut source_word1 = unsafe { *set1.get_unchecked(1) };
    while source_word1 != 0 {
        let bit_index_in_word1 = (usize::BITS - 1 - source_word1.trailing_zeros()) as usize;
        source_word1 &= source_word1 - 1;
        let source_bit_position = W + bit_index_in_word1;
        // SAFETY: padding bits are zero, so source_bit_position < n <= 2 * W.
        let mapped_bit_position = unsafe { *permutation.get_unchecked(source_bit_position) } as usize;
        if mapped_bit_position < W {
            permuted_word0 |= 1usize << (WMASK - mapped_bit_position);
        } else {
            permuted_word1 |= 1usize << (2 * W - 1 - mapped_bit_position);
        }
    }
    // SAFETY: caller guarantees set2.len() >= 2.
    unsafe {
        *set2.get_unchecked_mut(0) = permuted_word0;
        *set2.get_unchecked_mut(1) = permuted_word1;
    }
}

/// permute_set for words_per_vertex > 2: same algorithm as before, but bit
/// clearing and output writes go through direct unchecked mask ops instead of
/// the generic `BitSetOpsMSB0` trait (bounds-checked `get_mut` per output bit).
#[inline(always)]
fn permute_set_multiword(
    set1: &[usize],
    set2: &mut [usize],
    words_per_vertex: u32,
    permutation: &[u32],
) {
    let word_count = words_per_vertex as usize;
    // SAFETY: caller guarantees set2.len() >= word_count.
    unsafe {
        for word in set2.get_unchecked_mut(..word_count) {
            *word = 0;
        }
    }
    for source_word_index in 0..word_count {
        // SAFETY: caller guarantees set1.len() >= word_count.
        let mut current_set_word = unsafe { *set1.get_unchecked(source_word_index) };
        while current_set_word != 0 {
            let bit_index = (usize::BITS - 1 - current_set_word.trailing_zeros()) as usize;
            current_set_word &= current_set_word - 1;
            // SAFETY: zero padding keeps the source vertex below permutation.len().
            let target_bit_position =
                unsafe { *permutation.get_unchecked((source_word_index << WSHIFT) + bit_index) } as usize;
            let target_word_index = target_bit_position >> WSHIFT;
            let target_bit_in_word = target_bit_position & WMASK;
            // SAFETY: target_word_index < word_count since target_bit_position < n <= word_count*64.
            unsafe {
                *set2.get_unchecked_mut(target_word_index) |= 1usize << (WMASK - target_bit_in_word);
            }
        }
    }
}

pub fn int_to_string(integer: u32, s: &mut String) {
    let int_str = integer.to_string();
    s.push_str(&int_str);
}

/// Computes the orbit decomposition (connected components) of a permutation defined by `permutation`.
///
/// `orbits` is an output buffer. After execution, `orbits[i]` will contain the
/// canonical representative (the smallest index) of the cycle that `i` belongs to.
///
/// `permutation` must be a valid permutation of `0..vertex_count`.
///
/// Returns the total number of disjoint cycles (orbits).
///
/// # Examples
///
/// ```text
/// use canonaut::support::compute_orbits;
/// let permutation = [1, 2, 0, 4, 3];
/// // Cycles: (0 -> 1 -> 2 -> 0) and (3 -> 4 -> 3)
/// // Representatives: min(0,1,2)=0 and min(3,4)=3
///
/// let mut orbits = [0; 5];
/// let count = compute_orbits(&mut orbits, &permutation, 5);
///
/// assert_eq!(count, 2);
/// assert_eq!(orbits, [0, 0, 0, 3, 3]);
/// ```
pub fn compute_orbits(orbits: &mut [u32], permutation: &[u32], vertex_count: u32) -> u32 {
    let limit = vertex_count as usize;
    // Slice inputs to 'vertex_count'. This helps the optimizer elide bounds checks
    // inside the loop because it knows the range of 'vertex_index' matches the slice length.
    let orbits = &mut orbits[..limit];
    let permutation = &permutation[..limit];

    // Initialize with a sentinel (!0 is u32::MAX) to mark unvisited nodes.
    orbits.fill(!0);

    let mut orbit_count = 0;

    for vertex_index in 0..limit {
        // If the node is still marked with the sentinel, it hasn't been visited
        // by a smaller root yet. Thus, 'vertex_index' is the representative of a new cycle.
        if orbits[vertex_index] == !0 {
            orbit_count += 1;
            let root = vertex_index as u32;
            let mut cycle_cursor = vertex_index;

            // Traverse the cycle and mark all nodes with the root ID.
            // We stop when we hit a node that has already been visited (which will be 'vertex_index').
            while orbits[cycle_cursor] == !0 {
                orbits[cycle_cursor] = root;
                // Update cycle_cursor. Note: This will panic if permutation[cycle_cursor] >= vertex_count, ensuring safety.
                cycle_cursor = permutation[cycle_cursor] as usize;
            }
        }
    }

    orbit_count
}

/// Equivalent of C's fmperm(): compute (fixed_points, min_cycle_representatives) from a permutation.
/// fixed_points = vertices fixed by permutation; min_cycle_representatives = fixed_points ∪ {minimum element of each non-trivial cycle}.
#[inline]
pub fn fmperm_safe(
    permutation: &[u32],
    fixed_points: &mut [usize],
    min_cycle_representatives: &mut [usize],
    words_per_vertex: u32,
    vertex_count: u32,
) {
    let word_count = words_per_vertex as usize;
    fixed_points[..word_count].fill(0);
    min_cycle_representatives[..word_count].fill(0);

    if vertex_count as usize <= W {
        // Avoid heap allocation for the common word_count=1 case
        let mut visited: u64 = 0;
        for vertex_index in 0..vertex_count as usize {
            let mapped_vertex = unsafe { *permutation.get_unchecked(vertex_index) } as usize;
            if mapped_vertex == vertex_index {
                unsafe { *fixed_points.get_unchecked_mut(0) |= BIT_USIZE[vertex_index & WMASK] };
                unsafe { *min_cycle_representatives.get_unchecked_mut(0) |= BIT_USIZE[vertex_index & WMASK] };
            } else if (visited >> vertex_index) & 1 == 0 {
                let mut cycle_cursor = vertex_index;
                loop {
                    visited |= 1u64 << cycle_cursor;
                    cycle_cursor = unsafe { *permutation.get_unchecked(cycle_cursor) } as usize;
                    if cycle_cursor == vertex_index {
                        break;
                    }
                }
                unsafe {
                    *min_cycle_representatives.get_unchecked_mut(0) |= BIT_USIZE[vertex_index & WMASK]
                };
            }
        }
    } else {
        // n > 64 needs a visited array.  It is kept in a thread-local scratch
        // buffer that grows monotonically, so repeated canonicalizations do not
        // allocate: this routine runs once per automorphism found, which for
        // highly symmetric graphs is many times per `canonize` call.
        CYCLE_VISITED_SCRATCH.with_borrow_mut(|visited| {
            if visited.len() < vertex_count as usize {
                visited.resize(vertex_count as usize, false);
            }
            visited[..vertex_count as usize].fill(false);

            for vertex_index in 0..vertex_count as usize {
                let mapped_vertex = permutation[vertex_index] as usize;
                if mapped_vertex == vertex_index {
                    fixed_points[vertex_index >> WSHIFT] |= BIT_USIZE[vertex_index & WMASK];
                    min_cycle_representatives[vertex_index >> WSHIFT] |= BIT_USIZE[vertex_index & WMASK];
                } else if !visited[vertex_index] {
                    let mut cycle_cursor = vertex_index;
                    loop {
                        visited[cycle_cursor] = true;
                        cycle_cursor = permutation[cycle_cursor] as usize;
                        if cycle_cursor == vertex_index {
                            break;
                        }
                    }
                    min_cycle_representatives[vertex_index >> WSHIFT] |= BIT_USIZE[vertex_index & WMASK];
                }
            }
        });
    }
}

thread_local! {
    /// Visited-marker scratch for [`fmperm_safe`]'s `n > 64` path.
    static CYCLE_VISITED_SCRATCH: std::cell::RefCell<Vec<bool>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

/// Compute (fixed_points, min_cycle_representatives) from partition at the given level.
/// SAFETY: `lab[vertex_index] < vertex_count ≤ word_count*64`, partition has a sentinel ≤ level before index `vertex_count`,
///         so inner loops terminate before vertex_index=vertex_count; fixed_points/min_cycle_representatives have word_count words.
#[inline]
pub fn mark_partition_fixed_points(
    lab: &[u32],
    partition: &[u32],
    level: u32,
    fixed_points: &mut [usize],
    min_cycle_representatives: &mut [usize],
    words_per_vertex: u32,
    vertex_count: u32,
) {
    let mut vertex_index: u32 = 0;
    let mut min_label_in_cell: u32;
    fixed_points[..words_per_vertex as usize].fill(0);
    min_cycle_representatives[..words_per_vertex as usize].fill(0);
    while vertex_index < vertex_count {
        // SAFETY: vertex_index < vertex_count; lab[vertex_index] < vertex_count ≤ 64 for word_count=1; fixed_points/min_cycle_representatives have word_count=1 word.
        if unsafe { *partition.get_unchecked(vertex_index as usize) } <= level {
            let label_value = unsafe { *lab.get_unchecked(vertex_index as usize) };
            let word_index = (label_value >> WSHIFT) as usize;
            let bit_mask = BIT_USIZE[label_value as usize & WMASK];
            unsafe { *fixed_points.get_unchecked_mut(word_index) |= bit_mask };
            unsafe { *min_cycle_representatives.get_unchecked_mut(word_index) |= bit_mask };
        } else {
            min_label_in_cell = unsafe { *lab.get_unchecked(vertex_index as usize) };
            loop {
                vertex_index += 1;
                let label_value = unsafe { *lab.get_unchecked(vertex_index as usize) };
                if label_value < min_label_in_cell {
                    min_label_in_cell = label_value;
                }
                if unsafe { *partition.get_unchecked(vertex_index as usize) } <= level {
                    break;
                }
            }
            let word_index = (min_label_in_cell >> WSHIFT) as usize;
            unsafe {
                *min_cycle_representatives.get_unchecked_mut(word_index) |=
                    BIT_USIZE[min_label_in_cell as usize & WMASK]
            };
        }
        vertex_index += 1;
    }
}

pub fn sortparallel<T: Ord + Clone, U: Clone>(keys: &mut [T], values: &mut [U]) {
    assert_eq!(keys.len(), values.len());

    // 1. Combine keys and values into a new vector of OWNED tuples (requires cloning)
    let mut pairs: Vec<(T, U)> = keys.iter().cloned().zip(values.iter().cloned()).collect();

    // 2. Sort the combined vector
    pairs.sort_by(|left, right| left.0.cmp(&right.0));

    // 3. Write the sorted data back into the original slices
    for (pair_index, (key, value)) in pairs.into_iter().enumerate() {
        keys[pair_index] = key;
        values[pair_index] = value;
    }
}

#[inline]
pub fn apply_refinement(
    graph: &[usize],
    lab: &mut [u32],
    partition: &mut [u32],
    level: i32,
    numcells: &mut u32,
    qinvar: &mut u32,
    invar: &mut [u32],
    active: &mut [usize],
    code: &mut u32,
    refinement_procedure: Option<
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
    invariant_procedure: Option<
        fn(&[usize], &mut [u32], &mut [u32], i32, u32, u32, &mut [u32], u32, bool, u32, u32) -> (),
    >,
    minimum_invariant_level: i32,
    maximum_invariant_level: i32,
    invariant_argument: u32,
    is_directed_graph: bool,
    words_per_vertex: u32,
    number_of_vertices: u32,
    graph_context: &mut GraphContext,
) {
    let mut cell_work_permutation_value: u32;
    let mut split_index: u32;
    let mut cell_start: u32;
    let mut cell_end: u32;
    let previous_numcells: u32;
    
    
    let mut code_accumulator: i64;
    let mut cells_identical: bool;
    debug_assert!(graph_context.work_permutation.len() >= number_of_vertices as usize);
    refinement_procedure.expect("non-null function pointer")(
        graph,
        lab,
        partition,
        level,
        numcells,
        invar,
        active,
        code,
        words_per_vertex,
        number_of_vertices,
        graph_context,
    );
    let min_invariant_level: i32 = if minimum_invariant_level < 0 {
        -minimum_invariant_level
    } else {
        minimum_invariant_level
    };
    let max_invariant_level: i32 = if maximum_invariant_level < 0 {
        -maximum_invariant_level
    } else {
        maximum_invariant_level
    };
    let invariant_applies = *numcells < number_of_vertices
        && level >= min_invariant_level
        && level <= max_invariant_level;
    if let Some(procedure) = invariant_procedure.filter(|_| invariant_applies) {
        procedure(
            graph,
            lab,
            partition,
            level,
            *numcells,
            *qinvar,
            invar,
            invariant_argument,
            is_directed_graph,
            words_per_vertex,
            number_of_vertices,
        );
        active.clear_words(words_per_vertex as usize);
        for vertex_index in 0..number_of_vertices {
            graph_context.work_permutation[vertex_index as usize] =
                invar[lab[vertex_index as usize] as usize];
        }
        previous_numcells = *numcells;
        cell_start = 0;
        while cell_start < number_of_vertices {
            cell_work_permutation_value = graph_context.work_permutation[cell_start as usize];
            cells_identical = true;
            cell_end = cell_start;
            while partition[cell_end as usize] as i32 > level {
                if graph_context.work_permutation[(cell_end + 1) as usize] != cell_work_permutation_value
                {
                    cells_identical = false;
                }
                cell_end += 1;
            }
            if !cells_identical {
                let cell_slice_end = cell_end as usize + 1;
                sortparallel(
                    &mut graph_context.work_permutation[cell_start as usize..cell_slice_end],
                    &mut lab[cell_start as usize..cell_slice_end],
                );

                split_index = cell_start + 1;
                while split_index <= cell_end {
                    if graph_context.work_permutation[split_index as usize]
                        != graph_context.work_permutation[(split_index - 1) as usize]
                    {
                        assert!(level >= 0);
                        partition[(split_index - 1) as usize] = level as u32;
                        *numcells += 1;
                        active.set_bit(split_index as usize);
                    }
                    split_index += 1;
                }
            }
            cell_start = cell_end + 1;
        }
        if *numcells > previous_numcells {
            *qinvar = 2;
            code_accumulator = *code as i64;
            code_accumulator =
                ((code_accumulator ^ 0o65435_u32 as i64) + *code as i64) & 0o77777_u32 as i64;
            *code = (code_accumulator % 0o77777_u32 as i64) as u32;
        } else {
            *qinvar = 1;
        }
    } else {
        *qinvar = 0;
    };
}

#[inline]
pub fn select_target_cell(
    graph: &[usize],
    lab: &mut [u32],
    partition: &mut [u32],
    level: i32,
    target_cell: &mut [usize],
    target_cellsize: &mut u32,
    cellpos: &mut i32,
    target_cell_level: i32,
    digraph: bool,
    hint: i32,
    target_cell_fn: Option<
        fn(&[usize], &[u32], &[u32], i32, i32, bool, i32, u32, u32, &mut GraphContext) -> i32,
    >,
    words_per_vertex: u32,
    vertex_count: u32,
    graph_context: &mut GraphContext,
) {
    let mut fill_index: i32;
    let cell_start = target_cell_fn.expect("non-null function pointer")(
        graph,
        lab,
        partition,
        level,
        target_cell_level,
        digraph,
        hint,
        words_per_vertex,
        vertex_count,
        graph_context,
    );
    let mut cell_end = cell_start + 1;
    while partition[cell_end as usize] as i32 > level {
        cell_end += 1;
    }
    *target_cellsize = (cell_end - cell_start + 1) as u32;
    target_cell.clear_words(words_per_vertex as usize);
    fill_index = cell_start;
    while fill_index <= cell_end {
        // SAFETY: lab[fill_index] < vertex_count ≤ word_count*BITS, target_cell has word_count words
        let mapped_vertex = unsafe { *lab.get_unchecked(fill_index as usize) } as usize;
        unsafe {
            *target_cell.get_unchecked_mut(mapped_vertex / usize::BITS as usize) |=
                1usize << (usize::BITS - 1 - (mapped_vertex % usize::BITS as usize) as u32);
        }
        fill_index += 1;
    }
    *cellpos = cell_start;
}

pub fn intersect_prune(set1: &mut [usize], set2: &[usize]) {
    intersect_inplace(set1, set2, set1.len());
}

#[inline]
pub fn individualize_vertex(
    lab: &mut [u32],
    partition: &mut [u32],
    level: i32,
    target_cell_index: i32,
    target_vertex: u32,
    active: &mut [usize],
    words_per_vertex: u32,
) {
    active.clear_words(words_per_vertex as usize);
    active.set_bit(target_cell_index as usize);
    let mut shift_index = target_cell_index;
    let mut label_to_write = target_vertex;
    loop {
        // SAFETY: shift_index stays within the cell [target_cell_index, cell_end]; lab[shift_index] < n.
        let next_label = unsafe { *lab.get_unchecked(shift_index as usize) };
        unsafe {
            *lab.get_unchecked_mut(shift_index as usize) = label_to_write;
        }
        shift_index += 1;
        label_to_write = next_label;
        if label_to_write == target_vertex {
            break;
        }
    }
    unsafe {
        *partition.get_unchecked_mut(target_cell_index as usize) = level as u32;
    }
}

/// Port of nauty's `longprune`. Intersects `target_cell` with the second set of every
/// `(fixed_points_i, min_cycle_representatives_i)` pair in `generator_pairs` for which `fixed_points` is a subset of `fixed_points_i`.
#[inline]
pub fn prune_by_generators(
    target_cell: &mut [usize],
    fixed_points: &[usize],
    generator_pairs: &[usize],
    words_per_vertex: u32,
) {
    let word_count = words_per_vertex as usize;
    if word_count == 1 {
        // Fast path: single-word ops — avoid SIMD dispatch overhead.
        let fixed_points_word = unsafe { *fixed_points.get_unchecked(0) };
        let target_cell_word = unsafe { target_cell.get_unchecked_mut(0) };
        for pair in generator_pairs.chunks_exact(2) {
            let pair_fixed_points_word = unsafe { *pair.get_unchecked(0) };
            let pair_min_cycle_representatives_word = unsafe { *pair.get_unchecked(1) };
            if fixed_points_word & !pair_fixed_points_word == 0 {
                *target_cell_word &= pair_min_cycle_representatives_word;
            }
        }
        return;
    }
    for pair in generator_pairs.chunks_exact(2 * word_count) {
        let (pair_fixed_points, pair_min_cycle_representatives) = pair.split_at(word_count);
        if !andn_any(fixed_points, pair_fixed_points, word_count) {
            intersect_inplace(target_cell, pair_min_cycle_representatives, word_count);
        }
    }
}
fn write_formatted(output: &mut OutputStream, format: &'static str, args: &[&dyn std::fmt::Display]) {
    let formatted = match format {
        "{:.0}" => format!("{:.0}", args[0]),
        "{:.0}*{}" => format!("{:.0}*{}", args[0], args[1]),
        _ => panic!("Unsupported format string"),
    };
    let _ = output.write(formatted.as_bytes());
}

pub fn writegroupsize(output: &mut OutputStream, group_size_mantissa: f64, group_size_exponent: u32) {
    if group_size_exponent == 0 {
        write_formatted(output, "{:.0}", &[&(group_size_mantissa + 0.1f64)]);
    } else {
        write_formatted(
            output,
            "{:.0}*{}",
            &[&(group_size_mantissa + 0.1f64), &group_size_exponent],
        );
    }
}

/// Guards against a caller built against a different word size or a different
/// nauty version than this port targets (nauty's `nauty_check`).
pub fn check_build_consistency(wordsize: u32, version: u32) {
    if wordsize != WORDSIZE {
        eprintln!("canonaut: word-size mismatch (built for {WORDSIZE}, called with {wordsize})");
    }
    if version != 28090 {
        eprintln!("canonaut: version mismatch (this port targets nauty 2.8.9)");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The original implementation (for reference/comparison)
    fn orbjoin_original(orbits: &mut [u32], permutation: &[u32], vertex_count: u32) -> u32 {
        let mut vertex_index: u32;
        let mut cycle_cursor: u32;
        let mut orbit_count: u32 = 0;
        vertex_index = 0;
        while vertex_index < vertex_count {
            orbits[vertex_index as usize] = vertex_index;
            vertex_index += 1;
        }
        vertex_index = 0;
        while vertex_index < vertex_count {
            if orbits[vertex_index as usize] == vertex_index {
                orbit_count += 1;
                cycle_cursor = vertex_index;
                while {
                    cycle_cursor = permutation[cycle_cursor as usize];
                    cycle_cursor != vertex_index
                } {
                    orbits[cycle_cursor as usize] = vertex_index;
                }
            }
            vertex_index += 1;
        }
        orbit_count
    }

    // Helper: Fisher-Yates shuffle to generate a random permutation
    // Uses a simple LCG to avoid adding 'rand' dependency for this snippet.
    fn generate_random_permutation(vertex_count: usize, seed: u64) -> Vec<u32> {
        let mut rng_state = seed;
        let mut permutation: Vec<u32> = (0..vertex_count as u32).collect();

        for shuffle_index in (1..vertex_count).rev() {
            // Linear Congruential Generator
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let swap_index = (rng_state as usize) % (shuffle_index + 1);
            permutation.swap(shuffle_index, swap_index);
        }
        permutation
    }

    #[test]
    fn test_equivalence_fuzz() {
        let vertex_count = 100;
        let iterations = 1000;

        for seed in 0..iterations {
            let permutation = generate_random_permutation(vertex_count, seed as u64);

            // Run Original
            let mut orbits_orig = vec![0; vertex_count];
            let count_orig = orbjoin_original(&mut orbits_orig, &permutation, vertex_count as u32);

            // Run Optimized
            let mut orbits_opt = vec![0; vertex_count];
            let count_opt = compute_orbits(&mut orbits_opt, &permutation, vertex_count as u32);

            // Assert Correctness
            assert_eq!(
                count_orig, count_opt,
                "Cycle counts do not match for seed {}",
                seed
            );
            assert_eq!(
                orbits_orig, orbits_opt,
                "Orbit arrays do not match for seed {}",
                seed
            );
        }
    }

    #[test]
    fn test_edge_cases() {
        // Case 1: Identity (everyone points to self)
        let permutation = [0, 1, 2, 3];
        let mut orbits = [0; 4];
        assert_eq!(compute_orbits(&mut orbits, &permutation, 4), 4);
        assert_eq!(orbits, [0, 1, 2, 3]);

        // Case 2: One big cycle (0->1->2->3->0)
        let permutation = [1, 2, 3, 0];
        let mut orbits = [0; 4];
        assert_eq!(compute_orbits(&mut orbits, &permutation, 4), 1);
        assert_eq!(orbits, [0, 0, 0, 0]);

        // Case 3: Swap pairs (0<->1, 2<->3)
        let permutation = [1, 0, 3, 2];
        let mut orbits = [0; 4];
        assert_eq!(compute_orbits(&mut orbits, &permutation, 4), 2);
        assert_eq!(orbits, [0, 0, 2, 2]);
    }
}
