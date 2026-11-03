use crate::bit_ops::{W, WMASK, WSHIFT};
use crate::support::permute_set;
use crate::search::canonicalize;
use crate::simd_ops::{and_popcount, or_reduce_masked};
use crate::structs::manager::*;
use crate::structs::BitSetOpsMSB0;
use crate::utilities::{cleanup, mash, words_needed};

pub static DISPATCH_GRAPH: DispatchVec = {
    DispatchVec {
        is_automorphism: Some(is_automorphism),
        test_canonical_labeling: Some(test_canonical_labeling),
        update_canonical_graph: Some(updatecan),
        refine: Some(refine),
        refine1: Some(refine1),
        cheap_automorphism_check: Some(cheap_automorphism_check),
        target_cell: Some(target_cell),
        freedyn: Some(refine_freedyn),
        init: None,
        cleanup: None,
    }
};


pub struct GraphContext {
    pub work_bitset: Vec<usize>,
    pub work_permutation: Vec<u32>,
    pub bucket: Vec<u32>,
    pub dense_workspace: Vec<usize>,
}

impl GraphContext {
    pub fn new(number_of_vertices: usize) -> Self {
        let words_per_vertex = words_needed(number_of_vertices);
        GraphContext {
            work_bitset: vec![0; words_per_vertex],
            work_permutation: vec![0; number_of_vertices],
            bucket: vec![0; number_of_vertices + 2],  // +2 slack required by refine1
            dense_workspace: vec![0; words_per_vertex],
        }
    }
    pub fn ensure_workperm_size(&mut self, number_of_vertices: u32) {
        let vertex_count = number_of_vertices as usize;
        if self.work_permutation.len() < vertex_count {
            self.work_permutation.resize(vertex_count, 0);
        }
    }

    pub fn ensure_bucket_size(&mut self, number_of_vertices: u32) {
        let bucket_len = number_of_vertices as usize + 2;
        if self.bucket.len() < bucket_len {
            self.bucket.resize(bucket_len, 0);
        }
    }

    pub fn ensure_workset_size(&mut self, words_per_vertex: u32) {
        let words_per_vertex = words_per_vertex as usize;
        if self.work_bitset.len() < words_per_vertex {
            self.work_bitset.resize(words_per_vertex, 0);
        }
    }

    fn ensure_dnwork_size(&mut self, words_per_vertex: u32) {
        let words_per_vertex = words_per_vertex as usize;
        if self.dense_workspace.len() < 2 * 500 * words_per_vertex {
            self.dense_workspace.resize(2 * 500 * words_per_vertex, 0);
        }
    }
}

pub fn is_automorphism(
    graph: &[usize],
    permutation: &[u32],
    is_directed_graph: bool,
    words_per_vertex: u32,
    number_of_vertices: u32,
) -> bool {
    let words_per_vertex = words_per_vertex as usize;
    let number_of_vertices = number_of_vertices as usize;

    match words_per_vertex {
        1 => is_automorphism_single_word(graph, permutation, is_directed_graph, number_of_vertices),
        2 => is_automorphism_two_words(graph, permutation, is_directed_graph, number_of_vertices),
        _ => is_automorphism_multiword(graph, permutation, is_directed_graph, words_per_vertex, number_of_vertices),
    }
}

/// Automorphism check for the single-word case (n ≤ 64, words_per_vertex == 1).
/// Iterates set bits by clearing the lowest bit; no slice-indexing overhead.
#[inline(always)]
fn is_automorphism_single_word(
    graph: &[usize],
    permutation: &[u32],
    is_directed_graph: bool,
    number_of_vertices: usize,
) -> bool {
    for vertex in 0..number_of_vertices {
        // SAFETY: vertex < number_of_vertices = graph.len() = permutation.len()
        let adjacency_row = unsafe { *graph.get_unchecked(vertex) };
        let mapped_vertex = unsafe { *permutation.get_unchecked(vertex) } as usize;
        // SAFETY: mapped_vertex < n = graph.len() (permutation is a valid permutation)
        let mapped_adjacency_row = unsafe { *graph.get_unchecked(mapped_vertex) };

        // Undirected: check only j > vertex (upper triangle); directed: check all.
        let mut remaining_neighbors = if is_directed_graph {
            adjacency_row
        } else {
            // Clear the top (vertex+1) bits — those are columns j ≤ vertex in MSB-0 encoding.
            if vertex + 1 >= usize::BITS as usize {
                0
            } else {
                adjacency_row & (usize::MAX >> (vertex + 1))
            }
        };

        while remaining_neighbors != 0 {
            let neighbor = (usize::BITS - 1 - remaining_neighbors.trailing_zeros()) as usize;
            remaining_neighbors &= remaining_neighbors - 1;

            // SAFETY: neighbor < n = permutation.len()
            let mapped_neighbor = unsafe { *permutation.get_unchecked(neighbor) } as usize;
            if mapped_adjacency_row & (1usize << (usize::BITS as usize - 1 - mapped_neighbor)) == 0 {
                return false;
            }
        }
    }
    true
}

/// Automorphism check for words_per_vertex == 2 (64 < n ≤ 128).
/// Loads both adjacency words into registers per vertex — no Option overhead,
/// no slice bounds checks in the per-edge inner loop.
#[inline(always)]
fn is_automorphism_two_words(
    graph: &[usize],
    permutation: &[u32],
    is_directed_graph: bool,
    number_of_vertices: usize,
) -> bool {
    for vertex in 0..number_of_vertices {
        // SAFETY: vertex < number_of_vertices; graph has number_of_vertices*2 words; permutation has number_of_vertices entries.
        let row_base = vertex * 2;
        let row_word0 = unsafe { *graph.get_unchecked(row_base) };
        let row_word1 = unsafe { *graph.get_unchecked(row_base + 1) };

        let mapped_vertex = unsafe { *permutation.get_unchecked(vertex) } as usize;
        let mapped_base = mapped_vertex * 2;
        let mapped_word0 = unsafe { *graph.get_unchecked(mapped_base) };
        let mapped_word1 = unsafe { *graph.get_unchecked(mapped_base + 1) };

        // Undirected: check only j > vertex (upper triangle).
        // In MSB-0 encoding, bit j in word 0 is column j (0-63); word 1 is columns 64-127.
        let (mut remaining_word0, mut remaining_word1) = if !is_directed_graph {
            if vertex < W {
                (row_word0 & if vertex + 1 < W { usize::MAX >> (vertex + 1) } else { 0 }, row_word1)
            } else {
                let bit_in_word1 = vertex - W;
                (0, row_word1 & if bit_in_word1 + 1 < W { usize::MAX >> (bit_in_word1 + 1) } else { 0 })
            }
        } else {
            (row_word0, row_word1)
        };

        // Iterate bits in word 0 (columns 0..63)
        while remaining_word0 != 0 {
            let column = (usize::BITS - 1 - remaining_word0.trailing_zeros()) as usize;
            remaining_word0 &= remaining_word0 - 1;
            // SAFETY: column < 64 < n = permutation.len()
            let mapped_column = unsafe { *permutation.get_unchecked(column) } as usize;
            // Test bit mapped_column in mapped vertex's row (MSB-0): word mapped_column>>6, bit 63-(mapped_column&63)
            let (selected_word, selected_bit) = if mapped_column < W {
                (mapped_word0, mapped_column)
            } else {
                (mapped_word1, mapped_column - W)
            };
            if selected_word >> (usize::BITS as usize - 1 - selected_bit) & 1 == 0 {
                return false;
            }
        }

        // Iterate bits in word 1 (columns 64..127)
        while remaining_word1 != 0 {
            let bit_in_word1 = (usize::BITS - 1 - remaining_word1.trailing_zeros()) as usize;
            let column = W + bit_in_word1;
            remaining_word1 &= remaining_word1 - 1;
            // SAFETY: zero padding keeps column < n <= 2 * W.
            let mapped_column = unsafe { *permutation.get_unchecked(column) } as usize;
            let (selected_word, selected_bit) = if mapped_column < W {
                (mapped_word0, mapped_column)
            } else {
                (mapped_word1, mapped_column - W)
            };
            if selected_word >> (usize::BITS as usize - 1 - selected_bit) & 1 == 0 {
                return false;
            }
        }
    }
    true
}

/// Automorphism check for words_per_vertex > 2.
fn is_automorphism_multiword(
    graph: &[usize],
    permutation: &[u32],
    is_directed_graph: bool,
    words_per_vertex: usize,
    number_of_vertices: usize,
) -> bool {
    for vertex in 0..number_of_vertices {
        let current_edges = unsafe { graph.get_unchecked(vertex * words_per_vertex..(vertex + 1) * words_per_vertex) };
        let mapped_vertex = unsafe { *permutation.get_unchecked(vertex) } as usize;
        let mapped_edges = unsafe { graph.get_unchecked(mapped_vertex * words_per_vertex..(mapped_vertex + 1) * words_per_vertex) };

        let first_neighbor = if is_directed_graph { 0 } else { vertex + 1 };
        let first_word = first_neighbor >> WSHIFT;
        for (word_index, &word) in current_edges.iter().enumerate().skip(first_word) {
            let mut remaining = if word_index == first_word {
                word & (usize::MAX >> (first_neighbor & WMASK))
            } else { word };
            while remaining != 0 {
                let neighbor = word_index * W + WMASK - remaining.trailing_zeros() as usize;
                remaining &= remaining - 1;
                // SAFETY: input padding is zero and the permutation maps 0..n onto itself.
                let mapped_neighbor = unsafe { *permutation.get_unchecked(neighbor) } as usize;
                if !mapped_edges.is_bit_set(mapped_neighbor) { return false; }
            }
        }
    }
    true
}

pub fn test_canonical_labeling(
    graph: &[usize],
    canonical_graph: &[usize],
    labeling: &mut [u32],
    matching_rows: &mut u32,
    words_per_vertex: u32,
    number_of_vertices: u32,
    context: &mut GraphContext,
) -> i32 {
    let mut vertex_index: u32;
    let mut word_index: u32;
    let mut canon_cursor: &[usize];
    debug_assert!(context.work_permutation.len() >= number_of_vertices as usize);
    debug_assert!(context.work_bitset.len() >= words_per_vertex as usize);
    // Build the inverse permutation mapping
    for (position, &label) in labeling[..(number_of_vertices) as usize].iter().enumerate() {
        context.work_permutation[label as usize] = position as u32;
    }

    vertex_index = 0;
    canon_cursor = canonical_graph;
    while vertex_index < number_of_vertices {
        permute_set(
            &graph[(words_per_vertex as usize).wrapping_mul(labeling[vertex_index as usize] as usize)..],
            &mut context.work_bitset,
            words_per_vertex,
            &context.work_permutation,
        );
        word_index = 0;
        while word_index < words_per_vertex {
            if context.work_bitset[word_index as usize] < canon_cursor[word_index as usize] {
                *matching_rows = vertex_index;
                return -(1);
            } else if context.work_bitset[word_index as usize] > canon_cursor[word_index as usize] {
                *matching_rows = vertex_index;
                return 1;
            }
            word_index += 1;
        }
        vertex_index += 1;
        canon_cursor = &canon_cursor[words_per_vertex as usize..];
    }
    *matching_rows = number_of_vertices;
    0
}

pub fn updatecan(
    graph: &[usize],
    canonical_graph: &mut [usize],
    labeling: &[u32],
    matching_rows: u32,
    words_per_vertex: u32,
    number_of_vertices: u32,
    context: &mut GraphContext,
) {
    let mut vertex_index: u32;
    let mut canon_cursor: &mut [usize];
    debug_assert!(context.work_permutation.len() >= number_of_vertices as usize);

    vertex_index = 0;
    while vertex_index < number_of_vertices {
        context.work_permutation[labeling[vertex_index as usize] as usize] = vertex_index;
        vertex_index += 1;
    }
    vertex_index = matching_rows;
    canon_cursor = &mut canonical_graph[(words_per_vertex as usize).wrapping_mul(matching_rows as usize)..];
    while vertex_index < number_of_vertices {
        permute_set(
            &graph[(words_per_vertex as usize).wrapping_mul(labeling[vertex_index as usize] as usize)..],
            canon_cursor,
            words_per_vertex,
            &context.work_permutation,
        );
        vertex_index += 1;
        canon_cursor = &mut canon_cursor[words_per_vertex as usize..];
    }
}

pub fn refine(
    graph: &[usize],
    labeling: &mut [u32],
    partition: &mut [u32],
    level: i32,
    number_of_cells: &mut u32,
    neighbor_counts: &mut [u32],
    active_cells: &mut [usize],
    hash_code: &mut u32,
    words_per_vertex: u32,
    number_of_vertices: u32,
    context: &mut GraphContext,
) {
    // All index variables use i32 (like refine1) — avoids Option<usize> overhead.
    // Sentinel value -1 is used where Option would have been None.
    let mut scan_cursor: i32;
    let mut cell_split_low: i32;
    let mut cell_split_high: i32;
    let mut candidate_label: u32;
    let mut neighbor_row: &[usize];
    let mut pivot_cell_start: i32;
    let mut pivot_cell_end: i32;
    let mut scan_cell_start: i32;
    let mut scan_cell_end: i32;
    let mut popcount_value: u32;
    let mut min_count: u32;
    let mut max_count: u32;
    let mut hash_accumulator: i64;
    let mut pivot_row: &[usize];
    let mut largest_bucket_size: i32;
    let mut largest_bucket_start: i32 = 0;
    let mut hint: i32;
    let word_count = words_per_vertex as usize;
    let vertex_count = number_of_vertices as i32;

    // Safety invariants throughout:
    //   labeling.len() >= vertex_count, partition.len() >= vertex_count, neighbor_counts.len() >= vertex_count
    //   all labeling[i] are valid vertex indices 0..vertex_count
    //   bucket.len() >= vertex_count+2, work_permutation.len() >= vertex_count, work_bitset.len() >= word_count
    //   partition has a cell-end marker (value <= level) at some index < vertex_count,
    //   so partition-scanning loops terminate before reaching vertex_count.
    debug_assert!(context.work_permutation.len() >= number_of_vertices as usize);
    debug_assert!(context.work_bitset.len() >= words_per_vertex as usize);
    debug_assert!(context.bucket.len() >= number_of_vertices as usize + 2);

    hash_accumulator = *number_of_cells as i64;
    hint = 0;
    while *number_of_cells < number_of_vertices && {
        pivot_cell_start = hint;
        active_cells.is_bit_set(pivot_cell_start as usize)
            || {
                pivot_cell_start = active_cells.next_set_bit((pivot_cell_start + 1) as usize)
                    .map(|v| v as i32).unwrap_or(-1);
                pivot_cell_start >= 0
            }
            || {
                pivot_cell_start = active_cells.next_set_bit(0)
                    .map(|v| v as i32).unwrap_or(-1);
                pivot_cell_start >= 0
            }
    } {
        active_cells.unset_bit(pivot_cell_start as usize);

        pivot_cell_end = pivot_cell_start;
        // SAFETY: partition has cell-end marker before vertex_count; loop terminates.
        while unsafe { *partition.get_unchecked(pivot_cell_end as usize) } as i32 > level {
            pivot_cell_end += 1;
        }
        hash_accumulator = mash(hash_accumulator, (pivot_cell_start + pivot_cell_end) as i64);

        if pivot_cell_start == pivot_cell_end {
            // SAFETY: labeling[pivot_cell_start] < vertex_count; row fits in graph slice.
            let row_start = unsafe { *labeling.get_unchecked(pivot_cell_start as usize) } as usize * word_count;
            pivot_row = unsafe { graph.get_unchecked(row_start..row_start + word_count) };

            scan_cell_start = 0;
            while scan_cell_start < vertex_count {
                scan_cell_end = scan_cell_start;
                // SAFETY: cell-end marker guarantees termination before vertex_count.
                while unsafe { *partition.get_unchecked(scan_cell_end as usize) } as i32 > level {
                    scan_cell_end += 1;
                }
                if scan_cell_start != scan_cell_end {
                    cell_split_low = scan_cell_start;
                    cell_split_high = scan_cell_end;
                    while cell_split_low <= cell_split_high {
                        // SAFETY: cell_split_low in [scan_cell_start, scan_cell_end] ⊂ [0, vertex_count).
                        candidate_label = unsafe { *labeling.get_unchecked(cell_split_low as usize) };
                        // Direct MSB-0 bit test — avoids is_bit_set dispatch overhead.
                        // SAFETY: candidate_label < vertex_count ≤ pivot_row.len()*64, so word_idx in bounds.
                        let candidate_vertex = candidate_label as usize;
                        if unsafe { *pivot_row.get_unchecked(candidate_vertex >> WSHIFT) }
                            >> (usize::BITS - 1 - (candidate_vertex & WMASK) as u32)
                            & 1 != 0
                        {
                            cell_split_low += 1;
                        } else {
                            // SAFETY: cell_split_high >= cell_split_low >= scan_cell_start >= 0.
                            unsafe {
                                *labeling.get_unchecked_mut(cell_split_low as usize) =
                                    *labeling.get_unchecked(cell_split_high as usize);
                                *labeling.get_unchecked_mut(cell_split_high as usize) = candidate_label;
                            }
                            if cell_split_high == 0 { cell_split_high = -1; break; }
                            cell_split_high -= 1;
                        }
                    }
                    if cell_split_high >= scan_cell_start && cell_split_low <= scan_cell_end {
                        // SAFETY: cell_split_high in [scan_cell_start, scan_cell_end] ⊂ [0, vertex_count).
                        unsafe { *partition.get_unchecked_mut(cell_split_high as usize) = level as u32 };
                        hash_accumulator = mash(hash_accumulator, cell_split_high as i64);
                        *number_of_cells += 1;
                        if active_cells.is_bit_set(scan_cell_start as usize)
                            || (cell_split_high - scan_cell_start) >= (scan_cell_end - cell_split_low)
                        {
                            active_cells.set_bit(cell_split_low as usize);
                            if cell_split_low == scan_cell_end { hint = cell_split_low; }
                        } else {
                            active_cells.set_bit(scan_cell_start as usize);
                            if cell_split_high == scan_cell_start { hint = scan_cell_start; }
                        }
                    }
                }
                scan_cell_start = scan_cell_end + 1;
            }
        } else {
            context.work_bitset.clear_words(word_count);
            for pivot_member_index in pivot_cell_start as usize..=pivot_cell_end as usize {
                // SAFETY: pivot_member_index in [pivot_cell_start, pivot_cell_end] ⊂ [0, vertex_count); labeling[pivot_member_index] < vertex_count.
                context.work_bitset.set_bit(unsafe { *labeling.get_unchecked(pivot_member_index) } as usize);
            }

            hash_accumulator = mash(hash_accumulator, (pivot_cell_end - pivot_cell_start + 1) as i64);
            scan_cell_start = 0;
            while scan_cell_start < vertex_count {
                scan_cell_end = scan_cell_start;
                // SAFETY: cell-end marker guarantees termination before vertex_count.
                while unsafe { *partition.get_unchecked(scan_cell_end as usize) } as i32 > level {
                    scan_cell_end += 1;
                }
                if scan_cell_start != scan_cell_end {
                    scan_cursor = scan_cell_start;
                    // SAFETY: labeling[scan_cursor] < vertex_count; row fits in graph (vertex_count*word_count words).
                    let row_start = unsafe { *labeling.get_unchecked(scan_cursor as usize) } as usize * word_count;
                    neighbor_row = unsafe { graph.get_unchecked(row_start..row_start + word_count) };
                    popcount_value = and_popcount(&context.work_bitset, neighbor_row, word_count);
                    max_count = popcount_value;
                    min_count = popcount_value;
                    // SAFETY: scan_cursor < vertex_count.
                    unsafe { *neighbor_counts.get_unchecked_mut(scan_cursor as usize) = popcount_value };
                    context.bucket[popcount_value as usize] = 1;
                    loop {
                        scan_cursor += 1;
                        if scan_cursor > scan_cell_end { break; }
                        // SAFETY: scan_cursor <= scan_cell_end < vertex_count; labeling[scan_cursor] < vertex_count.
                        let row_start = unsafe { *labeling.get_unchecked(scan_cursor as usize) } as usize * word_count;
                        neighbor_row = unsafe { graph.get_unchecked(row_start..row_start + word_count) };
                        popcount_value = and_popcount(&context.work_bitset, neighbor_row, word_count);
                        while min_count > popcount_value { min_count -= 1; context.bucket[min_count as usize] = 0; }
                        while max_count < popcount_value { max_count += 1; context.bucket[max_count as usize] = 0; }
                        context.bucket[popcount_value as usize] += 1;
                        // SAFETY: scan_cursor < vertex_count.
                        unsafe { *neighbor_counts.get_unchecked_mut(scan_cursor as usize) = popcount_value };
                    }
                    if min_count == max_count {
                        hash_accumulator = mash(hash_accumulator, (min_count + scan_cell_start as u32) as i64);
                    } else {
                        cell_split_low = scan_cell_start;
                        largest_bucket_size = -1;
                        for popcount_bucket in min_count as usize..=max_count as usize {
                            if context.bucket[popcount_bucket] != 0 {
                                cell_split_high = cell_split_low + context.bucket[popcount_bucket] as i32;
                                context.bucket[popcount_bucket] = cell_split_low as u32;
                                hash_accumulator = mash(hash_accumulator, (popcount_bucket as i32 + cell_split_low) as i64);
                                if cell_split_high - cell_split_low > largest_bucket_size {
                                    largest_bucket_size = cell_split_high - cell_split_low;
                                    largest_bucket_start = cell_split_low;
                                }
                                if cell_split_low != scan_cell_start {
                                    active_cells.set_bit(cell_split_low as usize);
                                    if cell_split_high - cell_split_low == 1 { hint = cell_split_low; }
                                    *number_of_cells += 1;
                                }
                                if cell_split_high <= scan_cell_end {
                                    // SAFETY: cell_split_high >= 1 (cell_split_low != scan_cell_start ≥ 0 → cell_split_low ≥ 1 → cell_split_high > cell_split_low ≥ 1).
                                    unsafe { *partition.get_unchecked_mut((cell_split_high - 1) as usize) = level as u32 };
                                }
                                cell_split_low = cell_split_high;
                            }
                        }

                        let reorder_start = scan_cell_start as usize;
                        let reorder_end   = scan_cell_end as usize;
                        for reorder_index in reorder_start..=reorder_end {
                            // SAFETY: reorder_index < vertex_count; neighbor_counts[reorder_index] ≤ vertex_count < bucket.len().
                            let bucket_index = unsafe { *neighbor_counts.get_unchecked(reorder_index) } as usize;
                            let position = unsafe { *context.bucket.get_unchecked(bucket_index) };
                            unsafe { *context.bucket.get_unchecked_mut(bucket_index) += 1 };
                            unsafe {
                                *context.work_permutation.get_unchecked_mut(position as usize) =
                                    *labeling.get_unchecked(reorder_index);
                            }
                        }
                        for reorder_index in reorder_start..=reorder_end {
                            // SAFETY: reorder_index < vertex_count; work_permutation[reorder_index] is valid.
                            unsafe {
                                *labeling.get_unchecked_mut(reorder_index) =
                                    *context.work_permutation.get_unchecked(reorder_index);
                            }
                        }
                        if !active_cells.is_bit_set(scan_cell_start as usize) {
                            active_cells.set_bit(scan_cell_start as usize);
                            active_cells.unset_bit(largest_bucket_start as usize);
                        }
                    }
                }
                scan_cell_start = scan_cell_end + 1;
            }
        }
    }
    hash_accumulator = mash(hash_accumulator, *number_of_cells as i64);
    *hash_code = cleanup(hash_accumulator) as u32;
}

pub fn refine1(
    graph: &[usize],
    labeling: &mut [u32],
    partition: &mut [u32],
    level: i32,
    number_of_cells: &mut u32,
    neighbor_counts: &mut [u32],
    active_cells: &mut [usize],
    hash_code: &mut u32,
    _words_per_vertex: u32,
    number_of_vertices: u32,
    context: &mut GraphContext,
) {
    let mut scan_cursor: i32;
    let mut cell_split_low: i32;
    let mut cell_split_high: i32;
    let mut candidate_label: u32;
    let mut intersection_word: usize;
    let mut pivot_cell_start: i32;
    let mut pivot_cell_end: i32;
    let mut scan_cell_start: i32;
    let mut scan_cell_end: i32;
    let mut popcount_value: u32;
    let mut min_count: i32;
    let mut max_count: i32;
    let mut hash_accumulator: i64;
    let mut pivot_row_slice: &[usize];
    let mut pivot_bitset_word: usize;
    let mut largest_bucket_size: i32;
    let mut largest_bucket_start: i32 = 0;
    let mut hint: i32;
    let active_cell_word: &mut usize = &mut active_cells[0];

    let vertex_count = number_of_vertices as usize;
    // Safety invariants that hold throughout (established by callers):
    //   graph.len() == vertex_count, labeling.len() >= vertex_count, partition.len() >= vertex_count, neighbor_counts.len() >= vertex_count
    //   all labeling[i] are valid vertex indices in 0..vertex_count
    //   bucket.len() >= vertex_count+2, work_permutation.len() >= vertex_count
    //   The partition always has a cell-end marker (value <= level) before index vertex_count,
    //   so partition-scanning loops terminate before reaching vertex_count.
    debug_assert!(graph.len() >= vertex_count);
    debug_assert!(labeling.len() >= vertex_count);
    debug_assert!(partition.len() >= vertex_count);
    debug_assert!(neighbor_counts.len() >= vertex_count);
    debug_assert!(context.bucket.len() >= vertex_count + 2);
    debug_assert!(context.work_permutation.len() >= vertex_count);

    // Inline bit-ops on single-word active set (word_count=1): no bounds-check branch.
    macro_rules! a_test   { ($x:expr) => { (*active_cell_word >> (usize::BITS - 1 - $x as u32)) & 1 != 0 } }
    macro_rules! a_set    { ($x:expr) => { *active_cell_word |= 1usize << (usize::BITS - 1 - $x as u32) } }
    macro_rules! a_unset  { ($x:expr) => { *active_cell_word &= !(1usize << (usize::BITS - 1 - $x as u32)) } }
    // next bit in active_cell_word at or after `pos+1` (MSB-0); returns -1 if none.
    macro_rules! a_next   { ($pos:expr) => {{
        let masked = *active_cell_word & if ($pos as usize) < WMASK { usize::MAX >> ($pos as u32 + 1) } else { 0 };
        if masked != 0 { masked.leading_zeros() as i32 } else { -1 }
    }};}
    macro_rules! a_first  { () => {{
        if *active_cell_word != 0 { active_cell_word.leading_zeros() as i32 } else { -1 }
    }};}
    hash_accumulator = *number_of_cells as i64;
    hint = 0;
    while *number_of_cells < number_of_vertices && {
        pivot_cell_start = hint;
        a_test!(pivot_cell_start)
            || { pivot_cell_start = a_next!(pivot_cell_start); pivot_cell_start >= 0 }
            || { pivot_cell_start = a_first!();      pivot_cell_start >= 0 }
    } {
        a_unset!(pivot_cell_start);
        pivot_cell_end = pivot_cell_start;
        // SAFETY: partition has a cell-end marker before vertex_count; loop terminates < vertex_count.
        while unsafe { *partition.get_unchecked(pivot_cell_end as usize) } as i32 > level {
            pivot_cell_end += 1;
        }
        hash_accumulator =
            ((hash_accumulator ^ 0o65435_u32 as i64) + (pivot_cell_start + pivot_cell_end) as i64) & 0o77777_u32 as i64;
        if pivot_cell_start == pivot_cell_end {
            // SAFETY: pivot_cell_start < vertex_count, labeling[pivot_cell_start] < vertex_count, graph.len() == vertex_count → pivot_row_slice non-empty.
            pivot_row_slice = unsafe { graph.get_unchecked(*labeling.get_unchecked(pivot_cell_start as usize) as usize..) };
            let pivot_row_word = pivot_row_slice[0];
            scan_cell_start = 0;
            while scan_cell_start < number_of_vertices as i32 {
                scan_cell_end = scan_cell_start;
                // SAFETY: partition cell-end marker guarantees termination before vertex_count.
                while unsafe { *partition.get_unchecked(scan_cell_end as usize) } as i32 > level {
                    scan_cell_end += 1;
                }
                if scan_cell_start != scan_cell_end {
                    cell_split_low = scan_cell_start;
                    cell_split_high = scan_cell_end;
                    while cell_split_low <= cell_split_high {
                        // SAFETY: cell_split_low <= scan_cell_end < vertex_count; labeling[cell_split_low] < vertex_count.
                        candidate_label = unsafe { *labeling.get_unchecked(cell_split_low as usize) };
                        // Direct MSB-0 bit test — avoids trait-method branch (candidate_label < 64).
                        if pivot_row_word >> (usize::BITS - 1 - candidate_label) & 1 != 0 {
                            cell_split_low += 1;
                        } else {
                            // SAFETY: cell_split_high <= scan_cell_end < vertex_count.
                            unsafe {
                                *labeling.get_unchecked_mut(cell_split_low as usize) =
                                    *labeling.get_unchecked(cell_split_high as usize);
                                *labeling.get_unchecked_mut(cell_split_high as usize) = candidate_label;
                            }
                            cell_split_high -= 1;
                        }
                    }
                    if cell_split_high >= scan_cell_start && cell_split_low <= scan_cell_end {
                        // SAFETY: cell_split_high < vertex_count (cell end marker).
                        unsafe { *partition.get_unchecked_mut(cell_split_high as usize) = level as u32 };
                        hash_accumulator =
                            ((hash_accumulator ^ 0o65435_u32 as i64) + cell_split_high as i64) & 0o77777_u32 as i64;
                        *number_of_cells += 1;
                        if a_test!(scan_cell_start) || cell_split_high - scan_cell_start >= scan_cell_end - cell_split_low {
                            a_set!(cell_split_low);
                            if cell_split_low == scan_cell_end {
                                hint = cell_split_low;
                            }
                        } else {
                            a_set!(scan_cell_start);
                            if cell_split_high == scan_cell_start {
                                hint = scan_cell_start;
                            }
                        }
                    }
                }
                scan_cell_start = scan_cell_end + 1;
            }
        } else {
            pivot_bitset_word = 0_usize;
            scan_cursor = pivot_cell_start;
            while scan_cursor <= pivot_cell_end {
                // SAFETY: scan_cursor < vertex_count, labeling[scan_cursor] < vertex_count.
                let pivot_member_label = unsafe { *labeling.get_unchecked(scan_cursor as usize) };
                pivot_bitset_word |= 1usize << (usize::BITS - 1 - pivot_member_label);
                scan_cursor += 1;
            }
            hash_accumulator = ((hash_accumulator ^ 0o65435_u32 as i64) + (pivot_cell_end - pivot_cell_start + 1) as i64)
                & 0o77777_u32 as i64;
            scan_cell_start = 0;
            while scan_cell_start < number_of_vertices as i32 {
                scan_cell_end = scan_cell_start;
                // SAFETY: partition cell-end marker.
                while unsafe { *partition.get_unchecked(scan_cell_end as usize) } as i32 > level {
                    scan_cell_end += 1;
                }
                if scan_cell_start != scan_cell_end {
                    scan_cursor = scan_cell_start;
                    // SAFETY: scan_cursor < vertex_count, labeling[scan_cursor] < vertex_count, graph.len() == vertex_count.
                    intersection_word = pivot_bitset_word & unsafe { *graph.get_unchecked(*labeling.get_unchecked(scan_cursor as usize) as usize) };
                    popcount_value = intersection_word.count_ones();
                    max_count = popcount_value as i32;
                    min_count = max_count;
                    // SAFETY: scan_cursor < vertex_count.
                    unsafe { *neighbor_counts.get_unchecked_mut(scan_cursor as usize) = popcount_value };
                    // SAFETY: popcount_value <= vertex_count < bucket.len().
                    unsafe { *context.bucket.get_unchecked_mut(popcount_value as usize) = 1 };
                    loop {
                        scan_cursor += 1;
                        if scan_cursor > scan_cell_end {
                            break;
                        }
                        // SAFETY: scan_cursor <= scan_cell_end < vertex_count.
                        intersection_word = pivot_bitset_word & unsafe { *graph.get_unchecked(*labeling.get_unchecked(scan_cursor as usize) as usize) };
                        popcount_value = intersection_word.count_ones();
                        while min_count > popcount_value as i32 {
                            min_count -= 1;
                            // SAFETY: min_count in 0..=vertex_count.
                            unsafe { *context.bucket.get_unchecked_mut(min_count as usize) = 0 };
                        }
                        while max_count < popcount_value as i32 {
                            max_count += 1;
                            // SAFETY: max_count in 0..=vertex_count.
                            unsafe { *context.bucket.get_unchecked_mut(max_count as usize) = 0 };
                        }
                        unsafe { *context.bucket.get_unchecked_mut(popcount_value as usize) += 1 };
                        unsafe { *neighbor_counts.get_unchecked_mut(scan_cursor as usize) = popcount_value };
                    }
                    if min_count == max_count {
                        hash_accumulator = ((hash_accumulator ^ 0o65435_u32 as i64) + (min_count + scan_cell_start) as i64)
                            & 0o77777_u32 as i64;
                    } else {
                        cell_split_low = scan_cell_start;
                        largest_bucket_size = -(1);
                        scan_cursor = min_count;
                        while scan_cursor <= max_count {
                            // SAFETY: scan_cursor (a count) <= vertex_count < bucket.len().
                            if unsafe { *context.bucket.get_unchecked(scan_cursor as usize) } != 0 {
                                cell_split_high = cell_split_low + unsafe { *context.bucket.get_unchecked(scan_cursor as usize) } as i32;
                                unsafe { *context.bucket.get_unchecked_mut(scan_cursor as usize) = cell_split_low as u32 };
                                hash_accumulator = ((hash_accumulator ^ 0o65435_u32 as i64) + (scan_cursor + cell_split_low) as i64)
                                    & 0o77777_u32 as i64;
                                if cell_split_high - cell_split_low > largest_bucket_size {
                                    largest_bucket_size = cell_split_high - cell_split_low;
                                    largest_bucket_start = cell_split_low;
                                }
                                if cell_split_low != scan_cell_start {
                                    a_set!(cell_split_low);
                                    if cell_split_high - cell_split_low == 1 {
                                        hint = cell_split_low;
                                    }
                                    *number_of_cells += 1;
                                }
                                if cell_split_high <= scan_cell_end {
                                    // SAFETY: cell_split_high-1 < vertex_count.
                                    unsafe { *partition.get_unchecked_mut((cell_split_high - 1) as usize) = level as u32 };
                                }
                                cell_split_low = cell_split_high;
                            }
                            scan_cursor += 1;
                        }
                        let reorder_start = scan_cell_start;
                        let reorder_end = scan_cell_end;
                        scan_cursor = reorder_start;
                        while scan_cursor <= reorder_end {
                            // SAFETY: scan_cursor < vertex_count, neighbor_counts[scan_cursor] <= vertex_count < bucket.len().
                            let bucket_index = unsafe { *context.bucket.get_unchecked(*neighbor_counts.get_unchecked(scan_cursor as usize) as usize) };
                            unsafe { *context.bucket.get_unchecked_mut(*neighbor_counts.get_unchecked(scan_cursor as usize) as usize) += 1 };
                            // SAFETY: bucket_index < vertex_count (it's a position index).
                            unsafe { *context.work_permutation.get_unchecked_mut(bucket_index as usize) = *labeling.get_unchecked(scan_cursor as usize) };
                            scan_cursor += 1;
                        }
                        scan_cursor = reorder_start;
                        while scan_cursor <= reorder_end {
                            // SAFETY: scan_cursor < vertex_count, work_permutation.len() >= vertex_count.
                            unsafe { *labeling.get_unchecked_mut(scan_cursor as usize) = *context.work_permutation.get_unchecked(scan_cursor as usize) };
                            scan_cursor += 1;
                        }

                        if !a_test!(scan_cell_start) {
                            a_set!(scan_cell_start);
                            a_unset!(largest_bucket_start);
                        }
                    }
                }
                scan_cell_start = scan_cell_end + 1;
            }
        }
    }
    hash_accumulator = ((hash_accumulator ^ 0o65435_u32 as i64) + *number_of_cells as i64) & 0o77777_u32 as i64;
    *hash_code = (hash_accumulator % 0o77777_u32 as i64) as u32;
}

pub fn cheap_automorphism_check(
    partition: &[u32],
    level: u32,
    is_directed_graph: bool,
    number_of_vertices: u32,
) -> bool {
    if is_directed_graph {
        return false;
    }
    let vertex_count = number_of_vertices as usize;
    // cell_count            = count of cell-end markers  (partition[i] <= level)  = total cells.
    // nontrivial_cell_starts = count of nontrivial-cell starts
    //          = count of i where partition[i] > level and
    //            (i == 0 or partition[i-1] <= level).
    // nontrivial_vertex_count = vertex_count - cell_count  (original variable `k`, same semantics).
    //
    // Both are independent reductions over two overlapping stride-1 windows —
    // LLVM auto-vectorizes via two overlapping NEON loads + compare + popcount.
    // SAFETY: vertex_count <= partition.len() (caller guarantee).
    let first = unsafe { *partition.get_unchecked(0) } > level;
    let mut cell_count = (!first) as u32;
    let mut nontrivial_cell_starts = first as u32;
    for cell_index in 1..vertex_count {
        let curr = unsafe { *partition.get_unchecked(cell_index)     } > level;
        let prev = unsafe { *partition.get_unchecked(cell_index - 1) } > level;
        cell_count += (!curr) as u32;
        nontrivial_cell_starts += (curr & !prev) as u32;
    }
    let nontrivial_vertex_count = vertex_count as u32 - cell_count;
    nontrivial_vertex_count <= nontrivial_cell_starts + 1 || nontrivial_vertex_count <= 4
}
fn bestcell(
    graph: &[usize],
    labeling: &[u32],
    partition: &[u32],
    level: i32,
    _target_cell_level: i32,
    words_per_vertex: u32,
    number_of_vertices: u32,
    context: &mut GraphContext,
) -> u32 {
    let mut vertex_index: u32;
    let mut current_graph_index;
    let mut first_cell_index: u32;
    let mut second_cell_index: u32;
    let mut nontrivial_cell_count: u32;
    debug_assert!(context.work_permutation.len() >= number_of_vertices as usize);
    debug_assert!(context.work_bitset.len() >= words_per_vertex as usize);
    debug_assert!(context.bucket.len() >= number_of_vertices as usize + 2);
    // SAFETY: all vertex_index < number_of_vertices; partition invariant ensures inner loop terminates before number_of_vertices;
    // nontrivial_cell_count < number_of_vertices so bucket/work_permutation indices are valid; labeling[vertex_index] < number_of_vertices; first_cell_index,second_cell_index < nontrivial_cell_count.
    nontrivial_cell_count = 0;
    vertex_index = nontrivial_cell_count;
    while vertex_index < number_of_vertices {
        if unsafe { *partition.get_unchecked(vertex_index as usize) } as i32 > level {
            let nontrivial_index = nontrivial_cell_count;
            nontrivial_cell_count += 1;
            unsafe { *context.work_permutation.get_unchecked_mut(nontrivial_index as usize) = vertex_index };
            while unsafe { *partition.get_unchecked(vertex_index as usize) } as i32 > level {
                vertex_index += 1;
            }
        }
        vertex_index += 1;
    }
    if nontrivial_cell_count == 0 {
        return number_of_vertices;
    }
    context.bucket[..(nontrivial_cell_count as usize)].fill(0);
    second_cell_index = 1;
    while second_cell_index < nontrivial_cell_count {
        context.work_bitset[..words_per_vertex as usize].fill(0);
        vertex_index = unsafe { *context.work_permutation.get_unchecked(second_cell_index as usize) } - 1;
        loop {
            vertex_index += 1;
            let member_vertex = unsafe { *labeling.get_unchecked(vertex_index as usize) } as usize;
            unsafe { *context.work_bitset.get_unchecked_mut(member_vertex / usize::BITS as usize) |= 1usize << (usize::BITS - 1 - (member_vertex % usize::BITS as usize) as u32); }
            if unsafe { *partition.get_unchecked(vertex_index as usize) } as i32 <= level {
                break;
            }
        }
        first_cell_index = 0;
        while first_cell_index < second_cell_index {
            let first_cell_vertex = unsafe { *context.work_permutation.get_unchecked(first_cell_index as usize) };
            current_graph_index = words_per_vertex as usize * unsafe { *labeling.get_unchecked(first_cell_vertex as usize) } as usize;
            let (any1, any2) = or_reduce_masked(
                &context.work_bitset,
                unsafe { graph.get_unchecked(current_graph_index..current_graph_index + words_per_vertex as usize) },
                words_per_vertex as usize,
            );
            if any1 && any2 {
                unsafe { *context.bucket.get_unchecked_mut(first_cell_index as usize) += 1 };
                unsafe { *context.bucket.get_unchecked_mut(second_cell_index as usize) += 1 };
            }
            first_cell_index += 1;
        }
        second_cell_index += 1;
    }
    first_cell_index = 0;
    second_cell_index = unsafe { *context.bucket.get_unchecked(0) };
    vertex_index = 1;
    while vertex_index < nontrivial_cell_count {
        let bucket_value = unsafe { *context.bucket.get_unchecked(vertex_index as usize) };
        if bucket_value > second_cell_index {
            first_cell_index = vertex_index;
            second_cell_index = bucket_value;
        }
        vertex_index += 1;
    }
    unsafe { *context.work_permutation.get_unchecked(first_cell_index as usize) }
}

pub fn target_cell(
    graph: &[usize],
    labeling: &[u32],
    partition: &[u32],
    level: i32,
    target_cell_level: i32,
    mut _digraph: bool,
    hint: i32,
    words_per_vertex: u32,
    number_of_vertices: u32,
    context: &mut GraphContext,
) -> i32 {
    // SAFETY: cell_index < number_of_vertices in loop
    let mut cell_index: i32;
    debug_assert!(hint < 0 || (hint as u32) < number_of_vertices, "hint={hint} out of bounds for number_of_vertices={number_of_vertices}");
    if hint >= 0
        && (hint as u32) < number_of_vertices
        // SAFETY: hint >= 0 && hint < number_of_vertices checked above
        && unsafe { *partition.get_unchecked(hint as usize) } as i32 > level
        && (hint == 0 || unsafe { *partition.get_unchecked((hint - 1) as usize) } as i32 <= level)
    {
        hint
    } else if level <= target_cell_level {
        bestcell(
            graph,
            labeling,
            partition,
            level,
            target_cell_level,
            words_per_vertex,
            number_of_vertices,
            context,
        ) as i32
    } else {
        cell_index = 0;
        while cell_index < number_of_vertices as i32 && unsafe { *partition.get_unchecked(cell_index as usize) } as i32 <= level {
            cell_index += 1;
        }
        if cell_index == number_of_vertices as i32 { 0 } else { cell_index }
    }
}

pub fn dense_canonize(
    graph: &mut [usize],
    labeling: &mut [u32],
    partition: &mut [u32],
    orbits: &mut [u32],
    options: &mut CanonautOptions,
    stats: &mut Statistics,
    words_per_vertex: u32,
    number_of_vertices: u32,
    canonical_graph: Option<&mut [usize]>,
    context: &mut GraphContext,
) {
    context.ensure_dnwork_size(words_per_vertex);
    // Split borrow: dense_workspace needs &mut, rest of context needed by canonicalize.
    // Use a temporary to avoid simultaneous mutable borrows.
    let mut dense_workspace = std::mem::take(&mut context.dense_workspace);
    canonicalize(
        graph,
        labeling,
        partition,
        None,
        orbits,
        options,
        stats,
        &mut dense_workspace,
        words_per_vertex,
        number_of_vertices,
        canonical_graph,
        context,
    );
    context.dense_workspace = dense_workspace;
}

pub fn refine_freedyn() {
    // No dynamic buffers to free: GraphContext's Vec fields deallocate themselves on drop.
}
