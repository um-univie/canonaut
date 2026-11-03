//! Low-level bit-manipulation utilities and hash helpers.

const WORD_SIZE: usize = usize::BITS as usize;
const WORD_SHIFT: u32 = WORD_SIZE.trailing_zeros();

/// Union-find orbit join: merges orbits identified by `permutation` and returns the orbit count.
///
/// Implements `orbjoin` from the original nauty C source.  After the call,
/// `orbits[v]` contains the canonical representative of `v`'s orbit (path
/// compressed).  Returns the number of distinct orbits.
///
/// # Arguments
/// * `orbits` — orbit partition array; `orbits[v]` is the parent of `v`.
/// * `permutation` — a permutation; elements `i` and `permutation[i]` are merged into one orbit.
/// * `vertex_count` — number of elements; both slices must be at least `vertex_count` long.
#[inline]
pub(crate) fn safe_orbjoin(orbits: &mut [u32], permutation: &[u32], vertex_count: u32) -> u32 {
    if vertex_count == 0 || orbits.len() < vertex_count as usize || permutation.len() < vertex_count as usize {
        return 0;
    }

    // SAFETY: all orbit values and permutation values are valid vertex indices in 0..vertex_count;
    // orbit roots converge because orbits[x] <= x for all non-root x (union by rank).
    for vertex_index in 0..vertex_count as usize {
        let mapped_vertex = unsafe { *permutation.get_unchecked(vertex_index) };
        if mapped_vertex != vertex_index as u32 {
            let mut first_root = unsafe { *orbits.get_unchecked(vertex_index) };
            while unsafe { *orbits.get_unchecked(first_root as usize) } != first_root {
                first_root = unsafe { *orbits.get_unchecked(first_root as usize) };
            }
            let mut second_root = unsafe { *orbits.get_unchecked(mapped_vertex as usize) };
            while unsafe { *orbits.get_unchecked(second_root as usize) } != second_root {
                second_root = unsafe { *orbits.get_unchecked(second_root as usize) };
            }
            if first_root != second_root {
                if first_root < second_root {
                    unsafe { *orbits.get_unchecked_mut(second_root as usize) = first_root };
                } else {
                    unsafe { *orbits.get_unchecked_mut(first_root as usize) = second_root };
                }
            }
        }
    }

    let mut orbit_count = 0u32;
    for vertex_index in 0..vertex_count as usize {
        let mut root = unsafe { *orbits.get_unchecked(vertex_index) };
        while unsafe { *orbits.get_unchecked(root as usize) } != root {
            root = unsafe { *orbits.get_unchecked(root as usize) };
        }
        unsafe { *orbits.get_unchecked_mut(vertex_index) = root };
        if root == vertex_index as u32 {
            orbit_count += 1;
        }
    }

    orbit_count
}

/// Finds the next set bit in `bitset` strictly after position `position` (MSB-0 ordering).
///
/// Pass `position = -1` to start from bit 0.  Returns `None` when no further set bit exists.
/// This is the safe Rust equivalent of nauty's `nextelement` macro.
///
/// # Arguments
/// * `bitset` — slice of `usize` words, MSB-0 encoded.
/// * `words_per_vertex` — number of words to scan; must be `≤ bitset.len()`.
/// * `position` — last returned position, or `-1` for the first call.
pub fn next_element(bitset: &[usize], words_per_vertex: u32, position: i32) -> Option<i32> {
    if words_per_vertex == 0 || bitset.len() < words_per_vertex as usize {
        return None;
    }

    let mut remaining_word: usize;
    let mut word_index: i32;

    // Handle single word case
    if words_per_vertex == 1 {
        if position < 0 {
            remaining_word = bitset[0];
        } else {
            if position as usize >= WORD_SIZE - 1 {
                // Changed from 64 to 63 to prevent overflow
                return None;
            }
            remaining_word = bitset[0] & ((usize::MAX >> 1) >> position);
        }

        if remaining_word == 0 {
            return None;
        } else {
            return Some(remaining_word.leading_zeros() as i32);
        }
    }

    // Handle multi-word case
    if position < 0 {
        word_index = 0;
        remaining_word = bitset[0];
    } else {
        word_index = position >> WORD_SHIFT; // position / word size
        if word_index >= words_per_vertex as i32 {
            return None;
        }
        remaining_word = bitset[word_index as usize] & ((usize::MAX >> 1) >> (position as usize & (WORD_SIZE - 1)));
    }

    loop {
        if remaining_word != 0 {
            let bit_position = remaining_word.leading_zeros() as i32;
            let absolute_position = (word_index << WORD_SHIFT) + bit_position;
            // Ensure we don't return a position beyond the valid range
            if absolute_position >= (words_per_vertex << WORD_SHIFT) as i32 {
                return None;
            }
            return Some(absolute_position);
        }

        word_index += 1;
        if word_index >= words_per_vertex as i32 {
            return None;
        }
        remaining_word = bitset[word_index as usize];
    }
}
/// Returns the number of `usize` words needed to store `number_of_bits` bits.
///
/// Equivalent to nauty's `SETWORDSNEEDED` macro: `⌈number_of_bits / 64⌉` on
/// 64-bit targets.  Returns `0` for zero bits.
///
/// # Example
/// ```rust
/// use canonaut::utilities::words_needed;
/// assert_eq!(words_needed(0),   0);
/// assert_eq!(words_needed(1),   1);
/// assert_eq!(words_needed(64),  1);
/// assert_eq!(words_needed(65),  2);
/// assert_eq!(words_needed(128), 2);
/// ```
#[inline(always)]
pub fn words_needed(number_of_bits: usize) -> usize {
    if number_of_bits == 0 {
        0
    } else {
        ((number_of_bits - 1) >> WORD_SHIFT) + 1
    }
}

const MASH_XOR_CONST: i64 = 0o65435;
const MASH_MASK: i64 = 0o77777;

/// Accumulates a 15-bit running hash: `((accumulator ^ 0o65435) + increment) & 0o77777`.
///
/// Direct port of nauty's `MASH` macro, used to produce invariant codes during
/// partition refinement.
#[inline(always)]
pub fn mash(accumulator: i64, increment: i64) -> i64 {
    ((accumulator ^ MASH_XOR_CONST) + increment) & MASH_MASK
}

/// Reduces a hash accumulator to the 15-bit range: `accumulator % 0o77777`.
///
/// Direct port of nauty's `CLEANUP` macro.
#[inline(always)]
pub fn cleanup(accumulator: i64) -> i64 {
    accumulator % MASH_MASK
}
