//! Schreier–Sims stabilizer-chain implementation.
//!
//! [`SchreierSims`] builds and queries a stabilizer chain for the automorphism
//! group discovered during a canonicalization search.  [`SchreierArenaContext`] holds the
//! reusable scratch buffers shared by all operations on a single chain.

use crate::{rng::RngState, structs::BitSetOpsMSB0};
static DEFAULT_SCHREIERFAILS: u32 = 10;

/// Reusable scratch buffers for Schreier–Sims operations.
///
/// Pre-allocated once (sized for `n` vertices) and reused across every filter,
/// expand, and orbit-update call to avoid hot-path allocation.
#[derive(Debug, Default)]
pub struct SchreierArenaContext {
    /// Primary working permutation buffer (`n` elements).
    pub working_permutation: Vec<u32>,
    /// Secondary working permutation buffer (`n` elements).
    pub working_permutation2: Vec<u32>,
    /// Auxiliary buffer `a` for [`SchreierSims::apply_perm`] (`n` elements).
    pub working_permutation_a: Vec<u32>,
    /// Auxiliary buffer `b` for [`SchreierSims::apply_perm`] (`n` elements).
    pub working_permutation_b: Vec<u32>,
    /// Bitset scratch space for orbit marking (`⌈n/64⌉` words).
    pub working_set: Vec<usize>,
    /// Second bitset scratch space (`⌈n/64⌉` words).
    pub working_set2: Vec<usize>,
    /// Number of consecutive random-word failures before `expand_schreier` stops.
    pub maximum_schreier_failures: u32,
    /// Cumulative count of filter calls (diagnostic).
    pub filter_count: u64,
    /// Cumulative count of permutation multiplications (diagnostic).
    pub mult_count: u64,
}

impl SchreierArenaContext {
    /// Creates scratch buffers sized for `vertex_count` vertices.
    pub fn new(vertex_count: usize) -> Self {
        Self {
            working_permutation: vec![0; vertex_count],
            working_permutation2: vec![0; vertex_count],
            working_permutation_a: vec![0; vertex_count],
            working_permutation_b: vec![0; vertex_count],
            // Size hint for bitsets: (vertex_count + 63) / 64
            working_set: vec![0; vertex_count.div_ceil(64)],
            working_set2: vec![0; vertex_count.div_ceil(64)],
            maximum_schreier_failures: 10,
            filter_count: 0,
            mult_count: 0,
        }
    }
    /// Sets the maximum consecutive failure count and returns the previous value.
    ///
    /// Passing `0` resets to the default (`10`).
    pub fn set_schreier_fails(&mut self, number_of_fails: u32) -> u32 {
        let previous_max = self.maximum_schreier_failures;
        if number_of_fails == 0 {
            self.maximum_schreier_failures = DEFAULT_SCHREIERFAILS;
        } else {
            self.maximum_schreier_failures = number_of_fails;
        }
        previous_max
    }
}
/// Typed index into a [`PermArena`].  `u32` saves 4 bytes per slot vs. a pointer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PermId(u32);

impl PermId {
    const NONE: Self = PermId(u32::MAX);
    const IDENTITY: Self = PermId(u32::MAX - 1);
    fn is_valid(self) -> bool {
        self != Self::NONE
    }
}

/// Bump-allocator arena for permutations of fixed length `permutation_length`.
///
/// All permutations are stored in a single contiguous `Vec<u32>` (stride =
/// `permutation_length`). Freed slots are recorded in `free_storage` and
/// reused before appending.
#[derive(Default)]
pub struct PermArena {
    data: Vec<u32>,
    permutation_length: usize,
    free_storage: Vec<u32>,
}

static ARENA_SIZE: usize = 100;
static ARENA_REUSE_SIZE: usize = 64;

impl PermArena {
    /// Creates an arena pre-allocated for `ARENA_SIZE` permutations of length `permutation_length`.
    pub fn new(permutation_length: usize) -> Self {
        Self {
            data: Vec::with_capacity(permutation_length * ARENA_SIZE), // Pre-allocate space
            permutation_length,
            free_storage: Vec::with_capacity(ARENA_REUSE_SIZE),
        }
    }

    /// Pushes a new permutation into the arena and returns its ID.
    /// Reuses memory from deleted permutations if available.
    pub fn alloc(&mut self, permutation: &[u32]) -> PermId {
        debug_assert_eq!(permutation.len(), self.permutation_length, "Permutation size mismatch");

        // STRATEGY 1: Reuse a hole
        if let Some(reused_id) = self.free_storage.pop() {
            let start = (reused_id as usize) * self.permutation_length;
            let end = start + self.permutation_length;

            // Overwrite the old data
            self.data[start..end].copy_from_slice(permutation);

            return PermId(reused_id);
        }

        // STRATEGY 2: Append to the end
        let id = (self.data.len() / self.permutation_length) as u32;
        self.data.extend_from_slice(permutation);
        PermId(id)
    }

    /// Marks a permutation ID as free, allowing its memory to be reused.
    /// Does not physically shrink the vector (capacity remains).
    pub fn free(&mut self, id: PermId) {
        if !id.is_valid() {
            return;
        }

        let slot_index = id.0;
        self.free_storage.push(slot_index);
    }

    #[inline(always)]
    pub fn get(&self, id: PermId) -> &[u32] {
        let start = (id.0 as usize) * self.permutation_length;
        &self.data[start..start + self.permutation_length]
    }

    #[inline(always)]
    pub fn get_mut(&mut self, id: PermId) -> &mut [u32] {
        let start = (id.0 as usize) * self.permutation_length;
        &mut self.data[start..start + self.permutation_length]
    }
}

/// A strongly typed index into the Permutation Arena.
/// Using u32 saves 50% memory on 64-bit systems compared to raw pointers/usize.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FixedPoint(u32);

/// One level of the Schreier–Sims stabilizer chain.
///
/// Fixes `fixed_point` and tracks the coset representatives (transversal) and
/// orbit structure for that stabilizer.
pub struct Level {
    /// The base point fixed at this level.
    pub fixed_point: u32,
    /// Orbit partition: `orbits[v]` is the canonical representative of `v`'s orbit.
    pub orbits: Vec<u32>,
    /// Transversal: a power of `transversal[v]` maps `v` toward the base point,
    /// or `PermId::NONE` if `v` is not yet in the orbit.
    pub transversal: Vec<PermId>,
    /// Exponent mapping `v` to an earlier vertex in the rooted Schreier tree.
    pub powers: Vec<u32>,
}
impl Level {
    fn new(vertex_count: usize, fixed_point: u32) -> Self {
        Self {
            fixed_point,
            orbits: (0..vertex_count as u32).collect(),
            transversal: vec![PermId::NONE; vertex_count],
            powers: vec![0; vertex_count], // <-- Initialize with 0
        }
    }
}

/// Schreier–Sims stabilizer chain for an automorphism group on `n` points.
///
/// Stores a sequence of [`Level`]s, each fixing one more base point.  New
/// generators are added via [`add_generator`](SchreierSims::add_generator) and
/// the chain is expanded probabilistically via
/// [`expand_schreier`](SchreierSims::expand_schreier).
#[derive(Default)]
pub struct SchreierSims {
    number_of_vertices: usize,
    /// Arena owning all stored permutations.
    pub arena: PermArena,
    /// Stabilizer chain (one level per base point).
    pub chain: Vec<Level>,
    /// Shared scratch buffers for all chain operations.
    pub context: SchreierArenaContext,
    /// IDs of the generating permutations in `arena`.
    pub generators: Vec<PermId>,
    essential: Vec<bool>,
    retained: Vec<bool>,
    /// The identity permutation `[0, 1, …, n-1]` for fast identity tests.
    pub identity_slice: Vec<u32>,
    /// RNG driving the random word walk in [`expand_schreier`](Self::expand_schreier)
    /// (same KISS generator as nauty's `ran_nextran`).
    random_generator: RngState,
}

impl SchreierSims {
    /// Creates an empty chain for a group acting on `vertex_count` points.
    pub fn new(vertex_count: usize) -> Self {
        SchreierSims {
            number_of_vertices: vertex_count,
            arena: PermArena::new(vertex_count),
            chain: vec![Level::new(vertex_count, u32::MAX)],
            context: SchreierArenaContext::new(vertex_count),
            generators: Vec::new(),
            essential: Vec::new(),
            retained: Vec::new(),
            identity_slice: (0..vertex_count as u32).collect::<Vec<u32>>(),
            random_generator: RngState::new(),
        }
    }

    /// Applies permutation `permutation` exactly `power` times to each element of `working_permutation` in place.
    ///
    /// Uses direct composition for small `power` (≤ 5), a precomputed cube stride for
    /// `power ≤ 19`, and cycle decomposition for large `power` to stay O(n) in all cases.
    // Nested macro calls like `perm_at!(perm_at!(x))` expand to nested `unsafe` blocks.
    #[allow(unused_unsafe)]
    pub fn apply_perm(
        working_permutation: &mut [u32],
        permutation: &[u32],
        mut power: u32,
        workperm_a: &mut [u32],
        workperm_b: &mut [u32],
        workset: &mut [usize], // Used as a bitset for cycle detection
    ) {
        let permutation_length = permutation.len();
        // SAFETY macros: all *x values are vertex indices < permutation_length; all arrays have len = permutation_length.
        macro_rules! perm_at  { ($v:expr) => { unsafe { *permutation.get_unchecked($v as usize) } } }
        macro_rules! perm_a_at { ($v:expr) => { unsafe { *workperm_a.get_unchecked($v as usize) } } }
        macro_rules! perm_b_at { ($v:expr) => { unsafe { *workperm_b.get_unchecked($v as usize) } } }

        if power <= 5 {
            match power {
                0 => (),
                1 => { for entry in working_permutation.iter_mut() { *entry = perm_at!(*entry); } }
                2 => { for entry in working_permutation.iter_mut() { *entry = perm_at!(perm_at!(*entry)); } }
                3 => { for entry in working_permutation.iter_mut() { *entry = perm_at!(perm_at!(perm_at!(*entry))); } }
                4 => { for entry in working_permutation.iter_mut() { *entry = perm_at!(perm_at!(perm_at!(perm_at!(*entry)))); } }
                5 => { for entry in working_permutation.iter_mut() { *entry = perm_at!(perm_at!(perm_at!(perm_at!(perm_at!(*entry))))); } }
                _ => unreachable!(),
            }
        } else if power <= 19 {
            // Precompute permutation^3 into workperm_a
            for vertex_index in 0..permutation_length {
                let cubed_value = perm_at!(perm_at!(perm_at!(vertex_index as u32)));
                unsafe { *workperm_a.get_unchecked_mut(vertex_index) = cubed_value; }
            }
            // Apply in chunks of 6
            while power >= 6 {
                for entry in working_permutation.iter_mut() { *entry = perm_a_at!(perm_a_at!(*entry)); }
                power -= 6;
            }
            match power {
                0 => {}
                1 => { for entry in working_permutation.iter_mut() { *entry = perm_at!(*entry); } }
                2 => { for entry in working_permutation.iter_mut() { *entry = perm_at!(perm_at!(*entry)); } }
                3 => { for entry in working_permutation.iter_mut() { *entry = perm_a_at!(*entry); } }
                4 => { for entry in working_permutation.iter_mut() { *entry = perm_at!(perm_a_at!(*entry)); } }
                5 => { for entry in working_permutation.iter_mut() { *entry = perm_at!(perm_at!(perm_a_at!(*entry))); } }
                _ => unreachable!(),
            }
        } else {
            // power >= 20: Cycle Decomposition
            workset[..permutation_length.div_ceil(64)].fill(0);
            for vertex_index in 0..permutation_length {
                let word_index = vertex_index / 64;
                let bit_index = vertex_index % 64;
                if (workset[word_index] & (1 << bit_index)) != 0 { continue; }
                if perm_at!(vertex_index as u32) == vertex_index as u32 {
                    unsafe { *workperm_b.get_unchecked_mut(vertex_index) = vertex_index as u32; }
                } else {
                    let mut cycle_length = 1;
                    unsafe { *workperm_a.get_unchecked_mut(0) = vertex_index as u32; }
                    let mut cycle_cursor = perm_at!(vertex_index as u32) as usize;
                    while cycle_cursor != vertex_index {
                        unsafe { *workperm_a.get_unchecked_mut(cycle_length) = cycle_cursor as u32; }
                        cycle_length += 1;
                        workset[cycle_cursor / 64] |= 1 << (cycle_cursor % 64);
                        cycle_cursor = perm_at!(cycle_cursor as u32) as usize;
                    }
                    let mut rotated_offset = (power as usize) % cycle_length;
                    for cycle_position in 0..cycle_length {
                        let original_vertex = unsafe { *workperm_a.get_unchecked(cycle_position) } as usize;
                        let shifted_value = unsafe { *workperm_a.get_unchecked(rotated_offset) };
                        unsafe { *workperm_b.get_unchecked_mut(original_vertex) = shifted_value; }
                        rotated_offset += 1;
                        if rotated_offset == cycle_length { rotated_offset = 0; }
                    }
                }
            }
            for entry in working_permutation.iter_mut() { *entry = perm_b_at!(*entry); }
        }
    }

    /// Start a new graph, retaining permutation and scratch capacity.
    pub fn reset(&mut self, n: usize) {
        self.number_of_vertices = n;
        self.arena.data.clear();
        self.arena.free_storage.clear();
        self.arena.permutation_length = n;
        self.generators.clear();
        self.essential.clear();
        self.retained.clear();
        self.identity_slice.clear();
        self.identity_slice.extend(0..n as u32);
        self.context.working_permutation.resize(n, 0);
        self.context.working_permutation2.resize(n, 0);
        self.context.working_permutation_a.resize(n, 0);
        self.context.working_permutation_b.resize(n, 0);
        self.context.working_set.resize(n.div_ceil(64), 0);
        self.context.working_set2.resize(n.div_ceil(64), 0);
        self.context.maximum_schreier_failures = DEFAULT_SCHREIERFAILS;
        self.context.filter_count = 0;
        self.context.mult_count = 0;
        self.random_generator = RngState::new();
        for level in &mut self.chain {
            level.fixed_point = u32::MAX;
            level.orbits.clear();
            level.orbits.extend(0..n as u32);
            level.transversal.resize(n, PermId::NONE);
            level.transversal.fill(PermId::NONE);
            level.powers.resize(n, 0);
        }
        if self.chain.is_empty() {
            self.chain.push(Level::new(n, u32::MAX));
        }
    }

    /// Adds `partition` as a generator, filtering it through the full stabilizer chain.
    pub fn add_generator(&mut self, partition: &[u32]) -> bool {
        self.filter(partition, false, -1)
    }

    pub fn is_identity(&self, permutation: &[u32]) -> bool {
        self.identity_slice == permutation
    }

    pub fn filter(&mut self, permutation: &[u32], mut ingroup: bool, maximum_level: i32) -> bool {
        // 1. Copy permutation to context.working_permutation
        self.context
            .working_permutation
            .copy_from_slice(permutation);
        self.context.filter_count += 1;

        let mut changed = false;
        let vertex_count = self.number_of_vertices;

        let limit: usize = if maximum_level < 0 {
            self.chain.len()
        } else {
            (maximum_level as usize + 1).min(self.chain.len())
        };

        // Iterate through stabilizer chain levels
        for level_index in 0..limit {
            // Check if working_permutation is identity (if so, we are done)
            if self.is_identity(&self.context.working_permutation) {
                return changed;
            }

            let temp_working_permutation = std::mem::take(&mut self.context.working_permutation);
            // Porting the Orbit Union-Find logic from C:
            // In C, this is the block: "j1 = orbits[i] ... if j1 != j2 { orbits[j2] = j1 }"
            if self.update_orbits(level_index, &temp_working_permutation) {
                changed = true;
            }
            self.context.working_permutation = temp_working_permutation;

            if self.chain[level_index].fixed_point == u32::MAX {
                break;
            }
            let mut current_perm_id = None;

            // 2. Transversal "Hole Filling" Logic
            // SAFETY: vertex_index < vertex_count, working_permutation[vertex_index] < vertex_count, transversal.len() = vertex_count
            for vertex_index in 0..vertex_count {
                let mapped_vertex = unsafe { *self.context.working_permutation.get_unchecked(vertex_index) } as usize;
                if unsafe { *self.chain[level_index].transversal.get_unchecked(vertex_index) } != PermId::NONE
                    && unsafe { *self.chain[level_index].transversal.get_unchecked(mapped_vertex) } == PermId::NONE
                {
                    changed = true;

                    // If we haven't officially added this working_permutation to our storage, do it now
                    if current_perm_id.is_none() {
                        let new_perm_id = self.arena.alloc(&self.context.working_permutation);
                        self.essential.resize(self.arena.data.len() / vertex_count, false);
                        self.essential[new_perm_id.0 as usize] = !ingroup;
                        self.generators.push(new_perm_id);
                        current_perm_id = Some(new_perm_id);
                        ingroup = true; // Once added, it's effectively in the group
                    }

                    // 1. Count how many steps to reach an already-mapped orbit element
                    // SAFETY: orbit_cursor < vertex_count throughout (permutation image of valid vertex)
                    let mut power_needed = 0;
                    let mut orbit_cursor = mapped_vertex;
                    while unsafe { *self.chain[level_index].transversal.get_unchecked(orbit_cursor) } == PermId::NONE {
                        power_needed += 1;
                        orbit_cursor = unsafe { *self.context.working_permutation.get_unchecked(orbit_cursor) } as usize;
                    }

                    // 2. Fill the path backward with the calculated powers
                    orbit_cursor = mapped_vertex;
                    let resolved_perm_id = current_perm_id.unwrap();
                    while unsafe { *self.chain[level_index].transversal.get_unchecked(orbit_cursor) } == PermId::NONE {
                        unsafe { *self.chain[level_index].transversal.get_unchecked_mut(orbit_cursor) = resolved_perm_id; }
                        unsafe { *self.chain[level_index].powers.get_unchecked_mut(orbit_cursor) = power_needed; }
                        power_needed -= 1;
                        orbit_cursor = unsafe { *self.context.working_permutation.get_unchecked(orbit_cursor) } as usize;
                    }
                }
            }
            let fixed_point_index = self.chain[level_index].fixed_point as usize;
            let mut current_image = self.context.working_permutation[fixed_point_index] as usize;

            while current_image != fixed_point_index {
                // SAFETY: current_image < vertex_count, transversal/powers have len = vertex_count
                let transversal_id = unsafe { *self.chain[level_index].transversal.get_unchecked(current_image) };
                debug_assert!(transversal_id != PermId::NONE && transversal_id != PermId::IDENTITY);
                let transversal_power = unsafe { *self.chain[level_index].powers.get_unchecked(current_image) };

                // Read out of arena before splitting context's fields below, or the
                // borrow checker rejects the simultaneous borrows.
                let transversal_permutation = self.arena.get(transversal_id);

                // Borrow split the context arrays safely
                let (working_permutation_slice, working_permutation_a_slice, working_permutation_b_slice, working_set2_slice) = (
                    &mut self.context.working_permutation,
                    &mut self.context.working_permutation_a,
                    &mut self.context.working_permutation_b,
                    &mut self.context.working_set2,
                );

                // O(1) mathematical jump back to the fixed point!
                Self::apply_perm(working_permutation_slice, transversal_permutation, transversal_power, working_permutation_a_slice, working_permutation_b_slice, working_set2_slice);
                self.context.mult_count += 1;
                // SAFETY: fixed_point_index < vertex_count, working_permutation[fixed_point_index] < vertex_count
                current_image = unsafe { *self.context.working_permutation.get_unchecked(fixed_point_index) } as usize;
            }
            // 3. Stripping happens above by repeatedly composing with the transversal element.
        }

        if !self.is_identity(&self.context.working_permutation) && !ingroup {
            // Store the original generator, as in nauty: the stripped residue
            // need not generate the same group without its transversal factors.
            let id = self.arena.alloc(permutation);
            self.essential.resize(self.arena.data.len() / vertex_count, false);
            self.essential[id.0 as usize] = true;
            self.generators.push(id);
            changed = true;
        }

        changed
    }

    /// Generates random words from the current generators and filters them.
    /// Stops when it hits `schreier_fails` consecutive failures.
    /// Returns `true` if the stabilizer chain was expanded.
    pub fn expand_schreier(&mut self) -> bool {
        // If we have no generators, there's nothing to expand.
        if self.generators.is_empty() {
            return false;
        }

        let vertex_count = self.number_of_vertices;
        let mut changed = false;
        let mut failure_count = 0;

        // 1. Detach working_permutation2 from self to satisfy the borrow checker.
        let mut working_permutation2 = std::mem::take(&mut self.context.working_permutation2);

        // 2. Simulate Nauty's random walk on the linked list.
        let mut generator_index = (self.random_generator.next_random() % 17) % self.generators.len();

        let first_generator_id = self.generators[generator_index];
        let first_generator = self.arena.get(first_generator_id);
        working_permutation2.copy_from_slice(first_generator);

        let failure_limit = self.context.maximum_schreier_failures;

        // 3. Keep generating words until we fail too many times
        while failure_count < failure_limit {
            // Pick a word length of 1, 2, or 3
            let word_length = 1 + self.random_generator.next_random() % 3;

            for _ in 0..word_length {
                generator_index = (generator_index + self.random_generator.next_random() % 17) % self.generators.len();
                let generator_id = self.generators[generator_index];

                // Get the generator slice from the arena
                let generator = self.arena.get(generator_id);

                // SAFETY: working_permutation2[vertex_index] < vertex_count = generator.len()
                for vertex_index in 0..vertex_count {
                    let current_image = unsafe { *working_permutation2.get_unchecked(vertex_index) } as usize;
                    unsafe { *working_permutation2.get_unchecked_mut(vertex_index) = *generator.get_unchecked(current_image); }
                }
            }

            // 4. Filter the generated random word.
            // ingroup = true because we constructed it entirely from known generators.
            // max_level = -1 means filter all the way to the bottom.
            let added = self.filter(&working_permutation2, true, -1);

            if added {
                changed = true;
                failure_count = 0; // Reset failures on success!
            } else {
                failure_count += 1;
            }
        }

        // 5. Put working_permutation2 back into the context
        self.context.working_permutation2 = working_permutation2;

        changed
    }
    /// Merges orbits at `level_index` based on the permutation `permutation`.
    /// Returns true if any orbits were merged (changes made).
    /// Port of the Union-Find logic inside `filterschreier`.
    fn update_orbits(&mut self, level_index: usize, permutation: &[u32]) -> bool {
        let vertex_count = self.number_of_vertices;
        let orbits = &mut self.chain[level_index].orbits;
        let mut changed = false;

        // SAFETY: all orbit values are valid vertex indices < vertex_count; permutation is a permutation on 0..vertex_count
        for vertex_index in 0..vertex_count {
            let mut root_of_vertex = unsafe { *orbits.get_unchecked(vertex_index) } as usize;
            while unsafe { *orbits.get_unchecked(root_of_vertex) } as usize != root_of_vertex {
                root_of_vertex = unsafe { *orbits.get_unchecked(root_of_vertex) } as usize;
            }
            let mapped_vertex = unsafe { *permutation.get_unchecked(vertex_index) } as usize;
            let mut root_of_mapped_vertex = unsafe { *orbits.get_unchecked(mapped_vertex) } as usize;
            while unsafe { *orbits.get_unchecked(root_of_mapped_vertex) } as usize != root_of_mapped_vertex {
                root_of_mapped_vertex = unsafe { *orbits.get_unchecked(root_of_mapped_vertex) } as usize;
            }
            if root_of_vertex != root_of_mapped_vertex {
                changed = true;
                if root_of_vertex < root_of_mapped_vertex {
                    unsafe { *orbits.get_unchecked_mut(root_of_mapped_vertex) = root_of_vertex as u32; }
                } else {
                    unsafe { *orbits.get_unchecked_mut(root_of_vertex) = root_of_mapped_vertex as u32; }
                }
            }
        }

        if changed {
            for vertex_index in 0..vertex_count {
                let root_value = unsafe { *orbits.get_unchecked(vertex_index) } as usize;
                unsafe { *orbits.get_unchecked_mut(vertex_index) = *orbits.get_unchecked(root_value); }
            }
        }

        changed
    }

    /// Checks if the permutation `permutation` exists in the current set of generators.
    /// This replaces the C function `findpermutation`.
    ///
    /// Returns `true` if found, `false` otherwise.
    pub fn find_permutation(&self, permutation: &[u32]) -> Option<PermId> {
        // In Rust, slice comparison (a == b) checks contents, which is equivalent
        // to the while loop in C. It uses optimized memcmp under the hood.
        self.generators
            .iter()
            .find(|&&generator_id| self.arena.get(generator_id) == permutation)
            .copied()
    }

    /// Prune using only automorphisms fixing every point in `fixed_points`.
    /// Random expansion can miss stabilizer elements, which only reduces pruning.
    pub fn prune_set(&mut self, fixed_points: &[usize], cell: &mut [usize]) {
        let n = self.number_of_vertices;
        self.context.working_set.copy_from_slice(&fixed_points[..n.div_ceil(64)]);
        let mut depth = 0;
        while self.chain[depth].fixed_point != u32::MAX {
            let point = self.chain[depth].fixed_point as usize;
            if !self.context.working_set.is_bit_set(point) { break; }
            self.context.working_set.unset_bit(point);
            depth += 1;
        }
        if let Some(mut point) = self.context.working_set.next_set_bit(0) {
            // Preserve the compatible prefix and this level's orbits. Its
            // transversal and every deeper stabilizer depend on the new base.
            for level in &mut self.chain[depth..] {
                level.transversal.fill(PermId::NONE);
            }
            loop {
                self.chain[depth].fixed_point = point as u32;
                self.chain[depth].transversal[point] = PermId::IDENTITY;
                depth += 1;
                if depth == self.chain.len() {
                    self.chain.push(Level::new(n, u32::MAX));
                }
                let level = &mut self.chain[depth];
                level.fixed_point = u32::MAX;
                for (v, orbit) in level.orbits.iter_mut().enumerate() { *orbit = v as u32; }
                match self.context.working_set.next_set_bit(point + 1) {
                    Some(next) => point = next,
                    None => break,
                }
            }
            // Release derived permutations no longer used by any transversal.
            // Essential generators remain, so changing the base preserves the group.
            self.retained.clone_from(&self.essential);
            for level in &self.chain {
                for &id in &level.transversal {
                    if id != PermId::NONE && id != PermId::IDENTITY {
                        self.retained[id.0 as usize] = true;
                    }
                }
            }
            self.generators.retain(|&id| {
                if self.retained[id.0 as usize] { true }
                else { self.arena.free(id); false }
            });
            self.expand_schreier();
        }
        let orbits = &self.chain[depth].orbits;
        let mut next = cell.next_set_bit(0);
        while let Some(v) = next {
            if orbits[v] != v as u32 { cell.unset_bit(v); }
            next = cell.next_set_bit(v + 1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn closure(generators: &[Vec<u32>], n: usize) -> Vec<Vec<u32>> {
        let mut group = vec![(0..n as u32).collect::<Vec<_>>()];
        let mut seen = HashSet::from([group[0].clone()]);
        let mut i = 0;
        while i < group.len() {
            for g in generators {
                let p: Vec<_> = group[i].iter().map(|&v| g[v as usize]).collect();
                if seen.insert(p.clone()) { group.push(p); }
            }
            i += 1;
        }
        group
    }

    #[test]
    fn pruning_respects_every_pointwise_stabilizer() {
        let mut rng = RngState::new();
        let mut chain = SchreierSims::default();
        for n in [5, 3, 6, 1, 4, 0, 5] {
            for _ in 0..20 {
                chain.reset(n);
                let mut generators = Vec::new();
                for _ in 0..3 {
                    let mut g: Vec<_> = (0..n as u32).collect();
                    for i in 1..n { g.swap(i, rng.next_random() % (i + 1)); }
                    chain.add_generator(&g);
                    generators.push(g);
                    let group = closure(&generators, n);
                    // Deliberately switch incompatible bases and shorten them.
                    for mask in (0..1usize << n).rev() {
                        let mut fixed = vec![0usize; n.div_ceil(64)];
                        let mut cell = fixed.clone();
                        for v in 0..n {
                            cell.set_bit(v);
                            if mask & (1 << v) != 0 { fixed.set_bit(v); }
                        }
                        chain.prune_set(&fixed, &mut cell);
                        for v in 0..n {
                            if !cell.is_bit_set(v) {
                                assert!(group.iter().any(|p| {
                                    p[v] < v as u32 && (0..n).all(|w| mask & (1 << w) == 0 || p[w] == w as u32)
                                }), "unsound pruning: n={n}, mask={mask}, v={v}");
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn symmetric_group_actually_prunes_below_fixed_points() {
        let mut chain = SchreierSims::new(6);
        chain.add_generator(&[1, 2, 3, 4, 5, 0]);
        chain.add_generator(&[1, 0, 2, 3, 4, 5]);
        let mut fixed = vec![0usize];
        fixed.set_bit(0);
        fixed.set_bit(2);
        let mut cell = vec![0usize];
        for v in [1, 3, 4, 5] { cell.set_bit(v); }
        chain.prune_set(&fixed, &mut cell);
        assert_eq!(cell, vec![1usize << 62]);
        assert!(chain.context.mult_count > 0);
    }

    #[test]
    fn permutation_powers_match_repeated_composition() {
        let n = 71;
        let permutation: Vec<_> = (0..n).map(|i| ((i + 1) % n) as u32).collect();
        let mut context = SchreierArenaContext::new(n);
        for power in [0, 1, 2, 3, 4, 5, 6, 7, 11, 19, 20, 70, 71, 137] {
            let mut actual: Vec<_> = (0..n as u32).rev().collect();
            let mut expected = actual.clone();
            for _ in 0..power {
                for v in &mut expected { *v = permutation[*v as usize]; }
            }
            SchreierSims::apply_perm(&mut actual, &permutation, power,
                &mut context.working_permutation_a, &mut context.working_permutation_b,
                &mut context.working_set2);
            assert_eq!(actual, expected);
        }
    }
}
