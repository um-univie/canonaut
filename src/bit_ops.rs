pub const WORDSIZE: u32 = usize::BITS;
/// Bits per word.
pub const W: usize = usize::BITS as usize;
/// `x >> WSHIFT` is x's word; `x & WMASK` its bit in that word.
pub const WSHIFT: u32 = usize::BITS.trailing_zeros();
pub const WMASK: usize = W - 1;

/// Precomputed MSB-0 single-bit masks: `BIT_USIZE[i] == 1 << (63 - i)`.
/// Avoids repeated shifts in support.rs's hot loops.
pub static BIT_USIZE: [usize; WORDSIZE as usize] = {
    let mut bit_masks = [0; WORDSIZE as usize];
    let mut bit_index = 0;
    while bit_index < WORDSIZE as usize {
        bit_masks[bit_index] = 1usize << (WORDSIZE as usize - 1 - bit_index);
        bit_index += 1;
    }
    bit_masks
};
