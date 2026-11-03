/// State for nauty's KISS generator (Marsaglia), owned per-instance instead of
/// static mutables so multiple generators can run independently.
#[derive(Debug, Clone)]
pub struct RngState {
    kissx: u64,
    kissc: u64,
    kissy: u64,
    kissz: u64,
    kisst: u64,
}

impl RngState {
    pub fn new() -> Self {
        Self {
            kissx: 1234567890987654321_u64,
            kissc: 123456123456123456_u64,
            kissy: 362436362436362436_u64,
            kissz: 1066149217761810_u64,
            kisst: 0,
        }
    }

    pub fn init(&mut self, seed: usize) {
        self.init_with_two_seeds(seed, 0);
    }

    pub fn init_with_two_seeds(&mut self, seed1: usize, seed2: usize) {
        self.kissx = 1234567890987654321_u64.wrapping_add(seed1 as u64);
        self.kissy = 362436362436362436_u64.wrapping_add((seed2 as u64).wrapping_mul(997));
        self.kissc = 123456123456123456_u64;
        self.kissz = 1066149217761810_u64;

        // Warm up the generator
        for _ in 0..1000 {
            self.next_random();
        }
    }

    /// Seeds from the current time mixed with `extra`; returns the seed used.
    /// Doesn't need to match nauty's C `initialize` bit-for-bit, just comparable entropy.
    pub fn init_with_time(&mut self, extra: usize) -> usize {
        let time_val = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_micros() as f64
            * 1e-6;

        let seed = if time_val > 1660000000.0 {
            (time_val * 2100001.0) as usize
        } else {
            (time_val + 212300021.0) as usize
        };

        self.init_with_two_seeds(seed, extra);
        seed
    }

    pub fn next_random(&mut self) -> usize {
        self.kisst = (self.kissx << 58).wrapping_add(self.kissc);
        self.kissc = self.kissx >> 6;
        self.kissx = self.kissx.wrapping_add(self.kisst);
        self.kissc = self.kissc.wrapping_add((self.kissx < self.kisst) as u64);
        self.kissy ^= self.kissy << 13;
        self.kissy ^= self.kissy >> 17;
        self.kissy ^= self.kissy << 43;
        self.kissz = 6906969069_u64
            .wrapping_mul(self.kissz)
            .wrapping_add(1234567);

        self.kissx.wrapping_add(self.kissy).wrapping_add(self.kissz) as usize
    }
}

impl Default for RngState {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rng_state_initialization() {
        let generator = RngState::new();
        assert_eq!(generator.kissx, 1234567890987654321);
        assert_eq!(generator.kissc, 123456123456123456);
        assert_eq!(generator.kissy, 362436362436362436);
        assert_eq!(generator.kissz, 1066149217761810);
        assert_eq!(generator.kisst, 0);
    }

    #[test]
    fn test_rng_state_seeding() {
        let mut generator = RngState::new();
        generator.init(12345);

        // After initialization with seed, values should be different from defaults
        assert_ne!(generator.kissx, 1234567890987654321);
        assert_ne!(generator.kissy, 362436362436362436);
    }

    #[test]
    fn test_rng_randomness() {
        let mut generator = RngState::new();
        generator.init(42);

        let first = generator.next_random();
        let second = generator.next_random();
        let third = generator.next_random();

        // Should generate different values
        assert_ne!(first, second);
        assert_ne!(second, third);
        assert_ne!(first, third);
    }

    #[test]
    fn test_rng_deterministic() {
        let mut first_generator = RngState::new();
        let mut second_generator = RngState::new();

        first_generator.init(12345);
        second_generator.init(12345);

        // Same seed should produce same sequence
        for _ in 0..10 {
            assert_eq!(first_generator.next_random(), second_generator.next_random());
        }
    }

    #[test]
    fn test_two_seed_initialization() {
        let mut generator = RngState::new();
        generator.init_with_two_seeds(111, 222);

        let first = generator.next_random();
        let second = generator.next_random();

        assert_ne!(first, second);
    }
}
