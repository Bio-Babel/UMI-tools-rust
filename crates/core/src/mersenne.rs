//! `CPython`'s Mersenne Twister, reproducing `random.random()` exactly
//! (Python dependency: the standard library `random` module).
//!
//! `Utilities.Start` calls `random.seed(options.random_seed)` when
//! `--random-seed` is given, and `sam_methods.get_bundles.update_dicts` then
//! calls `random.random()` to break ties between equally good reads. So the
//! choice of representative read is fully determined by the seed — and the
//! package's own test suite depends on it, comparing deduplicated output
//! byte for byte against checked-in references produced with
//! `--random-seed=123456789`.
//!
//! Matching those references therefore means matching `CPython`'s generator
//! exactly, not merely drawing from the same distribution. That is what this
//! module does: MT19937 with `CPython`'s `init_by_array` seeding and its
//! `genrand_res53` double conversion.
//!
//! # Why not a general-purpose RNG crate
//!
//! `rand`'s `StdRng` is a different algorithm, and even `rand`'s MT
//! implementations do not reproduce `CPython`'s *seeding*: `random.seed(n)` for
//! an integer splits `abs(n)` into 32-bit little-endian words and runs
//! `init_by_array` over them, which is a Python-specific convention rather
//! than part of the MT19937 reference. Reproducing it here is about sixty
//! lines and is verified against the live interpreter.

/// Degree of recurrence — the state size in 32-bit words.
const N: usize = 624;
/// Middle word offset.
const M: usize = 397;
/// Twist matrix coefficient.
const MATRIX_A: u32 = 0x9908_b0df;
/// Most significant `w - r` bits.
const UPPER_MASK: u32 = 0x8000_0000;
/// Least significant `r` bits.
const LOWER_MASK: u32 = 0x7fff_ffff;

/// `CPython`'s `random.Random`, restricted to seeding and `random()`.
///
/// # Examples
///
/// ```
/// use umi_tools_core::mersenne::MersenneTwister;
///
/// // The same values CPython produces for `random.seed(1)`.
/// let mut rng = MersenneTwister::seed_from_u64(1);
/// assert!((rng.random() - 0.134_364_244_112_401_22).abs() < 1e-15);
/// assert!((rng.random() - 0.847_433_736_937_232_7).abs() < 1e-15);
/// ```
#[derive(Debug, Clone)]
pub struct MersenneTwister {
    state: [u32; N],
    index: usize,
}

impl MersenneTwister {
    /// Seeds the generator as `CPython`'s `random.seed(n)` does for a
    /// non-negative integer.
    ///
    /// `CPython` converts the integer to its absolute value, splits it into
    /// 32-bit little-endian words, and runs `init_by_array` over them. A seed
    /// of zero yields the single-word key `[0]`, not an empty one.
    #[must_use]
    pub fn seed_from_u64(seed: u64) -> Self {
        let mut key: Vec<u32> = Vec::new();
        let mut remaining = seed;
        loop {
            key.push((remaining & 0xffff_ffff) as u32);
            remaining >>= 32;
            if remaining == 0 {
                break;
            }
        }
        Self::init_by_array(&key)
    }

    /// Seeds the generator as `NumPy`'s legacy `numpy.random.seed(n)` does.
    ///
    /// `NumPy` and `CPython` seed the *same* algorithm differently: `NumPy`'s
    /// `RandomState.seed` calls the MT19937 reference `init_genrand` on the
    /// low 32 bits, while `CPython`'s `random.seed` runs `init_by_array` over
    /// the seed's 32-bit words. `dedup --output-stats` uses both — the stdlib
    /// stream for the tie-break, `NumPy`'s for the null distribution — so both
    /// seedings are needed, and mixing them up silently changes the null.
    ///
    /// # Examples
    ///
    /// ```
    /// use umi_tools_core::mersenne::MersenneTwister;
    ///
    /// // The values NumPy produces for `np.random.seed(0)`.
    /// let mut rng = MersenneTwister::seed_numpy(0);
    /// assert!((rng.random() - 0.548_813_503_927_324_8).abs() < 1e-15);
    /// assert!((rng.random() - 0.715_189_366_372_419_5).abs() < 1e-15);
    /// ```
    #[must_use]
    pub fn seed_numpy(seed: u32) -> Self {
        Self::init_genrand(seed)
    }

    /// The MT19937 reference `init_genrand`.
    fn init_genrand(seed: u32) -> Self {
        let mut state = [0u32; N];
        state[0] = seed;
        for i in 1..N {
            state[i] = 1_812_433_253u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Self { state, index: N }
    }

    /// The MT19937 reference `init_by_array`, which `CPython` uses for integer
    /// seeds.
    #[must_use]
    pub fn init_by_array(key: &[u32]) -> Self {
        let mut rng = Self::init_genrand(19_650_218);
        let mut i = 1usize;
        let mut j = 0usize;
        let mut k = N.max(key.len());

        while k > 0 {
            let prev = rng.state[i - 1];
            rng.state[i] = (rng.state[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                rng.state[0] = rng.state[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
            k -= 1;
        }

        let mut k = N - 1;
        while k > 0 {
            let prev = rng.state[i - 1];
            rng.state[i] = (rng.state[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1_566_083_941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                rng.state[0] = rng.state[N - 1];
                i = 1;
            }
            k -= 1;
        }

        // Guarantees a non-zero initial state.
        rng.state[0] = 0x8000_0000;
        rng.index = N;
        rng
    }

    /// The next 32-bit output — MT19937's `genrand_uint32`.
    #[must_use]
    pub fn next_u32(&mut self) -> u32 {
        if self.index >= N {
            self.twist();
        }
        let mut y = self.state[self.index];
        self.index += 1;
        // Tempering.
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    fn twist(&mut self) {
        for i in 0..N {
            let y = (self.state[i] & UPPER_MASK) | (self.state[(i + 1) % N] & LOWER_MASK);
            let mut next = self.state[(i + M) % N] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= MATRIX_A;
            }
            self.state[i] = next;
        }
        self.index = 0;
    }

    /// The next double in `[0, 1)` — `CPython`'s `random.random()`.
    ///
    /// Uses `genrand_res53`: two 32-bit outputs combined into 53 bits of
    /// mantissa, which is what gives a Python float its full precision.
    #[must_use]
    pub fn random(&mut self) -> f64 {
        let a = f64::from(self.next_u32() >> 5);
        let b = f64::from(self.next_u32() >> 6);
        (a * 67_108_864.0 + b) / 9_007_199_254_740_992.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference values taken from the live `CPython` interpreter:
    /// `random.seed(n); [random.random() for _ in range(5)]`.
    #[test]
    fn matches_cpython_for_known_seeds() {
        let cases: &[(u64, [f64; 5])] = &[
            (
                123_456_789,
                [
                    0.641_400_616_185_872_6,
                    0.542_189_268_096_949_5,
                    0.993_175_066_283_272_1,
                    0.843_252_136_686_916_6,
                    0.811_733_928_337_940_6,
                ],
            ),
            (
                0,
                [
                    0.844_421_851_525_048_1,
                    0.757_954_402_940_302_5,
                    0.420_571_580_830_845,
                    0.258_916_750_292_963_35,
                    0.511_274_721_368_608_5,
                ],
            ),
            (
                1,
                [
                    0.134_364_244_112_401_22,
                    0.847_433_736_937_232_7,
                    0.763_774_618_976_614,
                    0.255_069_025_739_421_7,
                    0.495_435_087_091_940_95,
                ],
            ),
            (
                1_099_511_627_783,
                [
                    0.613_703_777_993_651_1,
                    0.814_916_297_330_948_7,
                    0.945_011_508_759_287_3,
                    0.570_390_699_100_456_1,
                    0.475_924_098_870_489_36,
                ],
            ),
        ];

        for (seed, expected) in cases {
            let mut rng = MersenneTwister::seed_from_u64(*seed);
            for (i, want) in expected.iter().enumerate() {
                let got = rng.random();
                assert!(
                    (got - want).abs() < 1e-15,
                    "seed {seed}, draw {i}: {got} != {want}"
                );
            }
        }
    }

    /// Reference values from `np.random.seed(n); np.random.random_sample(5)`.
    #[test]
    fn numpy_seeding_matches_numpy() {
        let cases: &[(u32, [f64; 5])] = &[
            (
                123_456_789,
                [
                    0.532_833_024_789_759,
                    0.534_136_600_890_416_6,
                    0.509_553_035_526_467,
                    0.713_564_031_832_761_6,
                    0.256_998_954_756_975_43,
                ],
            ),
            (
                0,
                [
                    0.548_813_503_927_324_8,
                    0.715_189_366_372_419_5,
                    0.602_763_376_071_643_9,
                    0.544_883_182_996_896_9,
                    0.423_654_799_338_904_7,
                ],
            ),
        ];
        for (seed, expected) in cases {
            let mut rng = MersenneTwister::seed_numpy(*seed);
            for (i, want) in expected.iter().enumerate() {
                let got = rng.random();
                assert!(
                    (got - want).abs() < 1e-15,
                    "numpy seed {seed}, draw {i}: {got} != {want}"
                );
            }
        }
    }

    /// The two seedings must not coincide, or one of them is wrong.
    #[test]
    fn the_two_seedings_differ() {
        let mut cpython = MersenneTwister::seed_from_u64(123_456_789);
        let mut numpy = MersenneTwister::seed_numpy(123_456_789);
        assert!((cpython.random() - numpy.random()).abs() > 1e-9);
    }

    #[test]
    fn draws_stay_in_the_unit_interval() {
        let mut rng = MersenneTwister::seed_from_u64(42);
        for _ in 0..10_000 {
            let x = rng.random();
            assert!((0.0..1.0).contains(&x), "{x} outside [0, 1)");
        }
    }

    #[test]
    fn the_same_seed_replays_the_same_stream() {
        let a: Vec<f64> = (0..100)
            .scan(MersenneTwister::seed_from_u64(7), |rng, _| {
                Some(rng.random())
            })
            .collect();
        let b: Vec<f64> = (0..100)
            .scan(MersenneTwister::seed_from_u64(7), |rng, _| {
                Some(rng.random())
            })
            .collect();
        assert_eq!(a, b);
        let c: Vec<f64> = (0..100)
            .scan(MersenneTwister::seed_from_u64(8), |rng, _| {
                Some(rng.random())
            })
            .collect();
        assert_ne!(a, c);
    }

    #[test]
    fn the_stream_survives_a_state_refill() {
        // 624 outputs exhaust one block; the 625th forces a twist.
        let mut rng = MersenneTwister::seed_from_u64(1);
        for _ in 0..2000 {
            let _ = rng.next_u32();
        }
        assert!((0.0..1.0).contains(&rng.random()));
    }

    #[test]
    fn the_mean_of_many_draws_is_about_half() {
        let mut rng = MersenneTwister::seed_from_u64(99);
        let n = 100_000;
        let mean: f64 = (0..n).map(|_| rng.random()).sum::<f64>() / f64::from(n);
        assert!((mean - 0.5).abs() < 0.01, "mean = {mean}");
    }
}
