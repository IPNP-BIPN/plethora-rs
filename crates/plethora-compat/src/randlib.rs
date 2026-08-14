//! RANDLIB, as Perl's `Math::Random` compiles it.
//!
//! `merge_pairs.pl` reseeds the generator from an MD5 of every input line and
//! draws one normal deviate to decide how far to extend an unpaired read:
//!
//! ```text
//! my $eed = md5_hex($line);
//! random_set_seed_from_phrase($eed);
//! my $inner_distance = sprintf("%.0f", random_normal(1, $mean, $sd) / 2);
//! ```
//!
//! That makes the extension deterministic per read, and it makes reproducing
//! plethora's BED output a question of reproducing RANDLIB exactly. The chain
//! is `random_set_seed_from_phrase` -> `salfph` -> `phrtsd` + `setall`, then
//! `random_normal` -> `gennor` -> `sd * snorm() + av`.
//!
//! Transcribed from Math-Random-0.75 (sha256
//! 72f9c94c32fcc6dcda16bd2a5e50f6562386d2afc528e5efb369aa2e9596025b):
//! `randlib.c` for [`Randlib::phrtsd`], [`Randlib::snorm`], [`Randlib::ranf`]
//! and `mltmod`; `com.c` for [`Randlib::setall`], `initgn` and
//! [`Randlib::ignlgi`]; `helper.c` for `salfph`.
//!
//! The generator itself is L'Ecuyer and Cote, "Implementing a Random Number
//! Package with Splitting Facilities", ACM TOMS 17:98-111 (1991). The normal
//! deviate is Ahrens and Dieter, "Extensions of Forsythe's Method for Random
//! Sampling from the Normal Distribution", Math. Comput. 27:927-937 (1973),
//! algorithm FL with M = 5.

/// Number of virtual generators RANDLIB maintains. Plethora only ever uses the
/// first, but `setall` walks all of them, so the walk is reproduced.
const NUMG: usize = 32;

const M1: i64 = 2147483563;
const M2: i64 = 2147483399;
const A1: i64 = 40014;
const A2: i64 = 40692;
/// `A1**(2**(V+W)) mod M1` with V = 20, W = 30, precomputed in `inrgcm`.
const A1VW: i64 = 2082007225;
/// `A2**(2**(V+W)) mod M2`.
const A2VW: i64 = 784306273;

/// `1/M1`, at the precision `ranf` was widened to in randlib.c ("WGR, 2/12/01:
/// increased precision"). Using the older 4.656613057E-10 shifts every draw.
const INV_M1: f64 = 4.656_613_057_391_77e-10;

const TWOP30: i64 = 1073741824;

/// Which `phrtsd` the local `Math::Random` was built with.
///
/// `randlib.c` carries two implementations behind `#ifdef PHRTSD_ORIG`, and
/// `Makefile.PL` only defines it when the build was invoked as
/// `perl Makefile.PL phrtsd_orig`. A plain `cpanm Math::Random` therefore
/// compiles [`Phrtsd::New`], but a distribution packaged by hand may not, and
/// the two disagree on every phrase. The parity test resolves which one is live
/// rather than assuming.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Phrtsd {
    /// The `#else` branch: eight magic multipliers, and a loop bound of
    /// `lphr - 1` that silently ignores the final character of the phrase.
    #[default]
    New,
    /// The `#ifdef` branch: each character is looked up in an 88-glyph table
    /// and folded into the seeds through five shifts.
    Orig,
}

/// The character table the original `phrtsd` indexes into.
///
/// Written in C as four continued string literals; `\\\"` there is a single
/// backslash followed by a single double quote. The trailing space is
/// deliberate (`WGR added space, 5/19/1999`) and is the last glyph, which is
/// what makes the `if (!table[ix]) ix = 0;` check below reachable.
const TABLE: &[u8] =
    br#"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$%^&*()_+[];:'\"<>?,./ "#;

const SHIFT: [i64; 5] = [1, 64, 4096, 262144, 16777216];

/// The `#else` branch's multipliers.
const VALUES_NEW: [i64; 8] = [
    8521739, 5266711, 3254959, 2011673, 1243273, 768389, 474899, 293507,
];

/// RANDLIB's global state, made explicit.
///
/// In C this lives in a common block shared by the whole process, which is why
/// `merge_pairs.pl` can reseed it from anywhere. Here it is an owned value, so
/// parallel callers cannot interfere with each other. Sequential behaviour is
/// unchanged: plethora reseeds before every draw.
#[derive(Debug, Clone)]
pub struct Randlib {
    /// Initial seeds (`Xig1`, `Xig2`).
    ig: [(i64, i64); NUMG],
    /// Block-start seeds (`Xlg1`, `Xlg2`).
    lg: [(i64, i64); NUMG],
    /// Current seeds (`Xcg1`, `Xcg2`).
    cg: [(i64, i64); NUMG],
    /// Current generator, 1-based exactly as `gscgn` reports it.
    curntg: usize,
    /// Which `phrtsd` to use for [`Randlib::set_seed_from_phrase`].
    phrtsd: Phrtsd,
}

impl Default for Randlib {
    fn default() -> Self {
        Self::new()
    }
}

impl Randlib {
    /// A generator in the state `ignlgi` would auto-initialise itself to.
    ///
    /// `ignlgi` calls `setall(1234567890, 123456789)` on first use if nothing
    /// has seeded it yet, so an unseeded [`Randlib`] must start there too.
    #[must_use]
    pub fn new() -> Self {
        Self::with_phrtsd(Phrtsd::default())
    }

    /// As [`Randlib::new`], selecting which `phrtsd` variant is in force.
    #[must_use]
    pub fn with_phrtsd(phrtsd: Phrtsd) -> Self {
        let mut rng = Self {
            ig: [(0, 0); NUMG],
            lg: [(0, 0); NUMG],
            cg: [(0, 0); NUMG],
            curntg: 1,
            phrtsd,
        };
        rng.setall(1234567890, 123456789);
        rng
    }

    /// `salfph`: seed every generator from a phrase.
    ///
    /// This is what `Math::Random::random_set_seed_from_phrase` calls.
    pub fn set_seed_from_phrase(&mut self, phrase: &[u8]) {
        let (seed1, seed2) = self.phrtsd_seeds(phrase);
        self.setall(seed1, seed2);
    }

    /// `phrtsd`: derive two seeds from a phrase without touching generator state.
    ///
    /// Equivalent to `Math::Random::random_seed_from_phrase`.
    #[must_use]
    pub fn phrtsd_seeds(&self, phrase: &[u8]) -> (i64, i64) {
        match self.phrtsd {
            Phrtsd::New => phrtsd_new(phrase),
            Phrtsd::Orig => phrtsd_orig(phrase),
        }
    }

    /// `setall`: set generator 1 to the given seeds and derive the other 31.
    ///
    /// The C version flips the current generator inside the loop and restores
    /// it afterwards; since the saved value is always 1 here, the restore is
    /// written directly.
    pub fn setall(&mut self, iseed1: i64, iseed2: i64) {
        let ocgn = self.curntg;

        self.ig[0] = (iseed1, iseed2);
        self.curntg = 1;
        self.initgn_initial();

        for g in 2..=NUMG {
            self.ig[g - 1] = (
                mltmod(A1VW, self.ig[g - 2].0, M1),
                mltmod(A2VW, self.ig[g - 2].1, M2),
            );
            self.curntg = g;
            self.initgn_initial();
        }

        self.curntg = ocgn;
    }

    /// `initgn(-1)`: reset the current generator to its initial seed.
    ///
    /// Only the `isdtyp == -1` case is reachable from plethora, so the block
    /// and antithetic paths are not carried over.
    fn initgn_initial(&mut self) {
        let g = self.curntg - 1;
        self.lg[g] = self.ig[g];
        self.cg[g] = self.lg[g];
    }

    /// `ignlgi`: the next integer, uniform over `1..=2147483562`.
    pub fn ignlgi(&mut self) -> i64 {
        let g = self.curntg - 1;
        let (mut s1, mut s2) = self.cg[g];

        let k = s1 / 53668;
        s1 = A1 * (s1 - k * 53668) - k * 12211;
        if s1 < 0 {
            s1 += M1;
        }

        let k = s2 / 52774;
        s2 = A2 * (s2 - k * 52774) - k * 3791;
        if s2 < 0 {
            s2 += M2;
        }

        self.cg[g] = (s1, s2);

        let mut z = s1 - s2;
        if z < 1 {
            z += M1 - 1;
        }
        // `Xqanti` is all zeroes unless `setant` was called, which plethora
        // never does, so the antithetic branch is omitted.
        z
    }

    /// `ranf`: the next double, uniform over the open interval (0, 1).
    pub fn ranf(&mut self) -> f64 {
        self.ignlgi() as f64 * INV_M1
    }

    /// `gennor`: one normal deviate with the given mean and standard deviation.
    ///
    /// This is `Math::Random::random_normal(1, $av, $sd)` in scalar context.
    ///
    /// # Panics
    /// Panics if `sd` is negative, matching randlib's `SD < 0 in GENNOR - ABORT`.
    pub fn gennor(&mut self, av: f64, sd: f64) -> f64 {
        assert!(
            sd >= 0.0,
            "SD < 0 in GENNOR - ABORT (value of SD: {sd:16.6E})"
        );
        sd * self.snorm() + av
    }

    /// `snorm`: one standard normal deviate, algorithm FL with M = 5.
    ///
    /// Transcribed with the original statement labels preserved as an explicit
    /// state machine. The C source is a web of `goto`s whose control flow is
    /// not expressible as structured loops without duplicating blocks; keeping
    /// the labels means this can be diffed against randlib.c line by line, and
    /// it keeps the number and order of `ranf` calls identical, which is the
    /// only thing that actually has to match.
    #[allow(clippy::float_cmp)]
    pub fn snorm(&mut self) -> f64 {
        /// `a[32]`, `d[31]`, `t[31]` and `h[31]` from the Ahrens-Dieter paper.
        const A: [f64; 32] = [
            0.0,
            0.03917608550309,
            0.07841241273311,
            0.11776987457909,
            0.15731068461017,
            0.19709908429430,
            0.23720210932878,
            0.27769043982157,
            0.31863936396437,
            0.36012989178957,
            0.40225006532172,
            0.44509652498551,
            0.48877641111466,
            0.53340970624127,
            0.57913216225555,
            0.62609901234641,
            0.67448975019607,
            0.72451438349236,
            0.77642176114792,
            0.83051087820539,
            0.88714655901887,
            0.94678175630104,
            1.00999016924958,
            1.07751556704027,
            1.15034938037600,
            1.22985875921658,
            1.31801089730353,
            1.41779713799625,
            1.53412054435253,
            1.67593972277344,
            1.86273186742164,
            2.15387469406144,
        ];
        const D: [f64; 31] = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.26368432217502,
            0.24250845238097,
            0.22556744380930,
            0.21163416577204,
            0.19992426749317,
            0.18991075842246,
            0.18122518100691,
            0.17360140038056,
            0.16684190866667,
            0.16079672918053,
            0.15534971747692,
            0.15040938382813,
            0.14590257684509,
            0.14177003276856,
            0.13796317369537,
            0.13444176150074,
            0.13117215026483,
            0.12812596512583,
            0.12527909006226,
            0.12261088288608,
            0.12010355965651,
            0.11774170701949,
            0.11551189226063,
            0.11340234879117,
            0.11140272044119,
            0.10950385201710,
        ];
        const T: [f64; 31] = [
            7.6738283767E-4,
            2.30687039764E-3,
            3.86061844387E-3,
            5.43845406707E-3,
            7.05069876857E-3,
            8.70839582019E-3,
            1.042356984914E-2,
            1.220953194966E-2,
            1.408124734637E-2,
            1.605578804548E-2,
            1.815290075142E-2,
            2.039573175398E-2,
            2.281176732513E-2,
            2.543407332319E-2,
            2.830295595118E-2,
            3.146822492920E-2,
            3.499233438388E-2,
            3.895482964836E-2,
            4.345878381672E-2,
            4.864034918076E-2,
            5.468333844273E-2,
            6.184222395816E-2,
            7.047982761667E-2,
            8.113194985866E-2,
            9.462443534514E-2,
            0.11230007889456,
            0.13649799954975,
            0.17168856004707,
            0.22762405488269,
            0.33049802776911,
            0.58470309390507,
        ];
        const H: [f64; 31] = [
            3.920617164634E-2,
            3.932704963665E-2,
            3.950999486086E-2,
            3.975702679515E-2,
            4.007092772490E-2,
            4.045532602655E-2,
            4.091480886081E-2,
            4.145507115859E-2,
            4.208311051344E-2,
            4.280748137995E-2,
            4.363862733472E-2,
            4.458931789605E-2,
            4.567522779560E-2,
            4.691571371696E-2,
            4.833486978119E-2,
            4.996298427702E-2,
            5.183858644724E-2,
            5.401138183398E-2,
            5.654656186515E-2,
            5.953130423884E-2,
            6.308488965373E-2,
            6.737503494905E-2,
            7.264543556657E-2,
            7.926471414968E-2,
            8.781922325338E-2,
            9.930398323927E-2,
            0.11555994154118,
            0.14043438342816,
            0.18361418337460,
            0.27900163464163,
            0.70104742502766,
        ];

        /// The statement labels of algorithm FL, kept verbatim so this can be
        /// diffed against randlib.c. Fallthrough in the C source (S70 into S80,
        /// S110 into S120, S150 into S160) becomes an explicit transition.
        enum St {
            S40,
            S50,
            S60,
            S70,
            S80,
            S110,
            S120,
            S140,
            S150,
            S160,
        }

        let mut u = self.ranf();
        let mut s = 0.0_f64;
        if u > 0.5 {
            s = 1.0;
        }
        u += u - s;
        u = 32.0 * u;
        let mut i = u as i64;
        if i == 32 {
            i = 31;
        }

        let mut aa;
        let mut ustar = 0.0_f64;
        let mut w = 0.0_f64;
        let mut tt = 0.0_f64;

        let mut state = if i == 0 {
            // S100: start tail.
            i = 6;
            aa = A[31];
            St::S120
        } else {
            // Start center.
            ustar = u - i as f64;
            aa = A[(i - 1) as usize];
            St::S40
        };

        loop {
            state = match state {
                St::S40 => {
                    if ustar <= T[(i - 1) as usize] {
                        St::S60
                    } else {
                        w = (ustar - T[(i - 1) as usize]) * H[(i - 1) as usize];
                        St::S50
                    }
                }
                // Exit, reached from both the center and the tail.
                St::S50 => {
                    let y = aa + w;
                    return if s == 1.0 { -y } else { y };
                }
                // Center continued.
                St::S60 => {
                    u = self.ranf();
                    w = u * (A[i as usize] - aa);
                    tt = (0.5 * w + aa) * w;
                    St::S80
                }
                St::S70 => {
                    tt = u;
                    ustar = self.ranf();
                    St::S80
                }
                St::S80 => {
                    if ustar > tt {
                        St::S50
                    } else {
                        u = self.ranf();
                        if ustar >= u {
                            St::S70
                        } else {
                            ustar = self.ranf();
                            St::S40
                        }
                    }
                }
                // Tail.
                St::S110 => {
                    aa += D[(i - 1) as usize];
                    i += 1;
                    St::S120
                }
                St::S120 => {
                    u += u;
                    if u < 1.0 {
                        St::S110
                    } else {
                        u -= 1.0;
                        St::S140
                    }
                }
                St::S140 => {
                    w = u * D[(i - 1) as usize];
                    tt = (0.5 * w + aa) * w;
                    St::S160
                }
                St::S150 => {
                    tt = u;
                    St::S160
                }
                St::S160 => {
                    ustar = self.ranf();
                    if ustar > tt {
                        St::S50
                    } else {
                        u = self.ranf();
                        if ustar >= u {
                            St::S150
                        } else {
                            u = self.ranf();
                            St::S140
                        }
                    }
                }
            };
        }
    }
}

/// `mltmod`: `(a * s) mod m`, by L'Ecuyer and Cote's decomposition.
///
/// The decomposition exists so the intermediates fit in a 32-bit `long`. It is
/// kept rather than replaced by `i128` arithmetic because `setall` chains 31
/// calls through it, and any difference in an intermediate would propagate.
///
/// # Panics
/// Panics on arguments outside `0 < a < m`, `0 < s < m`, as randlib aborts.
#[must_use]
fn mltmod(a: i64, s: i64, m: i64) -> i64 {
    /// `2**((b-2)/2)` for a 32-bit `long`.
    const H: i64 = 32768;

    assert!(
        a > 0 && a < m && s > 0 && s < m,
        "a, m, s out of order in mltmod - ABORT! a = {a} s = {s} m = {m}"
    );

    let a0;
    let mut p;

    if a < H {
        a0 = a;
        p = 0;
    } else {
        let mut a1 = a / H;
        a0 = a - H * a1;
        let qh = m / H;
        let rh = m - H * qh;

        if a1 >= H {
            // A2 = 1
            a1 -= H;
            let k = s / qh;
            p = H * (s - k * qh) - k * rh;
            while p < 0 {
                p += m;
            }
        } else {
            p = 0;
        }

        // p = (A2 * s * H) mod m
        if a1 != 0 {
            let q = m / a1;
            let k = s / q;
            p -= k * (m - a1 * q);
            if p > 0 {
                p -= m;
            }
            p += a1 * (s - k * q);
            while p < 0 {
                p += m;
            }
        }

        // p = ((A2 * H + A1) * s) mod m
        let k = p / qh;
        p = H * (p - k * qh) - k * rh;
        while p < 0 {
            p += m;
        }
    }

    // p = ((A2 * H + A1) * H * s) mod m
    if a0 != 0 {
        let q = m / a0;
        let k = s / q;
        p -= k * (m - a0 * q);
        if p > 0 {
            p -= m;
        }
        p += a0 * (s - k * q);
        while p < 0 {
            p += m;
        }
    }

    p
}

/// `lennob`: length ignoring trailing blanks, but not other whitespace.
#[must_use]
fn lennob(phrase: &[u8]) -> usize {
    let mut i_nb: isize = -1;
    for (i, &c) in phrase.iter().enumerate() {
        if c != b' ' {
            i_nb = i as isize;
        }
    }
    (i_nb + 1) as usize
}

/// The `#else` branch of `phrtsd`.
///
/// Note the loop bound: `for(i=0; i<(lphr-1); i++)` stops one short, so the
/// last character of the phrase never reaches the seeds. Every 32-character
/// MD5 that plethora feeds in is therefore hashed as its first 31 characters.
#[must_use]
fn phrtsd_new(phrase: &[u8]) -> (i64, i64) {
    let mut seed1: i64 = 1234567890;
    let mut seed2: i64 = 123456789;

    let lphr = lennob(phrase);
    if lphr < 1 {
        return (seed1, seed2);
    }

    for i in 0..lphr.saturating_sub(1) {
        // `ichr` is a C `char` widened to `long`, and plain `char` is signed on
        // x86-64 (where plethora ran) and on Apple ARM64 (where the oracle
        // runs), so bytes >= 0x80 arrive negative. The golden vectors carry a
        // UTF-8 sequence precisely to pin this down: with an unsigned widening
        // the seeds for "\xc3\xa9" come out 748823347 instead of 714741811.
        let ichr = i64::from(phrase[i] as i8);
        let j = i % 8;
        seed1 = (seed1 + VALUES_NEW[j] * ichr) % TWOP30;
        seed2 = (seed2 + VALUES_NEW[7 - j] * ichr) % TWOP30;
    }

    (seed1, seed2)
}

/// The `#ifdef PHRTSD_ORIG` branch of `phrtsd`.
#[must_use]
fn phrtsd_orig(phrase: &[u8]) -> (i64, i64) {
    let mut seed1: i64 = 1234567890;
    let mut seed2: i64 = 123456789;

    let lphr = lennob(phrase);
    if lphr < 1 {
        return (seed1, seed2);
    }

    for &c in &phrase[..lphr] {
        // Find the character, then step one past it to match Fortran's 1-based
        // index ("JJV added ix++").
        let mut ix = TABLE.iter().position(|&t| t == c).unwrap_or(TABLE.len());
        ix += 1;

        // In C this reads `table[ix]` after the increment. For the final glyph
        // (the space) that lands on the terminating NUL and resets ix to 0; for
        // a character absent from the table it reads one past the NUL, which is
        // out of bounds and formally undefined. Both are treated as "not a
        // glyph" here, which is what the defined case does and what every
        // compiler observed in practice does for the other.
        if ix >= TABLE.len() {
            ix = 0;
        }

        let mut ichr = (ix % 64) as i64;
        if ichr == 0 {
            ichr = 63;
        }

        let mut values = [0_i64; 5];
        for j in 1..=5_i64 {
            let mut v = ichr - j;
            if v < 1 {
                v += 63;
            }
            values[(j - 1) as usize] = v;
        }

        for j in 1..=5_usize {
            seed1 = (seed1 + SHIFT[j - 1] * values[j - 1]) % TWOP30;
            seed2 = (seed2 + SHIFT[j - 1] * values[5 - j]) % TWOP30;
        }
    }

    (seed1, seed2)
}
