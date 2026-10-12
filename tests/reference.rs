//! Integration tests on real, geometry-optimized molecules.
//!
//! Dataset 1 — `tests/5dendralene/`: 16 conformers of [5]dendralene
//! (C10H12, 22 atoms each, identical atom ordering).
//!
//! Dataset 2 — `tests/golden/`: optimized small molecules — two ethane
//! rotamers (`ethane_*`, C2H6), a rotated-frame pair (`mirror_*`, C3H6),
//! a rotated/translated water (`water_*`), and two ALA dipeptide
//! snapshots (`mini_*`, PDB, 8 atoms).
//!
//! Reference values were recorded from Python `calculate_rmsd.py`
//! (numpy 2.4 / scipy 1.17) and must match to 1e-9 — EXCEPT the cases
//! marked `SENSITIVE` below: `ethane` and `mirror` are symmetric tops
//! (degenerate inertia eigenvalues, e.g. ethane
//! `[6.33, 25.41, 25.41]`), so `inertia-hungarian` depends on LAPACK
//! eigenvector-basis details, and `distance` depends on tie-breaking
//! among equal centroid distances (`np.argsort` quicksort is unstable).
//! For those, the tests only assert the methods run deterministically
//! and return sane values; exact cross-implementation parity is not
//! achievable even between LAPACK builds. [5]dendralene has three
//! distinct moments (`[240.7, 557.1, 615.2]`) and is fully deterministic.
//!
//! NOTE: `brute` reorder is only exercised on the 8-atom ethane pair —
//! on 22-atom molecules the per-element factorial search never finishes.

use rmsd::io::{pdb, xyz, Molecule};
use rmsd::pipeline::{Filter, Options, Reflections, Rotation};
use rmsd::reorder::ReorderMethod;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn fixture(sub: &str) -> std::path::PathBuf {
    // Manifest-relative so tests work regardless of runner cwd.
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join(sub)
}

fn load_dend(n: u8) -> Molecule {
    let text = std::fs::read_to_string(fixture(&format!("5dendralene/5dendralene_{n:02}.xyz")))
        .unwrap();
    let m = xyz::parse_xyz(&text).unwrap();
    assert_eq!(m.len(), 22, "dendralene conformer must have 22 atoms");
    m
}

fn load_xyz(name: &str) -> Molecule {
    let text = std::fs::read_to_string(fixture(&format!("golden/{name}"))).unwrap();
    xyz::parse_xyz(&text).unwrap()
}

fn load_pdb(name: &str) -> Molecule {
    let text = std::fs::read_to_string(fixture(&format!("golden/{name}"))).unwrap();
    pdb::parse_pdb(&text, false).unwrap()
}

fn opts(
    rotation: Rotation,
    reorder: ReorderMethod,
    do_reorder: bool,
    reflections: Reflections,
    filter: Filter,
) -> Options {
    Options {
        rotation,
        reorder,
        do_reorder,
        reflections,
        filter,
        print_rotated: false,
        print_only_rmsd_atoms: false,
    }
}

fn plain(rotation: Rotation) -> Options {
    opts(
        rotation,
        ReorderMethod::None,
        false,
        Reflections::None,
        Filter::All,
    )
}

fn reordered(method: ReorderMethod) -> Options {
    opts(
        Rotation::Kabsch,
        method,
        true,
        Reflections::None,
        Filter::All,
    )
}

fn rmsd(a: &Molecule, b: &Molecule, o: &Options) -> f64 {
    rmsd::pipeline::compute_rmsd(a, b, o, ("a", "b"))
        .unwrap()
        .rmsd
}

/// Smoke assertion for LAPACK/tie-SENSITIVE methods: deterministic,
/// finite, non-negative. See module docs for why no golden value applies.
fn assert_sane(a: &Molecule, b: &Molecule, o: &Options) {
    let r1 = rmsd(a, b, o);
    let r2 = rmsd(a, b, o);
    assert!(r1.is_finite() && r1 >= 0.0, "rmsd={r1}");
    assert_eq!(r1, r2, "method must be deterministic");
}

// ---------------------------------------------------------------------------
// [5]dendralene: identity + pairwise Kabsch vs Python references
// ---------------------------------------------------------------------------

#[test]
fn dendralene_identity_is_zero() {
    let m = load_dend(1);
    assert!(rmsd(&m, &m, &plain(Rotation::Kabsch)) < 1e-9);
    assert!(rmsd(&m, &m, &plain(Rotation::Quaternion)) < 1e-9);
    assert_eq!(rmsd(&m, &m, &plain(Rotation::None)), 0.0);
}

#[test]
fn dendralene_all_conformers_load() {
    // Every conformer: 22 atoms, C10H12, same ordering as #01.
    let ref_atoms = load_dend(1).atoms;
    for n in 1..=16u8 {
        let m = load_dend(n);
        assert_eq!(m.atoms, ref_atoms, "ordering differs in conformer {n:02}");
        assert_eq!(m.atoms.iter().filter(|&&z| z == 6).count(), 10);
        assert_eq!(m.atoms.iter().filter(|&&z| z == 1).count(), 12);
    }
}

#[test]
fn dendralene_kabsch_matches_python() {
    // (a, b, python `calculate_rmsd.py a b`)
    let cases = [
        (1u8, 2u8, 1.4907239185611663),
        (1, 8, 1.8808021830696275),
        (5, 12, 1.9929059496143893),
        (3, 16, 1.5164113506765755),
    ];
    for (a, b, expected) in cases {
        let p = load_dend(a);
        let q = load_dend(b);
        approx::assert_abs_diff_eq!(
            rmsd(&p, &q, &plain(Rotation::Kabsch)),
            expected,
            epsilon = 1e-9
        );
    }
}

#[test]
fn dendralene_quaternion_matches_kabsch() {
    for (a, b) in [(1u8, 2u8), (1, 8), (5, 12), (3, 16)] {
        let p = load_dend(a);
        let q = load_dend(b);
        let k = rmsd(&p, &q, &plain(Rotation::Kabsch));
        let v = rmsd(&p, &q, &plain(Rotation::Quaternion));
        approx::assert_abs_diff_eq!(k, v, epsilon = 1e-9);
    }
}

#[test]
fn dendralene_no_rotation_matches_python() {
    let cases = [
        (1u8, 2u8, 4.908887499676395),
        (1, 8, 4.842532795100448),
        (5, 12, 3.150807104735286),
        (3, 16, 4.90925043865833),
    ];
    for (a, b, expected) in cases {
        let p = load_dend(a);
        let q = load_dend(b);
        approx::assert_abs_diff_eq!(
            rmsd(&p, &q, &plain(Rotation::None)),
            expected,
            epsilon = 1e-9
        );
    }
}

// ---------------------------------------------------------------------------
// [5]dendralene: reorder methods vs Python references (deterministic here)
// ---------------------------------------------------------------------------

#[test]
fn dendralene_hungarian_matches_python() {
    let cases = [
        (1u8, 2u8, 1.534953769263205),
        (1, 8, 1.3248090151227827),
        (5, 12, 1.2915218707305982),
        (3, 16, 1.7905093191924768),
    ];
    for (a, b, expected) in cases {
        let p = load_dend(a);
        let q = load_dend(b);
        approx::assert_abs_diff_eq!(
            rmsd(&p, &q, &reordered(ReorderMethod::Hungarian)),
            expected,
            epsilon = 1e-9
        );
    }
}

#[test]
fn dendralene_inertia_hungarian_matches_python() {
    let cases = [
        (1u8, 2u8, 2.227061331783059),
        (1, 8, 1.107871112150934),
        (5, 12, 2.0960900474410584),
        (3, 16, 2.3031953569453925),
    ];
    for (a, b, expected) in cases {
        let p = load_dend(a);
        let q = load_dend(b);
        approx::assert_abs_diff_eq!(
            rmsd(&p, &q, &reordered(ReorderMethod::InertiaHungarian)),
            expected,
            epsilon = 1e-9
        );
    }
}

#[test]
fn dendralene_distance_matches_python() {
    let cases = [
        (1u8, 2u8, 2.827650466383082),
        (1, 8, 2.809242189168922),
        (5, 12, 2.906316788688958),
        (3, 16, 3.2443000465628535),
    ];
    for (a, b, expected) in cases {
        let p = load_dend(a);
        let q = load_dend(b);
        approx::assert_abs_diff_eq!(
            rmsd(&p, &q, &reordered(ReorderMethod::Distance)),
            expected,
            epsilon = 1e-9
        );
    }
}

// ---------------------------------------------------------------------------
// [5]dendralene: reflections, filtering, mathematical properties
// ---------------------------------------------------------------------------

#[test]
fn dendralene_reflections_match_python() {
    // Mirrors CLI `--use-reflections`: the default inertia-hungarian
    // reorder runs *inside* the reflection search (do_reorder=false).
    let o = opts(
        Rotation::Kabsch,
        ReorderMethod::InertiaHungarian,
        false,
        Reflections::All,
        Filter::All,
    );
    let cases = [
        (1u8, 2u8, 0.8926868464424041),
        (1, 8, 1.107871112150934),
        (5, 12, 0.8136470352421988),
        (3, 16, 0.6309375797268141),
    ];
    for (a, b, expected) in cases {
        let p = load_dend(a);
        let q = load_dend(b);
        approx::assert_abs_diff_eq!(rmsd(&p, &q, &o), expected, epsilon = 1e-9);
    }
}

#[test]
fn dendralene_reflections_keep_stereo_never_below_full_search() {
    // keep-stereo explores a subset of the 48 trials, so its RMSD must be
    // greater than or equal to the full reflection search.
    let full = opts(
        Rotation::Kabsch,
        ReorderMethod::None,
        false,
        Reflections::All,
        Filter::All,
    );
    let keep = opts(
        Rotation::Kabsch,
        ReorderMethod::None,
        false,
        Reflections::KeepStereo,
        Filter::All,
    );
    for (a, b) in [(1u8, 2u8), (5, 12)] {
        let p = load_dend(a);
        let q = load_dend(b);
        assert!(rmsd(&p, &q, &keep) + 1e-12 >= rmsd(&p, &q, &full));
    }
}

#[test]
fn dendralene_ignore_hydrogen_matches_python() {
    let o = opts(
        Rotation::Kabsch,
        ReorderMethod::None,
        false,
        Reflections::None,
        Filter::NoHydrogen,
    );
    let cases = [
        (1u8, 2u8, 1.0001689591463656),
        (1, 8, 1.3565110444225699),
        (5, 12, 1.4304095082330701),
        (3, 16, 0.9823157554425462),
    ];
    for (a, b, expected) in cases {
        let p = load_dend(a);
        let q = load_dend(b);
        approx::assert_abs_diff_eq!(rmsd(&p, &q, &o), expected, epsilon = 1e-9);
    }
}

#[test]
fn dendralene_rmsd_is_symmetric() {
    // RMSD(a,b) == RMSD(b,a) for alignment-based methods.
    for (a, b) in [(1u8, 2u8), (5, 12)] {
        let p = load_dend(a);
        let q = load_dend(b);
        for o in [
            plain(Rotation::Kabsch),
            plain(Rotation::Quaternion),
            reordered(ReorderMethod::Hungarian),
        ] {
            approx::assert_abs_diff_eq!(rmsd(&p, &q, &o), rmsd(&q, &p, &o), epsilon = 1e-9);
        }
    }
}

#[test]
fn dendralene_kabsch_never_above_unrotated() {
    // Optimal rotation can only lower (or keep) the RMSD.
    for (a, b) in [(1u8, 2u8), (1, 8), (5, 12), (3, 16)] {
        let p = load_dend(a);
        let q = load_dend(b);
        assert!(
            rmsd(&p, &q, &plain(Rotation::Kabsch)) <= rmsd(&p, &q, &plain(Rotation::None)) + 1e-12
        );
    }
}

#[test]
fn dendralene_shuffled_order_recovered_by_hungarian() {
    // Same geometry with permuted atom order (within each element) must
    // realign to ~0 via Hungarian reorder. Pure property test, no reference.
    let m = load_dend(7);
    // Deterministic shuffle: reverse each element block.
    let mut perm: Vec<usize> = (0..m.len()).collect();
    for z in [6u8, 1u8] {
        let block: Vec<usize> = m
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, &a)| a == z)
            .map(|(i, _)| i)
            .rev()
            .collect();
        let mut it = block.into_iter();
        for (i, a) in m.atoms.iter().enumerate() {
            if *a == z {
                perm[i] = it.next().unwrap();
            }
        }
    }
    let shuffled = Molecule {
        atoms: perm.iter().map(|&i| m.atoms[i]).collect(),
        coords: perm.iter().map(|&i| m.coords[i]).collect(),
    };
    let unreordered = rmsd(&m, &shuffled, &plain(Rotation::Kabsch));
    assert!(unreordered > 1e-6, "shuffle must disturb naive RMSD");
    let fixed = rmsd(&m, &shuffled, &reordered(ReorderMethod::Hungarian));
    assert!(fixed < 1e-9, "hungarian must recover shuffled order");
}

// ---------------------------------------------------------------------------
// Golden: ethane rotamers (C2H6, 8 atoms)
// ---------------------------------------------------------------------------

#[test]
fn golden_ethane_kabsch_matches_python() {
    // Two optimized rotamers: heavy atoms nearly identical (see noH below),
    // hydrogens differ -> Kabsch 0.947.
    let p = load_xyz("ethane_a.xyz");
    let q = load_xyz("ethane_b.xyz");
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &plain(Rotation::Kabsch)),
        0.9466884271457764,
        epsilon = 1e-9
    );
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &plain(Rotation::Quaternion)),
        0.9466884271457765,
        epsilon = 1e-9
    );
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &plain(Rotation::None)),
        2.0726029480639867,
        epsilon = 1e-9
    );
}

#[test]
fn golden_ethane_hungarian_matches_python() {
    let p = load_xyz("ethane_a.xyz");
    let q = load_xyz("ethane_b.xyz");
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &reordered(ReorderMethod::Hungarian)),
        0.4587744955857676,
        epsilon = 1e-9
    );
}

#[test]
fn golden_ethane_heavy_atoms_nearly_identical() {
    // Ignoring hydrogens, the rotamers coincide to ~9e-4.
    let p = load_xyz("ethane_a.xyz");
    let q = load_xyz("ethane_b.xyz");
    let o = opts(
        Rotation::Kabsch,
        ReorderMethod::None,
        false,
        Reflections::None,
        Filter::NoHydrogen,
    );
    approx::assert_abs_diff_eq!(rmsd(&p, &q, &o), 0.0009205572465013836, epsilon = 1e-9);
}

#[test]
fn golden_ethane_reflections_match_python() {
    let p = load_xyz("ethane_a.xyz");
    let q = load_xyz("ethane_b.xyz");
    // As above: default reorder method active inside the search.
    for refl in [Reflections::All, Reflections::KeepStereo] {
        let o = opts(
            Rotation::Kabsch,
            ReorderMethod::InertiaHungarian,
            false,
            refl,
            Filter::All,
        );
        approx::assert_abs_diff_eq!(rmsd(&p, &q, &o), 0.45877449558576744, epsilon = 1e-9);
    }
}

#[test]
fn golden_ethane_brute_matches_python() {
    // Small enough (C2/H6) for factorial search.
    let p = load_xyz("ethane_a.xyz");
    let q = load_xyz("ethane_b.xyz");
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &reordered(ReorderMethod::Brute)),
        0.8953036049080758,
        epsilon = 1e-9
    );
}

#[test]
fn golden_ethane_symmetric_top_methods_are_sane() {
    // SENSITIVE (see module docs): ethane is a symmetric top, so
    // inertia-hungarian and distance tie-breaking are LAPACK/sort
    // dependent. Only determinism + sanity are asserted.
    let p = load_xyz("ethane_a.xyz");
    let q = load_xyz("ethane_b.xyz");
    assert_sane(&p, &q, &reordered(ReorderMethod::InertiaHungarian));
    assert_sane(&p, &q, &reordered(ReorderMethod::Distance));
}

// ---------------------------------------------------------------------------
// Golden: rotated-frame pair (C3H6, 9 atoms) — rotation must be recovered
// ---------------------------------------------------------------------------

#[test]
fn golden_mirror_rotation_recovered_without_reorder() {
    // Same optimized geometry in a different frame: Kabsch -> ~0.
    let p = load_xyz("mirror_a.xyz");
    let q = load_xyz("mirror_b.xyz");
    assert!(rmsd(&p, &q, &plain(Rotation::Kabsch)) < 1e-6);
    assert!(rmsd(&p, &q, &plain(Rotation::Quaternion)) < 1e-6);
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &plain(Rotation::None)),
        1.9219280059605537,
        epsilon = 1e-9
    );
}

#[test]
fn golden_mirror_hungarian_and_distance_match_python() {
    let p = load_xyz("mirror_a.xyz");
    let q = load_xyz("mirror_b.xyz");
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &reordered(ReorderMethod::Hungarian)),
        0.849853792637902,
        epsilon = 1e-9
    );
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &reordered(ReorderMethod::Distance)),
        1.4328894040512608,
        epsilon = 1e-9
    );
}

#[test]
fn golden_mirror_inertia_hungarian_is_sane() {
    // SENSITIVE (see module docs): degenerate inertia pair
    // [24.94, 24.94, 39.75] makes the inertia frame LAPACK-dependent.
    let p = load_xyz("mirror_a.xyz");
    let q = load_xyz("mirror_b.xyz");
    assert_sane(&p, &q, &reordered(ReorderMethod::InertiaHungarian));
}

#[test]
fn synthetic_mirror_image_needs_reflections() {
    // True mirror pair (y negated in memory): rotation alone cannot align,
    // reflection search recovers ~0. No fixture file needed.
    let m = load_dend(1);
    let mirrored = Molecule {
        atoms: m.atoms.clone(),
        coords: m.coords.iter().map(|v| [v[0], -v[1], v[2]]).collect(),
    };
    assert!(rmsd(&m, &mirrored, &plain(Rotation::Kabsch)) > 1e-3);
    let o = opts(
        Rotation::Kabsch,
        ReorderMethod::None,
        false,
        Reflections::All,
        Filter::All,
    );
    assert!(rmsd(&m, &mirrored, &o) < 1e-9);
}

// ---------------------------------------------------------------------------
// Golden: water + PDB snapshots
// ---------------------------------------------------------------------------

#[test]
fn golden_water_matches_python() {
    let p = load_xyz("water_a.xyz");
    let q = load_xyz("water_b.xyz");
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &plain(Rotation::Kabsch)),
        1.6825252341787082e-09,
        epsilon = 1e-6
    );
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &plain(Rotation::None)),
        0.6770102577292689,
        epsilon = 1e-9
    );
}

#[test]
fn golden_pdb_matches_python() {
    // OpenBabel ALA snapshots (ATOM+HETATM mix, element column present).
    let p = load_pdb("mini_a.pdb");
    let q = load_pdb("mini_b.pdb");
    assert_eq!(p.len(), 8);
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &plain(Rotation::Kabsch)),
        0.012437937649508802,
        epsilon = 1e-9
    );
    approx::assert_abs_diff_eq!(
        rmsd(&p, &q, &plain(Rotation::None)),
        0.21160547015613776,
        epsilon = 1e-9
    );
    let o = opts(
        Rotation::Kabsch,
        ReorderMethod::InertiaHungarian,
        false,
        Reflections::All,
        Filter::All,
    );
    approx::assert_abs_diff_eq!(rmsd(&p, &q, &o), 0.01238223045975489, epsilon = 1e-9);
}
