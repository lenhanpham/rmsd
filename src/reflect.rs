//! Reflection / axis-swap search.
//!
//! Port of `check_reflections`, `AXIS_SWAPS`, `AXIS_REFLECTIONS`
//! in `calculate_rmsd.py`.

use crate::geom::Coord;

/// All 6 axis permutations.
pub const AXIS_SWAPS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 1, 0],
    [2, 0, 1],
];

/// All 8 sign reflections.
pub const AXIS_REFLECTIONS: [[f64; 3]; 8] = [
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [-1.0, 1.0, -1.0],
    [1.0, -1.0, -1.0],
    [-1.0, -1.0, -1.0],
];

fn centroid_of(coords: &[Coord]) -> Coord {
    let mut s = [0.0; 3];
    for p in coords {
        s[0] += p[0];
        s[1] += p[1];
        s[2] += p[2];
    }
    let n = coords.len() as f64;
    [s[0] / n, s[1] / n, s[2] / n]
}

/// Minimize RMSD over reflections/axis-swaps.
///
/// * `reorder`: optional `fn(p_atoms,q_atoms,p_coord,q_coord) -> q-view`.
/// * `rmsd_fn`: e.g. Kabsch / quaternion / plain RMSD.
/// * `keep_stereo`: skip enantiomers (`i*j == -1`), mirroring Python masks.
///
/// Returns `(min_rmsd, min_swap, min_reflection, min_review)`.
pub fn check_reflections(
    p_atoms: &[u8],
    q_atoms: &[u8],
    p_coord: &[Coord],
    q_coord: &[Coord],
    reorder: Option<fn(&[u8], &[u8], &[Coord], &[Coord]) -> Vec<usize>>,
    rmsd_fn: fn(&[Coord], &[Coord]) -> f64,
    keep_stereo: bool,
) -> (f64, [usize; 3], [f64; 3], Vec<usize>) {
    // Parity with Python's swap/reflection masks.
    const SWAP_MASK: [i32; 6] = [1, -1, -1, 1, -1, 1];
    const REFL_MASK: [i32; 8] = [1, -1, -1, -1, 1, 1, 1, -1];

    let mut min_rmsd = f64::INFINITY;
    let mut min_swap = AXIS_SWAPS[0];
    let mut min_reflection = AXIS_REFLECTIONS[0];
    let mut min_review: Vec<usize> = (0..p_atoms.len()).collect();

    for (swap, &si) in AXIS_SWAPS.iter().zip(SWAP_MASK.iter()) {
        for (refl, &rj) in AXIS_REFLECTIONS.iter().zip(REFL_MASK.iter()) {
            if keep_stereo && si * rj == -1 {
                continue;
            }
            // q[:, swap] * reflection, then recenter.
            let mut tmp: Vec<Coord> = q_coord
                .iter()
                .map(|q| {
                    [
                        q[swap[0]] * refl[0],
                        q[swap[1]] * refl[1],
                        q[swap[2]] * refl[2],
                    ]
                })
                .collect();
            let c = centroid_of(&tmp);
            for v in tmp.iter_mut() {
                v[0] -= c[0];
                v[1] -= c[1];
                v[2] -= c[2];
            }
            let mut review: Vec<usize> = (0..p_atoms.len()).collect();
            let trial: Vec<Coord>;
            if let Some(f) = reorder {
                review = f(p_atoms, q_atoms, p_coord, &tmp);
                trial = review.iter().map(|&i| tmp[i]).collect();
            } else {
                // No reorder: require atoms already in order (Python asserts).
                debug_assert_eq!(p_atoms, q_atoms);
                trial = tmp;
            }
            let r = rmsd_fn(p_coord, &trial);
            if r < min_rmsd {
                min_rmsd = r;
                min_swap = *swap;
                min_reflection = *refl;
                min_review = review;
            }
        }
    }
    (min_rmsd, min_swap, min_reflection, min_review)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::kabsch_rmsd;

    fn kabsch_plain(p: &[Coord], q: &[Coord]) -> f64 {
        kabsch_rmsd(p, q, None, false)
    }

    #[test]
    fn finds_mirror() {
        let p_atoms = vec![6, 6];
        let q_atoms = vec![6, 6];
        let p = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        // Mirror in y; without reflections RMSD > 0 after centering? With search -> 0.
        let q = vec![[1.0, 0.0, 0.0], [0.0, -1.0, 0.0]];
        let pc = {
            let c = centroid_of(&p);
            p.iter()
                .map(|v| [v[0] - c[0], v[1] - c[1], v[2] - c[2]])
                .collect::<Vec<_>>()
        };
        let qc = {
            let c = centroid_of(&q);
            q.iter()
                .map(|v| [v[0] - c[0], v[1] - c[1], v[2] - c[2]])
                .collect::<Vec<_>>()
        };
        let (r, _, _, _) =
            check_reflections(&p_atoms, &q_atoms, &pc, &qc, None, kabsch_plain, false);
        assert!(r < 1e-9, "rmsd={r}");
    }
}
