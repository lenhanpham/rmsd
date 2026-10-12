//! Atom reordering.
//!
//! Port of `reorder_hungarian` (+ `_hungarian`), `reorder_distance`,
//! `reorder_inertia_hungarian`, `reorder_brute` (+ `brute_permutation`,
//! `generate_permutations`) in `calculate_rmsd.py`.
//!
//! NOTE: `reorder_similarity` (QML/FCHL19) is intentionally omitted.
//! See README for supported methods.

use crate::geom::{kabsch_rmsd, Coord};
use crate::inertia::{center_of_mass, inertia_tensor, sorted_eigenvectors};
use std::collections::BTreeSet;

/// Supported reorder methods (no `qml` — see module docs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReorderMethod {
    None,
    Hungarian,
    InertiaHungarian,
    Distance,
    Brute,
}

impl ReorderMethod {
    pub fn name(&self) -> &'static str {
        match self {
            ReorderMethod::None => "none",
            ReorderMethod::Hungarian => "hungarian",
            ReorderMethod::InertiaHungarian => "inertia-hungarian",
            ReorderMethod::Distance => "distance",
            ReorderMethod::Brute => "brute",
        }
    }

    pub fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "none" => Some(ReorderMethod::None),
            "hungarian" => Some(ReorderMethod::Hungarian),
            "inertia-hungarian" => Some(ReorderMethod::InertiaHungarian),
            "distance" => Some(ReorderMethod::Distance),
            "brute" => Some(ReorderMethod::Brute),
            _ => None,
        }
    }
}

/// Dispatch helper used by pipeline/reflections.
pub fn reorder_dispatch(
    method: ReorderMethod,
    p_atoms: &[u8],
    q_atoms: &[u8],
    p_coord: &[Coord],
    q_coord: &[Coord],
) -> Option<Vec<usize>> {
    match method {
        ReorderMethod::None => None,
        ReorderMethod::Hungarian => Some(reorder_hungarian(p_atoms, q_atoms, p_coord, q_coord)),
        ReorderMethod::InertiaHungarian => Some(reorder_inertia_hungarian(
            p_atoms, q_atoms, p_coord, q_coord,
        )),
        ReorderMethod::Distance => Some(reorder_distance(p_atoms, q_atoms, p_coord, q_coord)),
        ReorderMethod::Brute => Some(reorder_brute(p_atoms, q_atoms, p_coord, q_coord)),
    }
}

fn unique_atoms(p_atoms: &[u8]) -> Vec<u8> {
    let mut set = BTreeSet::new();
    for &a in p_atoms {
        set.insert(a);
    }
    set.into_iter().collect()
}

fn indices_of(atoms: &[u8], z: u8) -> Vec<usize> {
    atoms
        .iter()
        .enumerate()
        .filter(|(_, &a)| a == z)
        .map(|(i, _)| i)
        .collect()
}

// --- Hungarian (Munkres) on a square cost matrix, row -> col ---

/// Solve linear-sum assignment, minimizing total cost.
/// Returns `col_of_row[row]`. Square matrices only (callers always
/// build per-element square blocks).
pub fn hungarian_assign(cost: &[Vec<f64>]) -> Vec<usize> {
    let n = cost.len();
    assert!(n > 0 && cost.iter().all(|r| r.len() == n));
    // Munkres (James Munkres 1957) O(n^3)-ish implementation.
    let mut c: Vec<Vec<f64>> = cost.to_vec();
    // Step 1: row reduction.
    for row in c.iter_mut() {
        let m = row.iter().cloned().fold(f64::INFINITY, f64::min);
        for v in row.iter_mut() {
            *v -= m;
        }
    }
    // Step 2: col reduction.
    for j in 0..n {
        let mut m = f64::INFINITY;
        for i in 0..n {
            m = m.min(c[i][j]);
        }
        for i in 0..n {
            c[i][j] -= m;
        }
    }
    let mut star = vec![vec![false; n]; n];
    let mut prime = vec![vec![false; n]; n];
    let mut row_covered = vec![false; n];
    let mut col_covered = vec![false; n];

    // Star a zero in each row/col if possible.
    for i in 0..n {
        for j in 0..n {
            if c[i][j].abs() < 1e-12 && !row_covered[i] && !col_covered[j] {
                star[i][j] = true;
                row_covered[i] = true;
                col_covered[j] = true;
            }
        }
    }
    row_covered.fill(false);
    col_covered.fill(false);
    for i in 0..n {
        for j in 0..n {
            if star[i][j] {
                col_covered[j] = true;
            }
        }
    }

    let cover_count = |cc: &[bool]| cc.iter().filter(|&&x| x).count();
    let find_uncovered_zero =
        |c: &[Vec<f64>], rc: &[bool], cc: &[bool], n: usize| -> Option<(usize, usize)> {
            for i in 0..n {
                if rc[i] {
                    continue;
                }
                for j in 0..n {
                    if cc[j] {
                        continue;
                    }
                    if c[i][j].abs() < 1e-10 {
                        return Some((i, j));
                    }
                }
            }
            None
        };

    loop {
        if cover_count(&col_covered) == n {
            break;
        }
        // Step 4-6 loop.
        loop {
            let Some((zi, zj)) = find_uncovered_zero(&c, &row_covered, &col_covered, n) else {
                // Step 6: adjust matrix.
                let mut m = f64::INFINITY;
                for i in 0..n {
                    if row_covered[i] {
                        continue;
                    }
                    for j in 0..n {
                        if col_covered[j] {
                            continue;
                        }
                        m = m.min(c[i][j]);
                    }
                }
                for i in 0..n {
                    for j in 0..n {
                        if row_covered[i] {
                            c[i][j] += m;
                        }
                        if !col_covered[j] {
                            c[i][j] -= m;
                        }
                    }
                }
                continue;
            };
            prime[zi][zj] = true;
            // Find starred zero in row.
            let star_col = (0..n).find(|&j| star[zi][j]);
            if let Some(sj) = star_col {
                row_covered[zi] = true;
                col_covered[sj] = false;
            } else {
                // Augmenting path.
                let mut path = vec![(zi, zj)];
                loop {
                    let (li, lj) = *path.last().unwrap();
                    // starred zero in column lj
                    let sr = (0..n).find(|&i| star[i][lj]);
                    let Some(r) = sr else { break };
                    path.push((r, lj));
                    // primed zero in row r
                    let pc = (0..n).find(|&j| prime[r][j]).unwrap();
                    path.push((r, pc));
                    let _ = li;
                }
                for (k, &(i, j)) in path.iter().enumerate() {
                    if k % 2 == 0 {
                        star[i][j] = true;
                    } else {
                        star[i][j] = false;
                    }
                }
                prime = vec![vec![false; n]; n];
                row_covered.fill(false);
                col_covered.fill(false);
                for i in 0..n {
                    for j in 0..n {
                        if star[i][j] {
                            col_covered[j] = true;
                        }
                    }
                }
                break;
            }
        }
    }

    let mut assign = vec![0usize; n];
    for i in 0..n {
        assign[i] = (0..n).find(|&j| star[i][j]).unwrap_or(0);
    }
    assign
}

fn euclidean_cost(a: &[Coord], b: &[Coord]) -> Vec<Vec<f64>> {
    a.iter()
        .map(|p| {
            b.iter()
                .map(|q| {
                    let dx = p[0] - q[0];
                    let dy = p[1] - q[1];
                    let dz = p[2] - q[2];
                    (dx * dx + dy * dy + dz * dz).sqrt()
                })
                .collect()
        })
        .collect()
}

/// Hungarian on same-element coordinate blocks (`_hungarian`).
fn hungarian_block(a: &[Coord], b: &[Coord]) -> Vec<usize> {
    let cost = euclidean_cost(a, b);
    hungarian_assign(&cost)
}

/// Hungarian-distance reorder (`reorder_hungarian`). Returns Q-view.
pub fn reorder_hungarian(
    p_atoms: &[u8],
    q_atoms: &[u8],
    p_coord: &[Coord],
    q_coord: &[Coord],
) -> Vec<usize> {
    let mut view = vec![usize::MAX; q_atoms.len()];
    for z in unique_atoms(p_atoms) {
        let pi = indices_of(p_atoms, z);
        let qi = indices_of(q_atoms, z);
        let a: Vec<Coord> = pi.iter().map(|&i| p_coord[i]).collect();
        let b: Vec<Coord> = qi.iter().map(|&i| q_coord[i]).collect();
        let cols = hungarian_block(&a, &b);
        for (k, &pidx) in pi.iter().enumerate() {
            view[pidx] = qi[cols[k]];
        }
    }
    view
}

/// Distance-from-centroid reorder (`reorder_distance`).
pub fn reorder_distance(
    p_atoms: &[u8],
    q_atoms: &[u8],
    p_coord: &[Coord],
    q_coord: &[Coord],
) -> Vec<usize> {
    let mut view = vec![0usize; q_atoms.len()];
    for z in unique_atoms(p_atoms) {
        let pi = indices_of(p_atoms, z);
        let qi = indices_of(q_atoms, z);
        let mut pa: Vec<(f64, usize)> = pi
            .iter()
            .map(|&i| {
                let p = p_coord[i];
                (p[0] * p[0] + p[1] * p[1] + p[2] * p[2], i)
            })
            .collect();
        let mut qb: Vec<(f64, usize)> = qi
            .iter()
            .map(|&i| {
                let p = q_coord[i];
                (p[0] * p[0] + p[1] * p[1] + p[2] * p[2], i)
            })
            .collect();
        pa.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        qb.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        // Project order of P onto Q: order_A sorted, order_B sorted.
        // Python: translator = argsort(reorder_A); view = reorder_B[translator].
        // Equivalent: k-th smallest in P maps to k-th smallest in Q.
        for (k, &(_, pidx)) in pa.iter().enumerate() {
            view[pidx] = qb[k].1;
        }
    }
    view
}

fn apply_mat_row(p: &Coord, m: &nalgebra::Matrix3<f64>) -> Coord {
    [
        p[0] * m[(0, 0)] + p[1] * m[(1, 0)] + p[2] * m[(2, 0)],
        p[0] * m[(0, 1)] + p[1] * m[(1, 1)] + p[2] * m[(2, 1)],
        p[0] * m[(0, 2)] + p[1] * m[(1, 2)] + p[2] * m[(2, 2)],
    ]
}

/// Inertia + Hungarian reorder (`reorder_inertia_hungarian`).
pub fn reorder_inertia_hungarian(
    p_atoms: &[u8],
    q_atoms: &[u8],
    p_coord: &[Coord],
    q_coord: &[Coord],
) -> Vec<usize> {
    let mut pc: Vec<Coord> = p_coord.to_vec();
    let mut qc: Vec<Coord> = q_coord.to_vec();
    let pcm = center_of_mass(p_atoms, &pc);
    let qcm = center_of_mass(q_atoms, &qc);
    for v in pc.iter_mut() {
        v[0] -= pcm[0];
        v[1] -= pcm[1];
        v[2] -= pcm[2];
    }
    for v in qc.iter_mut() {
        v[0] -= qcm[0];
        v[1] -= qcm[1];
        v[2] -= qcm[2];
    }
    let ip = inertia_tensor(p_atoms, &pc);
    let iq = inertia_tensor(q_atoms, &qc);
    let ep = sorted_eigenvectors(&ip);
    let eq = sorted_eigenvectors(&iq);
    let pa: Vec<Coord> = pc.iter().map(|v| apply_mat_row(v, &ep)).collect();

    let reflections = crate::reflect::AXIS_REFLECTIONS;
    let mut best_rmsd = f64::INFINITY;
    let mut best: Vec<usize> = (0..p_atoms.len()).collect();
    for mirror in reflections {
        // tmp_eigvec = eigvec_q with columns scaled by mirror (broadcast on last axis).
        let mut tmp_eig = eq;
        for j in 0..3 {
            tmp_eig[(0, j)] *= mirror[j];
            tmp_eig[(1, j)] *= mirror[j];
            tmp_eig[(2, j)] *= mirror[j];
        }
        let qa: Vec<Coord> = qc.iter().map(|v| apply_mat_row(v, &tmp_eig)).collect();
        let review = reorder_hungarian(p_atoms, q_atoms, &pa, &qa);
        let qa_re: Vec<Coord> = review.iter().map(|&i| qa[i]).collect();
        let r = kabsch_rmsd(&qa_re, &pa, None, false);
        if r < best_rmsd {
            best_rmsd = r;
            best = review;
        }
    }
    best
}

// --- Brute force ---

/// Heap's algorithm permutations, yielding index orders in place.
fn heaps<T: Clone, F: FnMut(&[T])>(a: &mut [T], mut f: F) {
    let n = a.len();
    f(a);
    let mut c = vec![0usize; n];
    let mut i = 0;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 {
                a.swap(0, i);
            } else {
                let k = c[i];
                a.swap(k, i);
            }
            f(a);
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
}

fn brute_permutation(a: &[Coord], b: &[Coord]) -> Vec<usize> {
    let n = a.len();
    let mut best_rmsd = f64::INFINITY;
    let mut best: Vec<usize> = (0..n).collect();
    let mut order: Vec<usize> = (0..n).collect();
    // Heap's yields the initial (identity) order first, matching Python.
    heaps(&mut order, |perm| {
        let re: Vec<Coord> = perm.iter().map(|&i| b[i]).collect();
        let r = kabsch_rmsd(a, &re, None, false);
        if r < best_rmsd {
            best_rmsd = r;
            best = perm.to_vec();
        }
    });
    best
}

/// Brute-force reorder (`reorder_brute`). Reference only.
pub fn reorder_brute(
    p_atoms: &[u8],
    q_atoms: &[u8],
    p_coord: &[Coord],
    q_coord: &[Coord],
) -> Vec<usize> {
    let mut view = vec![usize::MAX; q_atoms.len()];
    for z in unique_atoms(p_atoms) {
        let pi = indices_of(p_atoms, z);
        let qi = indices_of(q_atoms, z);
        let a: Vec<Coord> = pi.iter().map(|&i| p_coord[i]).collect();
        let b: Vec<Coord> = qi.iter().map(|&i| q_coord[i]).collect();
        let cols = brute_permutation(&a, &b);
        for (k, &pidx) in pi.iter().enumerate() {
            view[pidx] = qi[cols[k]];
        }
    }
    view
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hungarian_identity() {
        let cost = vec![vec![0.0, 5.0], vec![5.0, 0.0]];
        assert_eq!(hungarian_assign(&cost), vec![0, 1]);
    }

    #[test]
    fn hungarian_swap() {
        let cost = vec![vec![9.0, 1.0], vec![1.0, 9.0]];
        assert_eq!(hungarian_assign(&cost), vec![1, 0]);
    }

    #[test]
    fn reorder_hungarian_recovers_shuffle() {
        let p_atoms = vec![6, 6, 8];
        let q_atoms = vec![6, 6, 8];
        let p = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let q = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [5.0, 5.0, 5.0]];
        let v = reorder_hungarian(&p_atoms, &q_atoms, &p, &q);
        // p[0] best matches q[1], p[1] matches q[0]
        assert_eq!(v, vec![1, 0, 2]);
    }
}
