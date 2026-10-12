//! Geometry: RMSD, Kabsch (plain + weighted), quaternion rotation.
//!
//! Port of `rmsd`, `kabsch_rmsd`, `kabsch_rotate`, `kabsch_fit`, `kabsch`,
//! `kabsch_weighted`, `kabsch_weighted_fit`, `kabsch_weighted_rmsd`,
//! `quaternion_rmsd`, `quaternion_transform`, `_makeW`, `_makeQ`,
//! `quaternion_rotate`, `centroid` in `calculate_rmsd.py`.
//!
//! Row-vector convention is preserved: rotation `U` applies as
//! `p_rot[j] = sum_i p[i] * U[i][j]`, i.e. `p_col -> U^T p_col`.

use nalgebra::{Matrix3, Matrix4, SymmetricEigen, Vector3};

pub type Coord = [f64; 3];
pub type Quat = [f64; 4];

/// Mean position (`centroid`).
pub fn centroid(coords: &[Coord]) -> Coord {
    assert!(!coords.is_empty(), "centroid of empty set");
    let mut s = [0.0; 3];
    for p in coords {
        s[0] += p[0];
        s[1] += p[1];
        s[2] += p[2];
    }
    let n = coords.len() as f64;
    [s[0] / n, s[1] / n, s[2] / n]
}

/// Plain RMSD, no alignment (`rmsd`).
pub fn rmsd(p: &[Coord], q: &[Coord]) -> f64 {
    assert_eq!(p.len(), q.len());
    let mut acc = 0.0;
    for (a, b) in p.iter().zip(q.iter()) {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        acc += dx * dx + dy * dy + dz * dz;
    }
    (acc / p.len() as f64).sqrt()
}

fn covariance(p: &[Coord], q: &[Coord]) -> Matrix3<f64> {
    let mut c = Matrix3::zeros();
    for (a, b) in p.iter().zip(q.iter()) {
        for i in 0..3 {
            for j in 0..3 {
                c[(i, j)] += a[i] * b[j];
            }
        }
    }
    c
}

fn optimal_rotation(mut u: Matrix3<f64>, vt: Matrix3<f64>) -> Matrix3<f64> {
    // Enforce right-handed system: flip last column of U if det(U)*det(Vt) < 0.
    // Mirrors `d = (det(V)*det(W)) < 0` in `kabsch`.
    if u.determinant() * vt.determinant() < 0.0 {
        u[(0, 2)] = -u[(0, 2)];
        u[(1, 2)] = -u[(1, 2)];
        u[(2, 2)] = -u[(2, 2)];
    }
    u * vt
}

/// Kabsch rotation matrix (`kabsch`). Inputs assumed centered.
pub fn kabsch(p: &[Coord], q: &[Coord]) -> Matrix3<f64> {
    let c = covariance(p, q);
    let svd = c.svd(true, true);
    let u = svd.u.expect("SVD U");
    let vt = svd.v_t.expect("SVD Vt");
    optimal_rotation(u, vt)
}

/// Apply row-vector rotation `p @ u`.
pub fn apply_rotation(p: &Coord, u: &Matrix3<f64>) -> Coord {
    [
        p[0] * u[(0, 0)] + p[1] * u[(1, 0)] + p[2] * u[(2, 0)],
        p[0] * u[(0, 1)] + p[1] * u[(1, 1)] + p[2] * u[(2, 1)],
        p[0] * u[(0, 2)] + p[1] * u[(1, 2)] + p[2] * u[(2, 2)],
    ]
}

/// Rotate P onto Q (`kabsch_rotate`).
pub fn kabsch_rotate(p: &[Coord], q: &[Coord]) -> Vec<Coord> {
    let u = kabsch(p, q);
    p.iter().map(|v| apply_rotation(v, &u)).collect()
}

/// Rotate + translate P onto Q (`kabsch_fit`).
pub fn kabsch_fit(p: &[Coord], q: &[Coord]) -> Vec<Coord> {
    let qc = centroid(q);
    let pc = centroid(p);
    let pc_c: Vec<Coord> = p
        .iter()
        .map(|v| [v[0] - pc[0], v[1] - pc[1], v[2] - pc[2]])
        .collect();
    let qc_c: Vec<Coord> = q
        .iter()
        .map(|v| [v[0] - qc[0], v[1] - qc[1], v[2] - qc[2]])
        .collect();
    let u = kabsch(&pc_c, &qc_c);
    pc_c.iter()
        .map(|v| {
            let r = apply_rotation(v, &u);
            [r[0] + qc[0], r[1] + qc[1], r[2] + qc[2]]
        })
        .collect()
}

/// Kabsch RMSD with optional weights / centering (`kabsch_rmsd`).
pub fn kabsch_rmsd(p: &[Coord], q: &[Coord], w: Option<&[f64]>, translate: bool) -> f64 {
    if translate {
        let pc = centroid(p);
        let qc = centroid(q);
        let pp: Vec<Coord> = p
            .iter()
            .map(|v| [v[0] - pc[0], v[1] - pc[1], v[2] - pc[2]])
            .collect();
        let qq: Vec<Coord> = q
            .iter()
            .map(|v| [v[0] - qc[0], v[1] - qc[1], v[2] - qc[2]])
            .collect();
        if w.is_some() {
            kabsch_weighted_rmsd(&pp, &qq, w)
        } else {
            let rot = kabsch_rotate(&pp, &qq);
            rmsd(&rot, &qq)
        }
    } else if w.is_some() {
        kabsch_weighted_rmsd(p, q, w)
    } else {
        let rot = kabsch_rotate(p, q);
        rmsd(&rot, q)
    }
}

/// Weighted Kabsch (`kabsch_weighted`). Returns (U, V_row, rmsd) with
/// `P' = P @ U + V`, using the exact Python formula (including its
/// column-style translation `V = (CMP - U @ CMQ) / wsum`).
pub fn kabsch_weighted(p: &[Coord], q: &[Coord], w: Option<&[f64]>) -> (Matrix3<f64>, Coord, f64) {
    let n = p.len();
    assert_eq!(n, q.len());
    assert!(n > 0);
    let wvec: Vec<f64> = match w {
        Some(v) => {
            assert_eq!(v.len(), n);
            v.to_vec()
        }
        None => vec![1.0 / n as f64; n],
    };
    let wsum: f64 = wvec.iter().sum();
    let iw = 1.0 / wsum;

    let mut cmp = Vector3::zeros();
    let mut cmq = Vector3::zeros();
    let mut c = Matrix3::zeros();
    let mut psq = 0.0;
    let mut qsq = 0.0;
    for ((a, b), wi) in p.iter().zip(q.iter()).zip(wvec.iter()) {
        cmp[0] += a[0] * wi;
        cmp[1] += a[1] * wi;
        cmp[2] += a[2] * wi;
        cmq[0] += b[0] * wi;
        cmq[1] += b[1] * wi;
        cmq[2] += b[2] * wi;
        psq += (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) * wi;
        qsq += (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) * wi;
        for i in 0..3 {
            for j in 0..3 {
                c[(i, j)] += a[i] * wi * b[j];
            }
        }
    }
    psq -= cmp.dot(&cmp) * iw;
    qsq -= cmq.dot(&cmq) * iw;
    let outer = cmp * cmq.transpose();
    c = (c - outer * iw) * iw;

    let svd = c.svd(true, true);
    let u = svd.u.expect("SVD U");
    let vt = svd.v_t.expect("SVD Vt");
    let s = svd.singular_values;
    let d = u.determinant() * vt.determinant();
    let rot = optimal_rotation(u, vt);
    // Mirror Python's `S[-1] = -S[-1]` on reflection for the RMSD trace term.
    let s_sum = if d < 0.0 {
        s[0] + s[1] - s[2]
    } else {
        s[0] + s[1] + s[2]
    };

    let mut msd = (psq + qsq) * iw - 2.0 * s_sum;
    if msd < 0.0 {
        msd = 0.0;
    }
    let rms = msd.sqrt();
    let v = (cmp - rot * cmq) * iw;
    (rot, [v[0], v[1], v[2]], rms)
}

/// Weighted RMSD (`kabsch_weighted_rmsd`).
pub fn kabsch_weighted_rmsd(p: &[Coord], q: &[Coord], w: Option<&[f64]>) -> f64 {
    kabsch_weighted(p, q, w).2
}

/// Fit P to Q with optional weights (`kabsch_weighted_fit`).
/// Mirrors Python's swapped call `kabsch_weighted(Q, P, W)` + `P @ R.T + T`.
pub fn kabsch_weighted_fit(p: &[Coord], q: &[Coord], w: Option<&[f64]>) -> (Vec<Coord>, f64) {
    let (r, t, rms) = kabsch_weighted(q, p, w);
    let rt = r.transpose();
    let out: Vec<Coord> = p
        .iter()
        .map(|v| {
            let col = rt * Vector3::new(v[0], v[1], v[2]);
            [col[0] + t[0], col[1] + t[1], col[2] + t[2]]
        })
        .collect();
    (out, rms)
}

// --- Quaternion (Horn 1991) ---

fn make_w(r1: f64, r2: f64, r3: f64, r4: f64) -> Matrix4<f64> {
    Matrix4::new(
        r4, r3, -r2, r1, //
        -r3, r4, r1, r2, //
        r2, -r1, r4, r3, //
        -r1, -r2, -r3, r4,
    )
}

fn make_q(r1: f64, r2: f64, r3: f64, r4: f64) -> Matrix4<f64> {
    Matrix4::new(
        r4, -r3, r2, r1, //
        r3, r4, -r1, r2, //
        -r2, r1, r4, r3, //
        -r1, -r2, -r3, r4,
    )
}

/// Optimal rotation from unit quaternion (`quaternion_transform`).
pub fn quaternion_transform(r: Quat) -> Matrix3<f64> {
    let w = make_w(r[0], r[1], r[2], r[3]).transpose();
    let q = make_q(r[0], r[1], r[2], r[3]);
    let rot4 = w * q;
    Matrix3::new(
        rot4[(0, 0)],
        rot4[(0, 1)],
        rot4[(0, 2)], //
        rot4[(1, 0)],
        rot4[(1, 1)],
        rot4[(1, 2)], //
        rot4[(2, 0)],
        rot4[(2, 1)],
        rot4[(2, 2)],
    )
}

/// Quaternion rotation matrix (`quaternion_rotate`).
pub fn quaternion_rotate(x: &[Coord], y: &[Coord]) -> Matrix3<f64> {
    let mut a = Matrix4::zeros();
    for (px, py) in x.iter().zip(y.iter()) {
        let q = make_q(px[0], px[1], px[2], 0.0);
        let w = make_w(py[0], py[1], py[2], 0.0);
        a += q.transpose() * w;
    }
    let eig = SymmetricEigen::new(a);
    let mut best = 0;
    for i in 1..4 {
        if eig.eigenvalues[i] > eig.eigenvalues[best] {
            best = i;
        }
    }
    let col = eig.eigenvectors.column(best);
    quaternion_transform([col[0], col[1], col[2], col[3]])
}

/// Rotate P onto Q and report RMSD (`quaternion_rmsd`).
pub fn quaternion_rmsd(p: &[Coord], q: &[Coord]) -> f64 {
    let rot = quaternion_rotate(p, q);
    let pr: Vec<Coord> = p.iter().map(|v| apply_rotation(v, &rot)).collect();
    rmsd(&pr, q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn sample() -> (Vec<Coord>, Vec<Coord>) {
        let p = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let q = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];
        (p, q)
    }

    #[test]
    fn rmsd_zero() {
        let (p, _) = sample();
        assert_abs_diff_eq!(rmsd(&p, &p), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn kabsch_aligns_90deg() {
        let (p, q) = sample();
        let pc: Vec<Coord> = {
            let c = centroid(&p);
            p.iter()
                .map(|v| [v[0] - c[0], v[1] - c[1], v[2] - c[2]])
                .collect()
        };
        let qc: Vec<Coord> = {
            let c = centroid(&q);
            q.iter()
                .map(|v| [v[0] - c[0], v[1] - c[1], v[2] - c[2]])
                .collect()
        };
        assert_abs_diff_eq!(kabsch_rmsd(&pc, &qc, None, false), 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(quaternion_rmsd(&pc, &qc), 0.0, epsilon = 1e-9);
    }

    #[test]
    fn weighted_uniform_matches_plain() {
        let (p, q) = sample();
        let n = p.len() as f64;
        let w = vec![1.0 / n; p.len()];
        let r1 = kabsch_rmsd(&p, &q, None, true);
        let r2 = kabsch_weighted_rmsd(&p, &q, Some(&w));
        assert_abs_diff_eq!(r1, r2, epsilon = 1e-9);
    }
}
