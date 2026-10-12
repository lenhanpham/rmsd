//! Center of mass, inertia tensor, principal axis.
//!
//! Port of `get_cm`, `get_inertia_tensor`, `get_principal_axis`,
//! `rotation_matrix_vectors` in `calculate_rmsd.py`.

use crate::elements::mass;
use crate::geom::Coord;
use nalgebra::{Matrix3, SymmetricEigen, Vector3};

/// Mass-weighted center (`get_cm`).
pub fn center_of_mass(atoms: &[u8], coords: &[Coord]) -> Coord {
    assert_eq!(atoms.len(), coords.len());
    let mut wsum = 0.0;
    let mut acc = [0.0; 3];
    for (z, p) in atoms.iter().zip(coords.iter()) {
        let w = mass(*z);
        wsum += w;
        acc[0] += p[0] * w;
        acc[1] += p[1] * w;
        acc[2] += p[2] * w;
    }
    // Python `np.average` divides by wsum; fall back to centroid if all masses 0.
    if wsum == 0.0 {
        let n = coords.len() as f64;
        return [
            acc[0] / n.max(1.0),
            acc[1] / n.max(1.0),
            acc[2] / n.max(1.0),
        ];
    }
    [acc[0] / wsum, acc[1] / wsum, acc[2] / wsum]
}

/// Inertia tensor (`get_inertia_tensor`).
pub fn inertia_tensor(atoms: &[u8], coords: &[Coord]) -> Matrix3<f64> {
    let cm = center_of_mass(atoms, coords);
    let mut helper = Matrix3::<f64>::zeros();
    for (z, p) in atoms.iter().zip(coords.iter()) {
        let w = mass(*z);
        let x = p[0] - cm[0];
        let y = p[1] - cm[1];
        let zc = p[2] - cm[2];
        helper[(0, 0)] += x * w * x;
        helper[(0, 1)] += x * w * y;
        helper[(0, 2)] += x * w * zc;
        helper[(1, 0)] += y * w * x;
        helper[(1, 1)] += y * w * y;
        helper[(1, 2)] += y * w * zc;
        helper[(2, 0)] += zc * w * x;
        helper[(2, 1)] += zc * w * y;
        helper[(2, 2)] += zc * w * zc;
    }
    Matrix3::identity() * helper.trace() - helper
}

/// Principal axis: eigenvector of largest eigenvalue (`get_principal_axis`).
pub fn principal_axis(atoms: &[u8], coords: &[Coord]) -> Coord {
    let inertia = inertia_tensor(atoms, coords);
    let eig = SymmetricEigen::new(inertia);
    let mut best = 0;
    for i in 1..3 {
        if eig.eigenvalues[i] > eig.eigenvalues[best] {
            best = i;
        }
    }
    let col = eig.eigenvectors.column(best);
    [col[0], col[1], col[2]]
}

/// Rotation matrix taking v1 onto v2 (Rodrigues).
/// Port of `rotation_matrix_vectors` (kept for API completeness;
/// unused by the CLI pipeline, as in Python).
pub fn rotation_matrix_vectors(v1: Coord, v2: Coord) -> Matrix3<f64> {
    if v1 == v2 {
        return Matrix3::identity();
    }
    if v1 == [-v2[0], -v2[1], -v2[2]] {
        return Matrix3::new(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0);
    }
    let a = Vector3::new(v1[0], v1[1], v1[2]);
    let b = Vector3::new(v2[0], v2[1], v2[2]);
    let v = a.cross(&b);
    let s = v.norm();
    let c: f64 = a.dot(&b);
    let vx = Matrix3::new(0.0, -v[2], v[1], v[2], 0.0, -v[0], -v[1], v[0], 0.0);
    Matrix3::identity() + vx + vx * vx * ((1.0 - c) / (s * s))
}

/// Eigenvectors of a symmetric 3x3, columns sorted by ascending eigenvalue.
/// Mirrors the `eigh` + `argsort` idiom in `reorder_inertia_hungarian`.
pub fn sorted_eigenvectors(m: &Matrix3<f64>) -> Matrix3<f64> {
    let eig = SymmetricEigen::new(*m);
    let mut idx = [0, 1, 2];
    idx.sort_by(|&a, &b| {
        eig.eigenvalues[a]
            .partial_cmp(&eig.eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let e = &eig.eigenvectors;
    Matrix3::from_columns(&[
        e.column(idx[0]).clone_owned(),
        e.column(idx[1]).clone_owned(),
        e.column(idx[2]).clone_owned(),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cm_water() {
        // O(8) at origin, 2x H(1): CM on the bisector.
        let atoms = [8, 1, 1];
        let coords = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let cm = center_of_mass(&atoms, &coords);
        assert!((cm[0]).abs() < 1e-9);
        assert!(cm[1] > 0.0 && cm[1] < 0.586);
    }
}
