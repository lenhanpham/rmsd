//! End-to-end pipeline — port of `main()` in `calculate_rmsd.py`.

use crate::geom::{kabsch, kabsch_rmsd, quaternion_rmsd, rmsd, Coord};
use crate::io::Molecule;
use crate::reflect::check_reflections;
use crate::reorder::{reorder_dispatch, ReorderMethod};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PipelineError {
    #[error("structures have different sizes ({0} vs {1})")]
    SizeMismatch(usize, usize),
    #[error("atoms are not in the same order; use reordering to align them")]
    AtomOrder,
    #[error("structure not aligned after reordering (internal error)")]
    NotAligned,
    #[error("cannot combine rotated output with atom filtering + reordering")]
    PrintFilterConflict,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rotation {
    Kabsch,
    Quaternion,
    None,
}

impl Rotation {
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "kabsch" => Some(Rotation::Kabsch),
            "quaternion" => Some(Rotation::Quaternion),
            "none" => Some(Rotation::None),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reflections {
    None,
    All,
    KeepStereo,
}

#[derive(Debug, Clone, Default)]
pub enum Filter {
    #[default]
    All,
    NoHydrogen,
    RemoveIdx(Vec<usize>),
    AddIdx(Vec<usize>),
}

#[derive(Debug, Clone)]
pub struct Options {
    pub rotation: Rotation,
    /// Reorder method; `ReorderMethod::None` disables reordering.
    pub reorder: ReorderMethod,
    /// Explicit `--reorder` switch (mirrors Python). If false, `reorder`
    /// is only applied inside reflection search.
    pub do_reorder: bool,
    pub reflections: Reflections,
    pub filter: Filter,
    /// Print rotated structure instead of just RMSD.
    pub print_rotated: bool,
    /// With filtering, print only the atoms used for RMSD.
    pub print_only_rmsd_atoms: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            rotation: Rotation::Kabsch,
            reorder: ReorderMethod::InertiaHungarian,
            do_reorder: false,
            reflections: Reflections::None,
            filter: Filter::All,
            print_rotated: false,
            print_only_rmsd_atoms: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Output {
    pub rmsd: f64,
    /// Transformed Q molecule for `--print` (None unless print_rotated).
    pub rotated: Option<Molecule>,
    pub title: String,
}

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

fn rmsd_fn_for(rotation: Rotation) -> fn(&[Coord], &[Coord]) -> f64 {
    match rotation {
        Rotation::Kabsch => |p, q| kabsch_rmsd(p, q, None, false),
        Rotation::Quaternion => quaternion_rmsd,
        Rotation::None => rmsd,
    }
}

fn reorder_fn_for(m: ReorderMethod) -> Option<fn(&[u8], &[u8], &[Coord], &[Coord]) -> Vec<usize>> {
    match m {
        ReorderMethod::None => None,
        ReorderMethod::Hungarian => Some(crate::reorder::reorder_hungarian),
        ReorderMethod::InertiaHungarian => Some(crate::reorder::reorder_inertia_hungarian),
        ReorderMethod::Distance => Some(crate::reorder::reorder_distance),
        ReorderMethod::Brute => Some(crate::reorder::reorder_brute),
    }
}

/// Core pipeline: filter → center → reorder/reflect → RMSD → optional transform.
pub fn compute_rmsd(
    p: &Molecule,
    q: &Molecule,
    opts: &Options,
    names: (&str, &str),
) -> Result<Output, PipelineError> {
    if p.len() != q.len() {
        return Err(PipelineError::SizeMismatch(p.len(), q.len()));
    }

    let reorder_active = opts.do_reorder && opts.reorder != ReorderMethod::None;
    let reflection_reorder = reorder_fn_for(opts.reorder);

    // Atom-order precheck (mirrors Python: error unless --reorder).
    // With reflections we still allow reorder via the reflection search.
    let same_order = p.atoms == q.atoms;
    if !same_order && !reorder_active && opts.reflections == Reflections::None {
        // If reflection search carries a reorder method it can still fix order.
        let can_fix = opts.reflections != Reflections::None && reflection_reorder.is_some();
        if !can_fix {
            return Err(PipelineError::AtomOrder);
        }
    }

    // --- views / filtering ---
    let (p_idx, q_idx, use_view): (Vec<usize>, Vec<usize>, bool) = match &opts.filter {
        Filter::All => ((0..p.len()).collect(), (0..q.len()).collect(), false),
        Filter::NoHydrogen => {
            let pi: Vec<usize> = p
                .atoms
                .iter()
                .enumerate()
                .filter(|(_, &z)| z != 1)
                .map(|(i, _)| i)
                .collect();
            let qi: Vec<usize> = q
                .atoms
                .iter()
                .enumerate()
                .filter(|(_, &z)| z != 1)
                .map(|(i, _)| i)
                .collect();
            (pi, qi, true)
        }
        Filter::RemoveIdx(ids) => {
            let drop: std::collections::BTreeSet<usize> = ids.iter().copied().collect();
            let keep: Vec<usize> = (0..p.len()).filter(|i| !drop.contains(i)).collect();
            (keep.clone(), keep, true)
        }
        Filter::AddIdx(ids) => (ids.clone(), ids.clone(), true),
    };

    let p_atoms: Vec<u8> = p_idx.iter().map(|&i| p.atoms[i]).collect();
    let mut q_atoms: Vec<u8> = q_idx.iter().map(|&i| q.atoms[i]).collect();
    let mut p_sub: Vec<Coord> = p_idx.iter().map(|&i| p.coords[i]).collect();
    let mut q_sub: Vec<Coord> = q_idx.iter().map(|&i| q.coords[i]).collect();

    // Recenter sub sets to centroid (mirrors Python main).
    let p_cent = centroid_of(&p_sub);
    let q_cent = centroid_of(&q_sub);
    for v in p_sub.iter_mut() {
        v[0] -= p_cent[0];
        v[1] -= p_cent[1];
        v[2] -= p_cent[2];
    }
    for v in q_sub.iter_mut() {
        v[0] -= q_cent[0];
        v[1] -= q_cent[1];
        v[2] -= q_cent[2];
    }

    let rmsd_fn = rmsd_fn_for(opts.rotation);

    let mut result_rmsd: Option<f64> = None;
    let mut q_swap: Option<[usize; 3]> = None;
    let mut q_refl: Option<[f64; 3]> = None;
    let mut q_review: Option<Vec<usize>> = None;

    match opts.reflections {
        Reflections::None => {}
        Reflections::All | Reflections::KeepStereo => {
            if opts.print_rotated && use_view && reorder_active {
                return Err(PipelineError::PrintFilterConflict);
            }
            let keep = opts.reflections == Reflections::KeepStereo;
            let (r, swap, refl, review) = check_reflections(
                &p_atoms,
                &q_atoms,
                &p_sub,
                &q_sub,
                reflection_reorder,
                rmsd_fn,
                keep,
            );
            result_rmsd = Some(r);
            q_swap = Some(swap);
            q_refl = Some(refl);
            q_review = Some(review);
        }
    }

    if result_rmsd.is_none() && reorder_active {
        if opts.print_rotated && use_view {
            return Err(PipelineError::PrintFilterConflict);
        }
        let review = reorder_dispatch(opts.reorder, &p_atoms, &q_atoms, &p_sub, &q_sub)
            .expect("reorder method is Some");
        q_review = Some(review);
    }

    if let Some(review) = q_review.clone() {
        q_atoms = review.iter().map(|&i| q_atoms[i]).collect();
        q_sub = review.iter().map(|&i| q_sub[i]).collect();
        if p_atoms != q_atoms {
            return Err(PipelineError::NotAligned);
        }
    }

    let rmsd_value = match result_rmsd {
        Some(r) => r,
        None => rmsd_fn(&p_sub, &q_sub),
    };

    // --- rotated output (mirrors Python print branch; always Kabsch for fit) ---
    let rotated = if opts.print_rotated {
        // Apply stored swap/reflection to sub set, then Kabsch fit.
        let mut qs = q_sub.clone();
        if let Some(swap) = q_swap {
            qs = qs
                .iter()
                .map(|v| [v[swap[0]], v[swap[1]], v[swap[2]]])
                .collect();
        }
        if let Some(refl) = q_refl {
            qs = qs
                .iter()
                .map(|v| [v[0] * refl[0], v[1] * refl[1], v[2] * refl[2]])
                .collect();
        }
        let u = kabsch(&qs, &p_sub);
        let apply = |v: &Coord| {
            [
                v[0] * u[(0, 0)] + v[1] * u[(1, 0)] + v[2] * u[(2, 0)],
                v[0] * u[(0, 1)] + v[1] * u[(1, 1)] + v[2] * u[(2, 1)],
                v[0] * u[(0, 2)] + v[1] * u[(1, 2)] + v[2] * u[(2, 2)],
            ]
        };

        if opts.print_only_rmsd_atoms || !use_view {
            let out: Vec<Coord> = qs
                .iter()
                .map(|v| {
                    let r = apply(v);
                    [r[0] + p_cent[0], r[1] + p_cent[1], r[2] + p_cent[2]]
                })
                .collect();
            Some(Molecule {
                atoms: q_atoms.clone(),
                coords: out,
            })
        } else {
            // Full set recentered on sub centroids.
            let mut qf: Vec<Coord> = q
                .coords
                .iter()
                .map(|v| [v[0] - q_cent[0], v[1] - q_cent[1], v[2] - q_cent[2]])
                .collect();
            if let Some(swap) = q_swap {
                qf = qf
                    .iter()
                    .map(|v| [v[swap[0]], v[swap[1]], v[swap[2]]])
                    .collect();
            }
            if let Some(refl) = q_refl {
                qf = qf
                    .iter()
                    .map(|v| [v[0] * refl[0], v[1] * refl[1], v[2] * refl[2]])
                    .collect();
            }
            let out: Vec<Coord> = qf
                .iter()
                .map(|v| {
                    let r = apply(v);
                    [r[0] + p_cent[0], r[1] + p_cent[1], r[2] + p_cent[2]]
                })
                .collect();
            Some(Molecule {
                atoms: q.atoms.clone(),
                coords: out,
            })
        }
    } else {
        None
    };

    // Silence unused warning for p_atoms reassignment parity.
    let _ = &p_atoms;

    Ok(Output {
        rmsd: rmsd_value,
        rotated,
        title: format!("Rotated '{}' to match '{}'", names.1, names.0),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::Molecule;

    fn mol(atoms: Vec<u8>, coords: Vec<Coord>) -> Molecule {
        Molecule { atoms, coords }
    }

    #[test]
    fn identical_is_zero() {
        let p = mol(vec![6, 1], vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let q = p.clone();
        let o = compute_rmsd(&p, &q, &Options::default(), ("a", "b")).unwrap();
        assert!(o.rmsd < 1e-12);
    }

    #[test]
    fn shuffled_needs_reorder() {
        // Same element types, shuffled positions: no error, but RMSD > 0
        // until Hungarian reorder is enabled.
        // Same element types, shuffled positions: asymmetric case so the
        // mislabeled correspondence cannot be fixed by rotation alone.
        let p = mol(
            vec![6, 6, 8],
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 10.0, 0.0]],
        );
        let q = mol(
            vec![6, 6, 8],
            vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 10.0, 0.0]],
        );
        let plain = compute_rmsd(&p, &q, &Options::default(), ("a", "b")).unwrap();
        assert!(plain.rmsd > 1e-6);
        let opts = Options {
            do_reorder: true,
            reorder: ReorderMethod::Hungarian,
            ..Default::default()
        };
        let o = compute_rmsd(&p, &q, &opts, ("a", "b")).unwrap();
        assert!(o.rmsd < 1e-9);
        // Different atom order without reorder must error.
        let pd = mol(vec![6, 8], vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]);
        let qd = mol(vec![8, 6], vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]]);
        assert!(compute_rmsd(&pd, &qd, &Options::default(), ("a", "b")).is_err());
    }
}
