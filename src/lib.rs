//! `rmsd` library — Rust port of `charnley/rmsd` (`calculate_rmsd.py`).
//!
//! Supported rotations: `kabsch`, `quaternion`, `none`.
//! Supported reorders: `none`, `hungarian`, `inertia-hungarian`,
//! `distance`, `brute`.
//!
//! QML similarity reorder is intentionally **not** ported (see README).

pub mod elements;
pub mod geom;
pub mod inertia;
pub mod io;
pub mod list;
pub mod matrix;
pub mod pipeline;
pub mod reflect;
pub mod reorder;

pub use elements::{atomic_number, mass, symbol};
pub use geom::{
    centroid, kabsch, kabsch_fit, kabsch_rmsd, kabsch_rotate, kabsch_weighted, kabsch_weighted_fit,
    kabsch_weighted_rmsd, quaternion_rmsd, quaternion_rotate, rmsd,
};
pub use pipeline::{compute_rmsd, Filter, Options, PipelineError, Rotation};
pub use reorder::ReorderMethod;
