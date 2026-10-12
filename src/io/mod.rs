//! File I/O: XYZ + PDB readers, XYZ writer.
//!
//! Port of `get_coordinates*`, `set_coordinates`, `_parse_pdb_*`
//! in `calculate_rmsd.py`.

pub mod pdb;
pub mod xyz;

use crate::elements::symbol;
use crate::geom::Coord;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum IoError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error: {0}")]
    Parse(String),
    #[error("unsupported format: '{0}' (expected xyz or pdb)")]
    Format(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Xyz,
    Pdb,
}

impl Format {
    pub fn parse(s: &str) -> Result<Self, IoError> {
        match s.to_ascii_lowercase().as_str() {
            "xyz" => Ok(Format::Xyz),
            "pdb" => Ok(Format::Pdb),
            other => Err(IoError::Format(other.to_string())),
        }
    }

    /// Detect from file suffixes, handling `.gz` (mirrors `parse_arguments`).
    pub fn detect(path: &Path) -> Result<(Self, bool), IoError> {
        let name = path.to_string_lossy().to_lowercase();
        let is_gzip = name.ends_with(".gz");
        let stem = if is_gzip {
            name.trim_end_matches(".gz")
        } else {
            &name
        };
        let ext = stem.rsplit('.').next().unwrap_or("");
        Ok((Format::parse(ext)?, is_gzip))
    }
}

/// Molecule: atomic numbers + coordinates.
#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: Vec<u8>,
    pub coords: Vec<Coord>,
}

impl Molecule {
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }
}

/// Read raw text, transparently handling gzip.
pub fn read_text(path: &Path, is_gzip: bool) -> Result<String, IoError> {
    use std::io::Read;
    if is_gzip {
        let f = std::fs::File::open(path)?;
        let mut dec = flate2::read::GzDecoder::new(f);
        let mut s = String::new();
        dec.read_to_string(&mut s)?;
        Ok(s)
    } else {
        Ok(std::fs::read_to_string(path)?)
    }
}

/// XYZ writer (`set_coordinates`), 8 decimals by default.
pub fn write_xyz(mol: &Molecule, title: &str, decimals: usize) -> String {
    let mut out = Vec::new();
    out.push(mol.len().to_string());
    out.push(title.to_string());
    for (z, p) in mol.atoms.iter().zip(mol.coords.iter()) {
        let sym = symbol(*z).unwrap_or("X");
        out.push(format!(
            "{:<2} {:>width$.prec$} {:>width$.prec$} {:>width$.prec$}",
            sym,
            p[0],
            p[1],
            p[2],
            width = 15,
            prec = decimals
        ));
    }
    out.join("\n")
}
