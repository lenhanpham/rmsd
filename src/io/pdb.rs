//! PDB reader — port of `get_coordinates_pdb` + `_parse_pdb_*`.
//!
//! Tolerant like Python: column slices first, token scan as fallback,
//! stops at TER/END, reads first chain block.

use super::{IoError, Molecule};
use crate::elements::atomic_number;
use crate::geom::Coord;

fn col(line: &str, a: usize, b: usize) -> &str {
    if line.len() < b {
        return "";
    }
    &line[a..b.min(line.len())]
}

fn clean2(s: &str) -> String {
    let no_digits: String = s.chars().filter(|c| !c.is_ascii_digit()).collect();
    let t = no_digits.trim();
    let mut chars = t.chars();
    match chars.next() {
        None => String::new(),
        Some(f) => f.to_ascii_uppercase().to_string() + &chars.as_str().to_ascii_lowercase(),
    }
}

fn known(sym: &str) -> bool {
    atomic_number(sym).is_ok()
}

/// Port of `_parse_pdb_atom_line`. Returns None if unparseable.
pub fn parse_atom_line(line: &str) -> Option<String> {
    let atom_col = col(line, 12, 16);
    let mut atom = clean2(atom_col.get(0..2).unwrap_or(""));
    // H* is almost surely hydrogen, not He/Hf/Hg (PDB v2.1 note).
    if atom.len() == 2 && atom.starts_with('H') {
        atom = "H".to_string();
    }
    if !atom.is_empty() && known(&atom) {
        return Some(atom);
    }
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 3 {
        return None;
    }
    let name = tokens[2];
    let c0: String = name.chars().take(1).collect();
    if known(&c0) {
        return Some(c0);
    }
    // e.g. 1HD1 -> second char H
    let c1: String = name.chars().skip(1).take(1).collect();
    if c1.eq_ignore_ascii_case("H") {
        return Some(c1);
    }
    None
}

/// Port of `_parse_pdb_alphacarbon_line`.
pub fn is_alpha_carbon(line: &str) -> bool {
    let atom_col = col(line, 12, 16);
    let atom = clean2(atom_col.get(0..2).unwrap_or(""));
    let loc = atom_col.chars().nth(2).unwrap_or(' ');
    atom == "C" && loc == 'A'
}

/// Port of `_parse_pdb_coord_line`.
pub fn parse_coord_line(line: &str) -> Option<Coord> {
    // Fixed columns first.
    let x = col(line, 30, 38).trim();
    let y = col(line, 38, 46).trim();
    let z = col(line, 46, 54).trim();
    if let (Ok(a), Ok(b), Ok(c)) = (x.parse::<f64>(), y.parse::<f64>(), z.parse::<f64>()) {
        return Some([a, b, c]);
    }
    // Fallback: first 3 consecutive tokens containing '.'.
    let tokens: Vec<&str> = line.split_whitespace().collect();
    for i in 0..tokens.len().saturating_sub(2) {
        if tokens[i].contains('.')
            && tokens[i + 1].contains('.')
            && tokens[i + 2].contains('.')
            && let (Ok(a), Ok(b), Ok(c)) = (
                tokens[i].parse::<f64>(),
                tokens[i + 1].parse::<f64>(),
                tokens[i + 2].parse::<f64>(),
            )
        {
            return Some([a, b, c]);
        }
    }
    None
}

/// Parse PDB text.
pub fn parse_pdb(text: &str, only_alpha_carbon: bool) -> Result<Molecule, IoError> {
    let mut symbols: Vec<String> = Vec::new();
    let mut coords: Vec<Coord> = Vec::new();
    let mut alpha: Vec<bool> = Vec::new();

    for line in text.lines() {
        if line.starts_with("TER") || line.starts_with("END") {
            break;
        }
        if !(line.starts_with("ATOM") || line.starts_with("HETATM")) {
            continue;
        }
        let atom = parse_atom_line(line)
            .ok_or_else(|| IoError::Parse(format!("error: Parsing for atom line: {line}")))?;
        let coord = parse_coord_line(line).ok_or_else(|| {
            IoError::Parse(format!("error: Parsing coordinates for line: {line}"))
        })?;
        symbols.push(atom);
        coords.push(coord);
        alpha.push(is_alpha_carbon(line));
    }

    let mut atoms: Vec<u8> = Vec::with_capacity(symbols.len());
    for s in &symbols {
        atoms
            .push(atomic_number(s).map_err(|_| IoError::Parse(format!("unknown element: '{s}'")))?);
    }

    if only_alpha_carbon {
        if !alpha.iter().any(|&x| x) {
            return Err(IoError::Parse(
                "Trying to filter for alpha carbons, but couldn't find any".into(),
            ));
        }
        let mut fa = Vec::new();
        let mut fc = Vec::new();
        for (i, &is_a) in alpha.iter().enumerate() {
            if is_a {
                fa.push(atoms[i]);
                fc.push(coords[i]);
            }
        }
        return Ok(Molecule {
            atoms: fa,
            coords: fc,
        });
    }

    Ok(Molecule { atoms, coords })
}

#[cfg(test)]
mod tests {
    use super::*;

    const PDB: &str = "\
ATOM      1  N   ALA A   1      11.104  13.207   2.100\n\
ATOM      2  CA  ALA A   1      12.560  13.207   2.100\n\
ATOM      5  C   ALA A   1      13.000  14.000   3.000\n\
TER\n";

    #[test]
    fn parses() {
        let m = parse_pdb(PDB, false).unwrap();
        assert_eq!(m.atoms, vec![7, 6, 6]);
        assert_eq!(m.coords.len(), 3);
    }

    #[test]
    fn alpha_filter() {
        let m = parse_pdb(PDB, true).unwrap();
        assert_eq!(m.len(), 1);
        assert_eq!(m.atoms, vec![6]);
    }
}
