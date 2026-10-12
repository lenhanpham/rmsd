//! XYZ reader — port of `get_coordinates_xyz[_lines]`.

use super::{IoError, Molecule};
use crate::elements::{atomic_number, symbol};
use crate::geom::Coord;
use regex::Regex;

fn capitalize(s: &str) -> String {
    let t = s.trim();
    let mut chars = t.chars();
    match chars.next() {
        None => String::new(),
        Some(f) => f.to_ascii_uppercase().to_string() + &chars.as_str().to_ascii_lowercase(),
    }
}

/// Parse XYZ from full text.
pub fn parse_xyz(text: &str) -> Result<Molecule, IoError> {
    let lines: Vec<&str> = text.lines().collect();
    parse_xyz_lines(&lines)
}

fn parse_xyz_lines(lines: &[&str]) -> Result<Molecule, IoError> {
    if lines.is_empty() {
        return Err(IoError::Parse(
            "Could not obtain the number of atoms in the .xyz file.".into(),
        ));
    }
    let n_atoms: usize = lines[0].trim().parse().map_err(|_| {
        IoError::Parse("Could not obtain the number of atoms in the .xyz file.".into())
    })?;

    let atom_re = Regex::new(r"[a-zA-Z]+").unwrap();
    let num_re = Regex::new(r"[-]?\d+\.\d*(?:[Ee][-\+]\d+)?").unwrap();

    let mut atoms: Vec<String> = Vec::new();
    let mut coords: Vec<Coord> = Vec::new();

    // Skip title line; stop after n_atoms records (mirrors Python).
    for (k, line) in lines.iter().skip(2).enumerate() {
        if k == n_atoms {
            break;
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let values: Vec<&str> = line.split_whitespace().collect();
        let (sym, numbers): (String, Vec<f64>) = if values.len() < 4 {
            let sym = atom_re
                .find(line)
                .map(|m| m.as_str().to_ascii_uppercase())
                .unwrap_or_default();
            let numbers: Vec<f64> = num_re
                .find_iter(line)
                .filter_map(|m| m.as_str().parse().ok())
                .collect();
            (sym, numbers)
        } else {
            let sym = values[0].to_string();
            let mut numbers = Vec::new();
            for v in &values[1..] {
                match v.parse::<f64>() {
                    Ok(x) => numbers.push(x),
                    Err(_) => break,
                }
            }
            (sym, numbers)
        };
        if numbers.len() >= 3 {
            coords.push([numbers[0], numbers[1], numbers[2]]);
            atoms.push(sym);
        } else {
            return Err(IoError::Parse(format!(
                "Reading the .xyz file failed in line {}. Please check the format.",
                k + 2
            )));
        }
    }

    // Integer atom types (e.g. "6") -> symbols; else capitalize.
    let all_int = atoms.iter().all(|a| a.parse::<u8>().is_ok());
    let mut z: Vec<u8> = Vec::with_capacity(atoms.len());
    if all_int && !atoms.is_empty() {
        for a in &atoms {
            let n: u8 = a.parse().unwrap();
            symbol(n).map_err(|_| IoError::Parse(format!("unknown atomic number: {n}")))?;
            z.push(n);
        }
    } else {
        for a in &atoms {
            let c = capitalize(a);
            z.push(
                atomic_number(&c)
                    .map_err(|_| IoError::Parse(format!("unknown element symbol: '{a}'")))?,
            );
        }
    }

    Ok(Molecule { atoms: z, coords })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_basic() {
        let m = parse_xyz("3\ntitle\nC 0 0 0\nH 1 0 0\nH 0 1 0\n").unwrap();
        assert_eq!(m.atoms, vec![6, 1, 1]);
        assert_eq!(m.coords.len(), 3);
    }

    #[test]
    fn parses_int_types() {
        let m = parse_xyz("2\nt\n6 0.0 0.0 0.0\n1 1.0 0.0 0.0\n").unwrap();
        assert_eq!(m.atoms, vec![6, 1]);
    }

    #[test]
    fn parses_nospace() {
        // e.g. "C0.0 0.0 0.0" fallback path
        let m = parse_xyz("1\nt\nC0.0 0.0 0.0\n").unwrap();
        assert_eq!(m.atoms, vec![6]);
    }
}
