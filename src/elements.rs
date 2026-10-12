//! Element tables: atomic number <-> symbol, atomic masses.
//!
//! Port of `ELEMENT_WEIGHTS`, `ELEMENT_NAMES`, `NAMES_ELEMENT`,
//! `_ELEMENT_WEIGHTS_ARRAY`, `str_atom`, `int_atom` in `calculate_rmsd.py`.

use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq)]
pub enum ElementError {
    #[error("unknown element symbol: '{0}'")]
    UnknownSymbol(String),
    #[error("unknown atomic number: {0}")]
    UnknownNumber(u8),
}

/// Atomic masses indexed by atomic number. Index 0 unused.
/// Missing entries (113, 115, 117+) are 0.0, matching Python's
/// zero-filled `_ELEMENT_WEIGHTS_ARRAY` behaviour.
pub const MASS: [f64; 117] = {
    let mut m = [0.0; 117];
    m[1] = 1.00797;
    m[2] = 4.00260;
    m[3] = 6.941;
    m[4] = 9.01218;
    m[5] = 10.81;
    m[6] = 12.011;
    m[7] = 14.0067;
    m[8] = 15.9994;
    m[9] = 18.998403;
    m[10] = 20.179;
    m[11] = 22.98977;
    m[12] = 24.305;
    m[13] = 26.98154;
    m[14] = 28.0855;
    m[15] = 30.97376;
    m[16] = 32.06;
    m[17] = 35.453;
    m[18] = 39.948;
    m[19] = 39.0983;
    m[20] = 40.08;
    m[21] = 44.9559;
    m[22] = 47.90;
    m[23] = 50.9415;
    m[24] = 51.996;
    m[25] = 54.9380;
    m[26] = 55.847;
    m[27] = 58.9332;
    m[28] = 58.70;
    m[29] = 63.546;
    m[30] = 65.38;
    m[31] = 69.72;
    m[32] = 72.59;
    m[33] = 74.9216;
    m[34] = 78.96;
    m[35] = 79.904;
    m[36] = 83.80;
    m[37] = 85.4678;
    m[38] = 87.62;
    m[39] = 88.9059;
    m[40] = 91.22;
    m[41] = 92.9064;
    m[42] = 95.94;
    m[43] = 98.0;
    m[44] = 101.07;
    m[45] = 102.9055;
    m[46] = 106.4;
    m[47] = 107.868;
    m[48] = 112.41;
    m[49] = 114.82;
    m[50] = 118.69;
    m[51] = 121.75;
    m[52] = 127.60;
    m[53] = 126.9045;
    m[54] = 131.30;
    m[55] = 132.9054;
    m[56] = 137.33;
    m[57] = 138.9055;
    m[58] = 140.12;
    m[59] = 140.9077;
    m[60] = 144.24;
    m[61] = 145.0;
    m[62] = 150.4;
    m[63] = 151.96;
    m[64] = 157.25;
    m[65] = 158.9254;
    m[66] = 162.50;
    m[67] = 164.9304;
    m[68] = 167.26;
    m[69] = 168.9342;
    m[70] = 173.04;
    m[71] = 174.967;
    m[72] = 178.49;
    m[73] = 180.9479;
    m[74] = 183.85;
    m[75] = 186.207;
    m[76] = 190.2;
    m[77] = 192.22;
    m[78] = 195.09;
    m[79] = 196.9665;
    m[80] = 200.59;
    m[81] = 204.37;
    m[82] = 207.2;
    m[83] = 208.9804;
    m[84] = 209.0;
    m[85] = 210.0;
    m[86] = 222.0;
    m[87] = 223.0;
    m[88] = 226.0254;
    m[89] = 227.0278;
    m[90] = 232.0381;
    m[91] = 231.0359;
    m[92] = 238.029;
    m[93] = 237.0482;
    m[94] = 242.0;
    m[95] = 243.0;
    m[96] = 247.0;
    m[97] = 247.0;
    m[98] = 251.0;
    m[99] = 252.0;
    m[100] = 257.0;
    m[101] = 258.0;
    m[102] = 250.0;
    m[103] = 260.0;
    m[104] = 261.0;
    m[105] = 262.0;
    m[106] = 263.0;
    m[107] = 262.0;
    m[108] = 255.0;
    m[109] = 256.0;
    m[110] = 269.0;
    m[111] = 272.0;
    m[112] = 277.0;
    m
};

/// Element symbols indexed by atomic number. Index 0 is "".
pub const SYMBOL: [&str; 117] = {
    let mut s = [""; 117];
    s[1] = "H";
    s[2] = "He";
    s[3] = "Li";
    s[4] = "Be";
    s[5] = "B";
    s[6] = "C";
    s[7] = "N";
    s[8] = "O";
    s[9] = "F";
    s[10] = "Ne";
    s[11] = "Na";
    s[12] = "Mg";
    s[13] = "Al";
    s[14] = "Si";
    s[15] = "P";
    s[16] = "S";
    s[17] = "Cl";
    s[18] = "Ar";
    s[19] = "K";
    s[20] = "Ca";
    s[21] = "Sc";
    s[22] = "Ti";
    s[23] = "V";
    s[24] = "Cr";
    s[25] = "Mn";
    s[26] = "Fe";
    s[27] = "Co";
    s[28] = "Ni";
    s[29] = "Cu";
    s[30] = "Zn";
    s[31] = "Ga";
    s[32] = "Ge";
    s[33] = "As";
    s[34] = "Se";
    s[35] = "Br";
    s[36] = "Kr";
    s[37] = "Rb";
    s[38] = "Sr";
    s[39] = "Y";
    s[40] = "Zr";
    s[41] = "Nb";
    s[42] = "Mo";
    s[43] = "Tc";
    s[44] = "Ru";
    s[45] = "Rh";
    s[46] = "Pd";
    s[47] = "Ag";
    s[48] = "Cd";
    s[49] = "In";
    s[50] = "Sn";
    s[51] = "Sb";
    s[52] = "Te";
    s[53] = "I";
    s[54] = "Xe";
    s[55] = "Cs";
    s[56] = "Ba";
    s[57] = "La";
    s[58] = "Ce";
    s[59] = "Pr";
    s[60] = "Nd";
    s[61] = "Pm";
    s[62] = "Sm";
    s[63] = "Eu";
    s[64] = "Gd";
    s[65] = "Tb";
    s[66] = "Dy";
    s[67] = "Ho";
    s[68] = "Er";
    s[69] = "Tm";
    s[70] = "Yb";
    s[71] = "Lu";
    s[72] = "Hf";
    s[73] = "Ta";
    s[74] = "W";
    s[75] = "Re";
    s[76] = "Os";
    s[77] = "Ir";
    s[78] = "Pt";
    s[79] = "Au";
    s[80] = "Hg";
    s[81] = "Tl";
    s[82] = "Pb";
    s[83] = "Bi";
    s[84] = "Po";
    s[85] = "At";
    s[86] = "Rn";
    s[87] = "Fr";
    s[88] = "Ra";
    s[89] = "Ac";
    s[90] = "Th";
    s[91] = "Pa";
    s[92] = "U";
    s[93] = "Np";
    s[94] = "Pu";
    s[95] = "Am";
    s[96] = "Cm";
    s[97] = "Bk";
    s[98] = "Cf";
    s[99] = "Es";
    s[100] = "Fm";
    s[101] = "Md";
    s[102] = "No";
    s[103] = "Lr";
    s[104] = "Rf";
    s[105] = "Db";
    s[106] = "Sg";
    s[107] = "Bh";
    s[108] = "Hs";
    s[109] = "Mt";
    s[110] = "Ds";
    s[111] = "Rg";
    s[112] = "Cn";
    // 113, 115 missing upstream; 114/116 keep legacy names
    s[114] = "Uuq";
    s[116] = "Uuh";
    s
};

/// Convert atomic number to symbol (`str_atom`).
pub fn symbol(z: u8) -> Result<&'static str, ElementError> {
    let s = *SYMBOL.get(z as usize).unwrap_or(&"");
    if s.is_empty() {
        Err(ElementError::UnknownNumber(z))
    } else {
        Ok(s)
    }
}

/// Convert symbol to atomic number (`int_atom`).
/// Mirrors Python `atom.capitalize().strip()`.
pub fn atomic_number(sym: &str) -> Result<u8, ElementError> {
    let t = sym.trim();
    if t.is_empty() {
        return Err(ElementError::UnknownSymbol(sym.to_string()));
    }
    let mut chars = t.chars();
    let first = chars.next().unwrap().to_ascii_uppercase();
    let rest: String = chars.map(|c| c.to_ascii_lowercase()).collect();
    let norm = format!("{first}{rest}");
    for (i, s) in SYMBOL.iter().enumerate() {
        if *s == norm {
            return Ok(i as u8);
        }
    }
    Err(ElementError::UnknownSymbol(sym.to_string()))
}

/// Mass lookup, 0.0 for unknown (matches zero-filled array).
pub fn mass(z: u8) -> f64 {
    MASS.get(z as usize).copied().unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        assert_eq!(symbol(6).unwrap(), "C");
        assert_eq!(atomic_number("c").unwrap(), 6);
        assert_eq!(atomic_number(" Cl ").unwrap(), 17);
        assert!(atomic_number("1H").is_err());
    }

    #[test]
    fn masses() {
        assert!((mass(1) - 1.00797).abs() < 1e-9);
        assert!((mass(6) - 12.011).abs() < 1e-9);
    }
}
