//! Structure-list files for pairwise matrix mode.
//!
//! A list file contains structure filenames in any of these forms
//! (all three `tests/5dendralene/list-*.txt` variants parse identically):
//!
//! - one filename per line (`list-1.txt`),
//! - whitespace-separated on one line (`list-2.txt`),
//! - comma-separated on one line, optional trailing comma (`list-3.txt`).
//!
//! Tokenizer: split on commas and any whitespace, drop empties.
//! Labels are the tokens verbatim ("as listed"), so names cannot contain
//! commas or whitespace.

use std::path::{Path, PathBuf};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ListError {
    #[error("I/O error reading list '{0}': {1}")]
    Io(String, #[source] std::io::Error),
    #[error("list '{0}' contains no structure names")]
    Empty(String),
}

/// Split list text into structure-name tokens.
pub fn parse_list_text(text: &str) -> Vec<String> {
    text.split(|c: char| c == ',' || c.is_whitespace())
        .filter(|t| !t.is_empty())
        .map(|t| t.to_string())
        .collect()
}

/// One entry: label as listed + filesystem path (relative entries resolve
/// against the list file's parent directory).
#[derive(Debug, Clone)]
pub struct ListEntry {
    pub label: String,
    pub path: PathBuf,
}

/// Read + tokenize a list file, resolving entries to paths.
pub fn load_list_file(list_path: &Path) -> Result<Vec<ListEntry>, ListError> {
    let text = std::fs::read_to_string(list_path)
        .map_err(|e| ListError::Io(list_path.to_string_lossy().into_owned(), e))?;
    let tokens = parse_list_text(&text);
    if tokens.is_empty() {
        return Err(ListError::Empty(list_path.to_string_lossy().into_owned()));
    }
    let parent = list_path.parent().filter(|p| !p.as_os_str().is_empty());
    Ok(tokens
        .into_iter()
        .map(|label| {
            let raw = PathBuf::from(&label);
            let path = match parent {
                Some(dir) if raw.is_relative() => dir.join(raw),
                _ => raw,
            };
            ListEntry { label, path }
        })
        .collect())
}

/// CSV path for a list file: same dir + basename, `.csv` extension
/// (`list-2.txt` -> `list-2.csv`).
pub fn csv_path_for(list_path: &Path) -> PathBuf {
    list_path.with_extension("csv")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_per_line() {
        let v = parse_list_text("a.xyz\nb.xyz\nc.xyz\n");
        assert_eq!(v, ["a.xyz", "b.xyz", "c.xyz"]);
    }

    #[test]
    fn space_separated_single_line() {
        let v = parse_list_text("a.xyz  b.xyz\tc.xyz\n");
        assert_eq!(v, ["a.xyz", "b.xyz", "c.xyz"]);
    }

    #[test]
    fn comma_separated_with_trailing_comma() {
        let v = parse_list_text("a.xyz, b.xyz,c.xyz,\n");
        assert_eq!(v, ["a.xyz", "b.xyz", "c.xyz"]);
    }

    #[test]
    fn blanks_and_mixed_separators() {
        let v = parse_list_text("\n  a.xyz,\n\nb.xyz  c.xyz,,\n");
        assert_eq!(v, ["a.xyz", "b.xyz", "c.xyz"]);
    }

    #[test]
    fn parses_real_list_files() {
        // All three repo list formats must yield the same 16 structures
        // (list-2/3 order odd-then-even, list-1 sequential — same set).
        let base = Path::new("tests/5dendralene");
        let l1 = load_list_file(&base.join("list-1.txt")).unwrap();
        let l2 = load_list_file(&base.join("list-2.txt")).unwrap();
        let l3 = load_list_file(&base.join("list-3.txt")).unwrap();
        for l in [&l1, &l2, &l3] {
            assert_eq!(l.len(), 16);
        }
        let mut s1: Vec<_> = l1.iter().map(|e| &e.label).collect();
        let mut s2: Vec<_> = l2.iter().map(|e| &e.label).collect();
        let mut s3: Vec<_> = l3.iter().map(|e| &e.label).collect();
        s1.sort();
        s2.sort();
        s3.sort();
        assert_eq!(s1, s2);
        assert_eq!(s1, s3);
        assert_eq!(l1[0].label, "5dendralene_01.xyz");
        assert_eq!(l1[15].label, "5dendralene_16.xyz");
        // Relative entries resolve against the list file's directory.
        assert_eq!(l1[0].path, base.join("5dendralene_01.xyz"));
    }

    #[test]
    fn csv_path_replaces_extension() {
        assert_eq!(
            csv_path_for(Path::new("tests/5dendralene/list-2.txt")),
            PathBuf::from("tests/5dendralene/list-2.csv")
        );
    }
}
