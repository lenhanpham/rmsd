//! Pairwise RMSD matrix mode.
//!
//! Computes every pair once (`j <= i`; RMSD is symmetric, so the mirrored
//! half needs no recompute) and serializes a lower-triangular CSV with
//! label header row/column for Excel heatmap visualisation.
//!
//! Scale design (HPC headnode + queue systems):
//! - rows compute in parallel via Rayon; pairs within a row stay serial,
//!   so per-pair results are bit-identical to serial runs;
//! - output streams in row blocks — peak RAM is molecules + one block,
//!   never the full matrix;
//! - `--rows` range chunks let queue jobs split the row space; partial
//!   CSVs merge back with [`merge_partials`];
//! - an existing output CSV is archived
//!   (`<stem>_<YYYYMMDD_HHMMSS>[_<n>]<ext>` beside the list file) and the
//!   matrix is recomputed from scratch by default; pass `--resume` to
//!   validate it, truncate a torn tail, and append only the missing rows.

use crate::io::Molecule;
use crate::pipeline::{Options, PipelineError, compute_rmsd};
use rayon::prelude::*;
use std::fmt::Write as _;
use std::ops::Range;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum MatrixError {
    #[error("structures have different sizes ({0} vs {1}); a matrix needs uniform sizes")]
    Ragged(usize, usize),
    #[error(transparent)]
    Pipeline(#[from] PipelineError),
    #[error(
        "molecule set needs ~{needed} bytes but the memory budget is {cap} bytes \
         (raise --max-mem or split the run with --rows)"
    )]
    OverBudget { needed: u64, cap: u64 },
}

/// Rough in-memory footprint of one molecule (atoms + coords + overhead).
pub fn estimate_molecule_bytes(m: &Molecule) -> u64 {
    m.len() as u64 * 32 + 64
}

/// Total estimate for a molecule slice.
pub fn estimate_set_bytes(mols: &[Molecule]) -> u64 {
    mols.iter().map(estimate_molecule_bytes).sum()
}

/// Enforce the `--max-mem` budget before expensive work starts.
pub fn check_budget(mols: &[Molecule], cap: u64) -> Result<(), MatrixError> {
    let needed = estimate_set_bytes(mols);
    if needed > cap {
        return Err(MatrixError::OverBudget { needed, cap });
    }
    Ok(())
}

/// Parse `--max-mem` values like `512MB`, `4GB`, `1024` (bytes).
/// Suffixes are powers of 1024, case-insensitive, optional trailing `B`.
pub fn parse_mem_size(s: &str) -> Result<u64, String> {
    let t = s.trim();
    let cut = t.find(|c: char| c.is_alphabetic()).unwrap_or(t.len());
    let (num, suf) = t.split_at(cut);
    let n: u64 = num
        .parse()
        .map_err(|_| format!("invalid memory size: '{s}'"))?;
    let mult = match suf.to_ascii_uppercase().as_str() {
        "" | "B" => 1,
        "K" | "KB" => 1024,
        "M" | "MB" => 1024 * 1024,
        "G" | "GB" => 1024 * 1024 * 1024,
        _ => return Err(format!("invalid memory size: '{s}'")),
    };
    n.checked_mul(mult)
        .ok_or_else(|| format!("memory size overflows: '{s}'"))
}

/// Parse `--rows START:END` (0-based, END-exclusive) against list length.
pub fn parse_row_range(s: &str, n: usize) -> Result<Range<usize>, String> {
    let (a, b) = s
        .split_once(':')
        .ok_or_else(|| format!("rows must look like START:END, got '{s}'"))?;
    let start: usize = a
        .trim()
        .parse()
        .map_err(|_| format!("invalid row range: '{s}'"))?;
    let end: usize = b
        .trim()
        .parse()
        .map_err(|_| format!("invalid row range: '{s}'"))?;
    if start >= end || end > n {
        return Err(format!("row range '{s}' out of bounds for {n} structures"));
    }
    Ok(start..end)
}

/// Rows per output block: targets ~1/64 of the memory budget per block,
/// clamped to [1, 4096]. `row_bytes` is the estimated CSV bytes of one row.
pub fn block_rows(row_bytes: u64, max_mem: u64) -> usize {
    let target = max_mem / 64;
    let rows = target / row_bytes.max(1);
    rows.clamp(1, 4096) as usize
}

/// Lower-triangular RMSD rows for `range` (global indices): result
/// position `r` holds global row `range.start + r` as a full-width
/// `Some(j<=i)` / `None(j>i)` vector. Diagonal is exactly `0.0`.
/// Rows run in parallel; `progress(done_pairs, total_pairs)` is invoked
/// once per finished row (callers throttle printing themselves).
pub fn pairwise_rows(
    mols: &[Molecule],
    labels: &[String],
    opts: &Options,
    range: Range<usize>,
    progress: &(dyn Fn(usize, usize) + Sync),
) -> Result<Vec<Vec<Option<f64>>>, MatrixError> {
    // Chunked runs load only a molecule prefix but keep full labels for
    // mergeable full-width output.
    assert!(mols.len() <= labels.len());
    let n = labels.len();
    assert!(
        range.end <= mols.len(),
        "row range exceeds loaded molecules"
    );
    // Uniform sizes over the needed prefix, checked before any alignment.
    for m in mols.iter().take(range.end) {
        if m.len() != mols[0].len() {
            return Err(MatrixError::Ragged(mols[0].len(), m.len()));
        }
    }
    let total = n * (n + 1) / 2;
    let s = range.start;
    let base_done = s * (s + 1) / 2;
    range
        .into_par_iter()
        .map(|i| {
            let mut row = vec![None; n];
            for j in 0..=i {
                row[j] = Some(if i == j {
                    0.0
                } else {
                    compute_rmsd(
                        &mols[i],
                        &mols[j],
                        opts,
                        (labels[i].as_str(), labels[j].as_str()),
                    )?
                    .rmsd
                });
            }
            progress(base_done + row_pairs(i, s), total);
            Ok(row)
        })
        .collect()
}

/// Pairs completed when global row `i` finishes, given block start `s`:
/// rows s..=i hold (i(i+1)-s(s+1))/2 + ... — triangular numbers difference.
fn row_pairs(i: usize, s: usize) -> usize {
    (i + 1) * (i + 2) / 2 - s * (s + 1) / 2
}

/// Full-matrix convenience wrapper (serial-equivalent results).
pub fn pairwise_matrix(
    mols: &[Molecule],
    labels: &[String],
    opts: &Options,
    progress: &(dyn Fn(usize, usize) + Sync),
) -> Result<Vec<Vec<Option<f64>>>, MatrixError> {
    pairwise_rows(mols, labels, opts, 0..mols.len(), progress)
}

/// Serialize rows starting at global `offset` to CSV fragment (no header).
/// Each row is padded to full `labels` width: values for `j<=i` (rows may
/// be shorter when only a column prefix was loaded, e.g. range chunks),
/// blanks after.
pub fn rows_to_csv(labels: &[String], rows: &[Vec<Option<f64>>], offset: usize) -> String {
    let mut out = String::new();
    for (r, row) in rows.iter().enumerate() {
        out.push_str(&labels[offset + r]);
        for cell in row {
            out.push(',');
            if let Some(v) = cell {
                let _ = write!(out, "{v:.6}");
            }
        }
        for _ in row.len()..labels.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out
}

/// CSV header line (with trailing newline) for a label set.
pub fn csv_header(labels: &[String]) -> String {
    let mut out = String::new();
    out.push(',');
    out.push_str(&labels.join(","));
    out.push('\n');
    out
}

/// Serialize a full matrix (offset 0) including header.
pub fn to_csv(labels: &[String], matrix: &[Vec<Option<f64>>]) -> String {
    csv_header(labels) + &rows_to_csv(labels, matrix, 0)
}

/// Count complete leading rows in existing CSV `text` for `--resume`:
/// cuts a torn trailing line, checks the header, and requires data row
/// `k` to carry `labels[k]` (else the file belongs to another run).
/// Returns `(done_rows, keep_bytes)` where `keep_bytes` is the prefix to
/// retain (torn tail dropped) before appending.
pub fn resume_state(text: &str, labels: &[String]) -> Result<(usize, usize), String> {
    let cut = match text.rfind('\n') {
        Some(p) => p + 1,
        None => return Err("existing output has no complete lines".into()),
    };
    let body = &text[..cut];
    let mut lines = body.lines();
    let header = lines
        .next()
        .ok_or_else(|| "existing output is empty".to_string())?;
    if header != csv_header(labels).trim_end() {
        return Err(
            "existing output header does not match this structure list (use --overwrite)".into(),
        );
    }
    let mut done = 0usize;
    for line in lines {
        let (label, _) = line
            .split_once(',')
            .ok_or_else(|| format!("malformed row in existing output: '{line}'"))?;
        if done >= labels.len() || label != labels[done] {
            return Err(
                "existing output rows do not match this structure list (use --overwrite)".into(),
            );
        }
        done += 1;
    }
    Ok((done, cut))
}

/// Count complete leading rows in existing CSV `text` for `--resume`.
/// See [`resume_state`].
pub fn completed_rows(text: &str, labels: &[String]) -> Result<usize, String> {
    resume_state(text, labels).map(|(done, _)| done)
}

/// Merge partial CSVs (each: same header + a subset of data rows) into
/// `(labels, ordered data lines)`. Row order follows the header; duplicate
/// or missing rows are errors. Lines pass through byte-identical.
pub fn merge_partials(texts: &[String]) -> Result<(Vec<String>, Vec<String>), String> {
    if texts.is_empty() {
        return Err("nothing to merge".into());
    }
    let header = texts[0]
        .lines()
        .next()
        .ok_or_else(|| "empty partial output".to_string())?;
    if !header.starts_with(',') {
        return Err("partial output missing label header".into());
    }
    let labels: Vec<String> = header[1..].split(',').map(|s| s.to_string()).collect();
    let pos = |label: &str| {
        labels
            .iter()
            .position(|l| l == label)
            .ok_or_else(|| format!("unknown row label '{label}'"))
    };
    let mut rows: Vec<Option<String>> = vec![None; labels.len()];
    for text in texts {
        let mut lines = text.lines();
        let h = lines.next().unwrap_or("");
        if h != header {
            return Err("partial outputs have mismatching headers".into());
        }
        for line in lines {
            if line.is_empty() {
                continue;
            }
            let (label, _) = line
                .split_once(',')
                .ok_or_else(|| format!("malformed row: '{line}'"))?;
            let idx = pos(label)?;
            if rows[idx].is_some() {
                return Err(format!("duplicate row '{label}' in partial outputs"));
            }
            rows[idx] = Some(line.to_string());
        }
    }
    let mut out = Vec::with_capacity(labels.len());
    for (i, row) in rows.into_iter().enumerate() {
        out.push(row.ok_or_else(|| format!("missing row '{}' in partial outputs", labels[i]))?);
    }
    Ok((labels, out))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::Rotation;
    use crate::reorder::ReorderMethod;

    fn mol(atoms: Vec<u8>, coords: Vec<crate::geom::Coord>) -> Molecule {
        Molecule { atoms, coords }
    }

    fn kabsch_opts() -> Options {
        Options {
            rotation: Rotation::Kabsch,
            reorder: ReorderMethod::None,
            ..Default::default()
        }
    }

    fn noop(_: usize, _: usize) {}

    #[test]
    fn triangle_shape_and_diagonal() {
        let mols = vec![
            mol(vec![6], vec![[0.0, 0.0, 0.0]]),
            mol(vec![6], vec![[1.0, 0.0, 0.0]]),
            mol(vec![6], vec![[0.0, 2.0, 0.0]]),
        ];
        let labels = vec!["a".into(), "b".into(), "c".into()];
        let mat = pairwise_matrix(&mols, &labels, &kabsch_opts(), &noop).unwrap();
        for (i, row) in mat.iter().enumerate() {
            assert_eq!(row[i], Some(0.0));
            for cell in row.iter().skip(i + 1) {
                assert_eq!(*cell, None);
            }
        }
        for (i, j) in [(1, 0), (2, 0), (2, 1)] {
            let direct = compute_rmsd(&mols[i], &mols[j], &kabsch_opts(), ("a", "b")).unwrap();
            approx::assert_abs_diff_eq!(mat[i][j].unwrap(), direct.rmsd, epsilon = 1e-12);
        }
    }

    #[test]
    fn chunked_csv_merges_to_full() {
        // Simulates queue chunks: prefix-loaded rows serialize narrow but
        // pad to full width, and merge_partials reassembles the full CSV.
        let all = vec![
            mol(vec![6], vec![[0.0, 0.0, 0.0]]),
            mol(vec![6], vec![[1.0, 0.0, 0.0]]),
            mol(vec![6], vec![[0.0, 2.0, 0.0]]),
            mol(vec![6], vec![[3.0, 1.0, 0.0]]),
        ];
        let labels: Vec<String> = ["a", "b", "c", "d"].iter().map(|s| s.to_string()).collect();
        let full = pairwise_matrix(&all, &labels, &kabsch_opts(), &noop).unwrap();
        // Chunk loader keeps full labels but a molecule prefix.
        let top = pairwise_rows(&all[..2], &labels, &kabsch_opts(), 0..2, &noop).unwrap();
        let bottom = pairwise_rows(&all, &labels, &kabsch_opts(), 2..4, &noop).unwrap();
        let t1 = csv_header(&labels) + &rows_to_csv(&labels, &top, 0);
        let t2 = csv_header(&labels) + &rows_to_csv(&labels, &bottom, 2);
        let (merged_labels, rows) = merge_partials(&[t1, t2]).unwrap();
        assert_eq!(merged_labels, labels);
        let mut merged = csv_header(&merged_labels);
        for row in rows {
            merged.push_str(&row);
            merged.push('\n');
        }
        assert_eq!(merged, to_csv(&labels, &full));
    }

    #[test]
    fn range_chunk_concatenates_to_full() {
        let mols = vec![
            mol(vec![6], vec![[0.0, 0.0, 0.0]]),
            mol(vec![6], vec![[1.0, 0.0, 0.0]]),
            mol(vec![6], vec![[0.0, 2.0, 0.0]]),
            mol(vec![6], vec![[3.0, 1.0, 0.0]]),
        ];
        let labels: Vec<String> = ["a", "b", "c", "d"].iter().map(|s| s.to_string()).collect();
        let full = pairwise_matrix(&mols, &labels, &kabsch_opts(), &noop).unwrap();
        let top = pairwise_rows(&mols, &labels, &kabsch_opts(), 0..2, &noop).unwrap();
        let bottom = pairwise_rows(&mols, &labels, &kabsch_opts(), 2..4, &noop).unwrap();
        assert_eq!(&top, &full[0..2]);
        assert_eq!(&bottom, &full[2..4]);
    }

    #[test]
    fn csv_layout() {
        let labels = vec!["a".into(), "b".into()];
        let matrix = vec![vec![Some(0.0), None], vec![Some(1.5), Some(0.0)]];
        assert_eq!(
            to_csv(&labels, &matrix),
            ",a,b\na,0.000000,\nb,1.500000,0.000000\n"
        );
    }

    #[test]
    fn ragged_sizes_rejected() {
        let mols = vec![
            mol(vec![6], vec![[0.0, 0.0, 0.0]]),
            mol(vec![6, 6], vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        ];
        let labels = vec!["a".into(), "b".into()];
        assert!(pairwise_matrix(&mols, &labels, &kabsch_opts(), &noop).is_err());
    }

    #[test]
    fn mem_and_range_parsing() {
        assert_eq!(parse_mem_size("4GB").unwrap(), 4 * 1024 * 1024 * 1024);
        assert_eq!(parse_mem_size("512mb").unwrap(), 512 * 1024 * 1024);
        assert_eq!(parse_mem_size("1024").unwrap(), 1024);
        assert!(parse_mem_size("big").is_err());
        assert_eq!(parse_row_range("2:8", 16).unwrap(), 2..8);
        assert!(parse_row_range("8:2", 16).is_err());
        assert!(parse_row_range("0:17", 16).is_err());
        assert!(parse_row_range("nope", 16).is_err());
    }

    #[test]
    fn budget_gate() {
        let mols = vec![mol(vec![6], vec![[0.0, 0.0, 0.0]])];
        assert!(check_budget(&mols, 1).is_err());
        assert!(check_budget(&mols, u64::MAX).is_ok());
    }

    #[test]
    fn resume_counts_and_rejects() {
        let labels = vec!["a".into(), "b".into(), "c".into()];
        let text = ",a,b,c\na,0.000000,,\nb,1.000000,0.000000,\n";
        assert_eq!(completed_rows(text, &labels).unwrap(), 2);
        // Torn tail is cut, not counted.
        let torn = text.to_string() + "c,2.000000,1.500";
        assert_eq!(completed_rows(&torn, &labels).unwrap(), 2);
        // Wrong header / wrong order rejected.
        assert!(completed_rows(",x,y\na,0.000000,\n", &labels).is_err());
        assert!(completed_rows(",a,b,c\nb,1.000000,0.000000,\n", &labels).is_err());
    }

    #[test]
    fn merge_roundtrip() {
        let h = ",a,b,c\n".to_string();
        let p1 = h.clone() + "a,0.000000,,\n";
        let p2 = h.clone() + "b,1.000000,0.000000,\nc,2.000000,1.500000,0.000000\n";
        let (labels, rows) = merge_partials(&[p1.clone(), p2]).unwrap();
        assert_eq!(labels, ["a", "b", "c"]);
        assert_eq!(rows.len(), 3);
        assert!(rows[0].starts_with("a,"));
        // Duplicate / missing / header mismatch rejected.
        assert!(merge_partials(&[]).is_err());
        let dup = h.clone() + "a,0.000000,,\n";
        assert!(merge_partials(&[p1.clone(), dup]).is_err());
        assert!(merge_partials(std::slice::from_ref(&p1)).is_err()); // missing b, c
        let bad = ",a,b\n".to_string() + "a,0.000000,\n";
        assert!(merge_partials(&[p1.clone(), bad]).is_err());
    }
}
