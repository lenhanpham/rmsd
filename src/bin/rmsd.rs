//! `rmsd` CLI — improved interface (see README for Python deltas).
//!
//! Two modes:
//! - pair mode: `rmsd FILE_A FILE_B` prints one RMSD (or rotated XYZ).
//! - matrix mode: `rmsd --from-list LIST` or `rmsd LIST` computes every
//!   pair once and writes a lower-triangular CSV beside the list file.

use clap::{ArgGroup, Parser};
use rmsd::io::{Format, Molecule, read_text, write_xyz};
use rmsd::pipeline::{Filter, Options, Reflections, Rotation};
use rmsd::reorder::ReorderMethod;
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
#[command(
    name = "rmsd",
    version,
    about = "Calculate RMSD between structure A and B (XYZ or PDB)"
)]
#[command(group(
    ArgGroup::new("filter")
        .args(["only_alpha_carbons", "ignore_hydrogen", "remove_idx", "add_idx"])
        .multiple(false)
))]
struct Cli {
    /// First structure (reference). With no FILE_B, treated as a structure list (matrix mode)
    file_a: Option<PathBuf>,
    /// Second structure (mobile). With --merge-csv, partial CSVs to merge
    file_b: Option<PathBuf>,
    /// Further partial CSVs for --merge-csv
    extra: Vec<PathBuf>,

    /// Structure list file for pairwise matrix mode
    #[arg(long, value_name = "LIST")]
    from_list: Option<PathBuf>,

    /// Rotation / alignment method
    #[arg(short, long, default_value = "kabsch", value_parser = ["kabsch", "quaternion", "none"])]
    rotation: String,

    /// Reorder atoms to find best correspondence
    #[arg(short = 'e', long)]
    reorder: bool,

    /// Reorder method (no `qml`; see README)
    #[arg(
        long,
        default_value = "inertia-hungarian",
        value_parser = ["none", "hungarian", "inertia-hungarian", "distance", "brute"]
    )]
    reorder_method: String,

    /// Scan reflections + axis swaps (affects stereochemistry)
    #[arg(long, conflicts_with = "use_reflections_keep_stereo")]
    use_reflections: bool,

    /// Scan reflections but keep stereochemistry
    #[arg(long, conflicts_with = "use_reflections")]
    use_reflections_keep_stereo: bool,

    /// Only alpha carbons (PDB only)
    #[arg(long)]
    only_alpha_carbons: bool,

    /// Ignore hydrogens
    #[arg(short = 'n', long = "ignore-hydrogen", alias = "no-hydrogen")]
    ignore_hydrogen: bool,

    /// Atom indices NOT to consider
    #[arg(long, num_args = 1..)]
    remove_idx: Option<Vec<usize>>,

    /// Atom indices to consider
    #[arg(long, num_args = 1..)]
    add_idx: Option<Vec<usize>>,

    /// Input format (auto-detected from extension by default)
    #[arg(long, value_parser = ["xyz", "pdb"])]
    format: Option<String>,

    /// Threads for matrix mode (default: half of available cores)
    #[arg(long)]
    threads: Option<usize>,

    /// Memory budget for matrix mode, e.g. 4GB, 512MB (default: 4GB)
    #[arg(long, default_value = "4GB")]
    max_mem: String,

    /// Matrix mode: only compute rows START..END (0-based, END-exclusive)
    /// for queue-system chunking
    #[arg(long, value_name = "START:END")]
    rows: Option<String>,

    /// Matrix mode: output CSV path (default: beside the list file)
    #[arg(long)]
    out: Option<PathBuf>,

    /// Matrix mode: discard existing output and recompute from scratch
    /// without archiving it
    #[arg(long, conflicts_with = "resume")]
    overwrite: bool,

    /// Matrix mode: resume an incomplete existing CSV instead of the
    /// default (archive the output, then recalculate from scratch).
    /// Validates the header and row labels against the structure list,
    /// truncates any torn trailing line, and appends only the remaining
    /// rows. Conflicts with `--overwrite`.
    #[arg(long, conflicts_with = "overwrite")]
    resume: bool,

    /// Merge partial matrix CSVs into OUT and exit
    #[arg(long, value_name = "OUT")]
    merge_csv: Option<PathBuf>,

    /// Print rotated structure B in XYZ format instead of just RMSD
    #[arg(short = 'p', long = "print", alias = "output")]
    print: bool,

    /// With filtering, print only atoms used for RMSD
    #[arg(long)]
    print_only_rmsd_atoms: bool,
}

fn load(path: &Path, fmt: Format, is_gzip: bool, only_alpha: bool) -> anyhow::Result<Molecule> {
    let text = read_text(path, is_gzip)?;
    match fmt {
        Format::Xyz => Ok(rmsd::io::xyz::parse_xyz(&text)?),
        Format::Pdb => Ok(rmsd::io::pdb::parse_pdb(&text, only_alpha)?),
    }
}

fn detect(path: &Path, forced: &Option<String>) -> anyhow::Result<(Format, bool)> {
    match forced {
        Some(f) => Ok((Format::parse(f)?, path.to_string_lossy().ends_with(".gz"))),
        None => Ok(Format::detect(path)?),
    }
}

fn build_options(cli: &Cli) -> anyhow::Result<Options> {
    let rotation = Rotation::parse(&cli.rotation).expect("clap-validated");
    let reorder = ReorderMethod::parse(&cli.reorder_method).expect("clap-validated");
    if cli.reorder && reorder == ReorderMethod::None {
        anyhow::bail!("--reorder requires a --reorder-method other than 'none'");
    }
    let reflections = if cli.use_reflections {
        Reflections::All
    } else if cli.use_reflections_keep_stereo {
        Reflections::KeepStereo
    } else {
        Reflections::None
    };
    let filter = if cli.ignore_hydrogen {
        Filter::NoHydrogen
    } else if let Some(ids) = cli.remove_idx.clone() {
        Filter::RemoveIdx(ids)
    } else if let Some(ids) = cli.add_idx.clone() {
        Filter::AddIdx(ids)
    } else {
        Filter::All
    };
    Ok(Options {
        rotation,
        reorder,
        do_reorder: cli.reorder,
        reflections,
        filter,
        print_rotated: cli.print,
        print_only_rmsd_atoms: cli.print_only_rmsd_atoms,
    })
}

fn run_pair(cli: &Cli, file_a: &Path, file_b: &Path) -> anyhow::Result<()> {
    let (fmt_a, gz_a) = detect(file_a, &cli.format)?;
    let (fmt_b, gz_b) = detect(file_b, &cli.format)?;
    if fmt_a != fmt_b {
        anyhow::bail!("both inputs must use the same format");
    }
    if cli.only_alpha_carbons && fmt_a != Format::Pdb {
        anyhow::bail!("--only-alpha-carbons only applies to PDB format");
    }
    let opts = build_options(cli)?;
    let p = load(file_a, fmt_a, gz_a, cli.only_alpha_carbons)?;
    let q = load(file_b, fmt_b, gz_b, cli.only_alpha_carbons)?;
    let out = rmsd::pipeline::compute_rmsd(
        &p,
        &q,
        &opts,
        (&file_a.to_string_lossy(), &file_b.to_string_lossy()),
    )?;
    if let Some(m) = out.rotated {
        let title = format!("{} with RMSD {:.8}", out.title, out.rmsd);
        println!("{}", write_xyz(&m, &title, 8));
    } else {
        println!("{}", out.rmsd);
    }
    Ok(())
}

/// Half the available cores (min 1): safe default on shared headnodes.
/// Overridable with `--threads` for queue allocations.
fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| (n.get() / 2).max(1))
        .unwrap_or(1)
}

fn setup_pool(threads: Option<usize>) -> anyhow::Result<()> {
    let n = threads.unwrap_or_else(default_threads);
    if n == 0 {
        anyhow::bail!("--threads must be at least 1");
    }
    // May already exist in tests; only the first build wins, which is fine.
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global();
    Ok(())
}

fn run_matrix(cli: &Cli, list_path: &Path) -> anyhow::Result<()> {
    if cli.print {
        anyhow::bail!("--print is only available in pair mode (two structures)");
    }
    if cli.merge_csv.is_some() {
        anyhow::bail!("--merge-csv cannot be combined with matrix computation");
    }
    setup_pool(cli.threads)?;
    let max_mem = rmsd::matrix::parse_mem_size(&cli.max_mem).map_err(|e| anyhow::anyhow!(e))?;
    let entries = rmsd::list::load_list_file(list_path)?;
    let opts = build_options(cli)?;

    let range_full = 0..entries.len();
    let range = match &cli.rows {
        Some(s) => {
            rmsd::matrix::parse_row_range(s, entries.len()).map_err(|e| anyhow::anyhow!(e))?
        }
        None => range_full,
    };
    // Rows need columns 0..=i: load the prefix the range depends on.
    let need = &entries[..range.end];
    let mut mols = Vec::with_capacity(need.len());
    let mut labels = Vec::with_capacity(need.len());
    for e in need {
        let (fmt, gz) = detect(&e.path, &cli.format)?;
        if cli.only_alpha_carbons && fmt != Format::Pdb {
            anyhow::bail!("--only-alpha-carbons only applies to PDB format");
        }
        mols.push(
            load(&e.path, fmt, gz, cli.only_alpha_carbons)
                .map_err(|err| anyhow::anyhow!("failed to load '{}': {err}", e.label))?,
        );
        labels.push(e.label.clone());
    }
    rmsd::matrix::check_budget(&mols, max_mem)?;
    let full_labels: Vec<String> = entries.iter().map(|e| e.label.clone()).collect();

    let out_path = cli
        .out
        .clone()
        .unwrap_or_else(|| rmsd::list::csv_path_for(list_path));
    // Create a missing output directory instead of failing with a bare
    // OS error (e.g. temp dirs on CI runners).
    if let Some(parent) = out_path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent).map_err(|e| {
            anyhow::anyhow!(
                "cannot create output directory '{}': {e}",
                parent.to_string_lossy()
            )
        })?;
    }
    let expected_header = rmsd::matrix::csv_header(&full_labels);

    let exists = out_path.exists();
    let mut start = range.start;
    // A file we (re)create always leads with the header; --resume appends to
    // the existing (torn-tail-truncated) file without one.
    let mut need_header = true;

    if exists {
        if cli.resume {
            // Opt-in resume: validate header/row labels, drop any torn
            // trailing line, then append only the rows still missing.
            let text = std::fs::read_to_string(&out_path)?;
            let (done, keep) = rmsd::matrix::resume_state(&text, &full_labels)
                .map_err(|e| anyhow::anyhow!("cannot resume '{out_path:?}': {e}"))?;
            if done < range.start || done > range.end {
                anyhow::bail!(
                    "existing output covers rows 0..{done} which does not extend this range (use --overwrite)"
                );
            }
            if keep as u64 != text.len() as u64 {
                std::fs::OpenOptions::new()
                    .write(true)
                    .open(&out_path)?
                    .set_len(keep as u64)?;
            }
            start = done;
            need_header = false;
        } else if cli.overwrite {
            // Full rewrite: discard the existing output without archiving.
            start = range.start;
        } else {
            // Default: always back up the existing output first, then
            // recalculate from scratch — never resume silently.
            let archived = archive_path_for(&out_path)?;
            std::fs::rename(&out_path, &archived).map_err(|e| {
                anyhow::anyhow!("cannot archive '{out_path:?}' to '{archived:?}': {e}")
            })?;
            println!(
                "rmsd: output '{}' already exists; archived to '{}' and recalculating from scratch",
                out_path.to_string_lossy(),
                archived.to_string_lossy(),
            );
            start = range.start;
        }
    }

    let total = entries.len() * (entries.len() + 1) / 2;
    let base_done = start * (start + 1) / 2;
    let shown = std::sync::atomic::AtomicUsize::new(base_done);
    let report = |done: usize, _total: usize| {
        // Rows finish out of order; display the high-water mark, throttled
        // to whole-percent steps.
        let prev = shown.fetch_max(done, std::sync::atomic::Ordering::Relaxed);
        let pct = |d: usize| d * 100 / total.max(1);
        if pct(done) > pct(prev) || done == total {
            eprintln!("[{done}/{total}] pairs computed");
        }
    };

    // A resumed range that is already complete has nothing left to append.
    if start >= range.end {
        println!(
            "output '{}' already complete ({}x{}, {total} pairs); nothing to resume",
            out_path.to_string_lossy(),
            entries.len(),
            entries.len(),
        );
        return Ok(());
    }
    let mut cur = start;
    // Row blocks bound output buffering: compute, append, free.
    // Fresh files (new path, --overwrite, or default re-calc) always start
    // with the header — including range chunks, which are standalone
    // mergeable partials. --resume appends to the truncated existing file.
    let row_bytes = (entries.len() as u64 + 1) * 14;
    let block = rmsd::matrix::block_rows(row_bytes, max_mem);
    while cur < range.end {
        let next = (cur + block).min(range.end);
        let rows = rmsd::matrix::pairwise_rows(&mols, &labels, &opts, cur..next, &report)?;
        if need_header {
            std::fs::write(&out_path, &expected_header)?;
            need_header = false;
        }
        let chunk = rmsd::matrix::rows_to_csv(&full_labels, &rows, cur);
        use std::io::Write as _;
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&out_path)?;
        f.write_all(chunk.as_bytes())?;
        cur = next;
    }
    println!(
        "wrote {} ({}/{} rows, {total} pairs)",
        out_path.to_string_lossy(),
        range.end - range.start,
        entries.len(),
    );
    Ok(())
}

fn run_merge(out: &Path, partials: &[PathBuf]) -> anyhow::Result<()> {
    if partials.len() < 2 {
        anyhow::bail!("--merge-csv needs at least two partial CSV files");
    }
    let mut texts = Vec::with_capacity(partials.len());
    for p in partials {
        texts.push(
            std::fs::read_to_string(p)
                .map_err(|e| anyhow::anyhow!("cannot read '{}': {e}", p.to_string_lossy()))?,
        );
    }
    let (labels, rows) =
        rmsd::matrix::merge_partials(&texts).map_err(|e| anyhow::anyhow!("cannot merge: {e}"))?;
    let mut csv = rmsd::matrix::csv_header(&labels);
    for row in rows {
        csv.push_str(&row);
        csv.push('\n');
    }
    std::fs::write(out, csv)?;
    println!(
        "merged {} partials into {} ({}x{})",
        partials.len(),
        out.to_string_lossy(),
        labels.len(),
        labels.len()
    );
    Ok(())
}

/// Build an archived path for `path` in the same directory:
/// `full.csv -> full_YYYYMMDD_HHMMSS.csv`, adding `_1`, `_2`, ... when a
/// same-timestamp archive already exists (same-second reruns).
fn archive_path_for(path: &Path) -> anyhow::Result<PathBuf> {
    let ts = timestamp_suffix();
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("output");
    let (stem, ext) = if let Some((s, e)) = name.rsplit_once('.') {
        (s, Some(e))
    } else {
        (name, None)
    };
    let build = |counter: Option<usize>| -> PathBuf {
        let archived = match counter {
            None => match ext {
                Some(e) => format!("{stem}_{ts}.{e}"),
                None => format!("{name}_{ts}"),
            },
            Some(c) => match ext {
                Some(e) => format!("{stem}_{ts}_{c}.{e}"),
                None => format!("{name}_{ts}_{c}"),
            },
        };
        path.with_file_name(archived)
    };
    let mut counter = None;
    for _ in 0..10_000 {
        let candidate = build(counter);
        if !candidate.exists() {
            return Ok(candidate);
        }
        counter = Some(counter.unwrap_or(0) + 1);
    }
    anyhow::bail!("cannot find a free archive name for '{path:?}'")
}

/// Current timestamp as `YYYYMMDD_HHMMSS` from UNIX seconds (pure std,
/// no `chrono` dependency).
fn timestamp_suffix() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = secs / 86400;
    let sod = secs % 86400;
    let hh = sod / 3600;
    let mm = (sod % 3600) / 60;
    let ss = sod % 60;
    // Civil calendar from days since 1970-01-01 (Howard Hinnant's algorithm).
    let z = days as i64 + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; // [0, 399]
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = doy - (153 * mp + 2) / 5 + 1; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 }; // [1, 12]
    let y = if m <= 2 { y + 1 } else { y };
    format!("{:04}{:02}{:02}_{:02}{:02}{:02}", y, m, d, hh, mm, ss)
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    if let Some(out) = &cli.merge_csv {
        // Merge mode: positionals are partial CSVs.
        if cli.from_list.is_some() {
            anyhow::bail!("--merge-csv cannot be combined with --from-list");
        }
        let mut partials: Vec<PathBuf> = Vec::new();
        partials.extend(cli.file_a.clone());
        partials.extend(cli.file_b.clone());
        partials.extend(cli.extra.clone());
        return run_merge(out, &partials);
    }
    if !cli.extra.is_empty() {
        anyhow::bail!("unexpected extra arguments (did you mean --merge-csv OUT PARTIALS...?)");
    }
    match (&cli.from_list, &cli.file_a, &cli.file_b) {
        (Some(_), Some(_), _) | (Some(_), _, Some(_)) => {
            anyhow::bail!("--from-list cannot be combined with FILE_A/FILE_B")
        }
        (Some(list), None, None) => run_matrix(&cli, list),
        (None, Some(a), Some(b)) => run_pair(&cli, a, b),
        (None, Some(list), None) => run_matrix(&cli, list),
        (None, None, _) => {
            anyhow::bail!("expected FILE_A FILE_B, a single LIST file, or --from-list LIST")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn timestamp_suffix_has_expected_format() {
        let ts = timestamp_suffix();
        let b = ts.as_bytes();
        // YYYYMMDD_HHMMSS = 15 chars.
        assert_eq!(b.len(), 15);
        assert_eq!(b[8], b'_');
        assert!(b[..8].iter().all(|x| x.is_ascii_digit()));
        assert!(b[9..].iter().all(|x| x.is_ascii_digit()));
    }

    #[test]
    fn archive_path_inserts_timestamp_before_extension() {
        let src = std::env::temp_dir().join("full.csv");
        let archived = archive_path_for(&src).unwrap();
        assert_eq!(archived.parent(), Some(std::env::temp_dir().as_path()));
        let b = archived.file_name().unwrap().to_str().unwrap().as_bytes();
        // name == "full_" (5) + ts (15) + ".csv" (4) = 24 bytes.
        assert_eq!(b.len(), 24);
        assert_eq!(&b[..5], b"full_");
        assert_eq!(&b[20..], b".csv");
        // timestamp's internal underscore at global index 5 + 8 = 13.
        assert_eq!(b[13], b'_');
        assert!(b[5..13].iter().all(|x| x.is_ascii_digit()));
        assert!(b[14..20].iter().all(|x| x.is_ascii_digit()));
    }

    #[test]
    fn archive_path_adds_counter_when_name_taken() {
        // Unique stem so this test never collides with other tests running
        // in the same second in the shared temp dir.
        let src = std::env::temp_dir().join("archive-counter-test.csv");
        // Drop any leftovers from a previous run.
        let _ = std::fs::remove_file(&src);
        let first = archive_path_for(&src).unwrap();
        let first_name = first.file_name().unwrap().to_str().unwrap().to_string();
        // Reserve the would-be archive; the next call must escalate to _1.
        std::fs::File::create(&first).unwrap();
        let second = archive_path_for(&src).unwrap();
        let second_name = second.file_name().unwrap().to_str().unwrap().to_string();
        assert_eq!(first.parent(), second.parent());
        assert!(second_name.starts_with(&first_name[..first_name.len() - 4])); // same "name_<ts>"
        assert!(second_name.ends_with("_1.csv"));
        // Reserve it too; the next call must escalate to _2.
        std::fs::File::create(&second).unwrap();
        let third = archive_path_for(&src).unwrap();
        let third_name = third.file_name().unwrap().to_str().unwrap().to_string();
        assert!(third_name.starts_with(&first_name[..first_name.len() - 4]));
        assert!(third_name.ends_with("_2.csv"));
        for f in [&first, &second, &third] {
            let _ = std::fs::remove_file(f);
        }
        let _ = std::fs::remove_file(&src);
    }
}
