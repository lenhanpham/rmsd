//! `rmsd` CLI — improved interface (see README for Python deltas).
//!
//! Two modes:
//! - pair mode: `rmsd FILE_A FILE_B` prints one RMSD (or rotated XYZ).
//! - matrix mode: `rmsd --from-list LIST` or `rmsd LIST` computes every
//!   pair once and writes a lower-triangular CSV beside the list file.

use clap::{ArgGroup, Parser};
use rmsd::io::{read_text, write_xyz, Format, Molecule};
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
    #[arg(long)]
    overwrite: bool,

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
    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| {
                anyhow::anyhow!(
                    "cannot create output directory '{}': {e}",
                    parent.to_string_lossy()
                )
            })?;
        }
    }
    let expected_header = rmsd::matrix::csv_header(&full_labels);

    // Resume: count complete leading rows already on disk, dropping any
    // torn trailing line by truncating the file to the resume point.
    let exists = out_path.exists();
    let mut start = range.start;
    if exists && !cli.overwrite {
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

    if start >= range.end {
        println!(
            "output '{}' already complete ({}x{}, {total} pairs)",
            out_path.to_string_lossy(),
            entries.len(),
            entries.len(),
        );
        return Ok(());
    }
    // Row blocks bound output buffering: compute, append, free.
    // Fresh files (new path, or --overwrite) always start with the header —
    // including range chunks, which are standalone mergeable partials.
    let row_bytes = (entries.len() as u64 + 1) * 14;
    let block = rmsd::matrix::block_rows(row_bytes, max_mem);
    let mut need_header = !exists || cli.overwrite;
    let mut cur = start;
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
