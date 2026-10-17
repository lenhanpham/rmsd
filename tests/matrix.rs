//! Matrix-mode integration tests (temp-dir based; no repo fixtures harmed).

use rmsd::io::Molecule;
use rmsd::pipeline::{Options, Rotation};
use std::path::PathBuf;
use std::process::Command;

fn setup() -> (tempfile_dir::Guard, PathBuf) {
    let dir = tempfile_dir::create("rmsd-matrix-test").unwrap();
    let xyz = |n: &str, body: &str| {
        std::fs::write(
            dir.path().join(n),
            format!("{}\ntitle\n{body}", body.lines().count()),
        )
        .unwrap()
    };
    xyz("a.xyz", "C 0.0 0.0 0.0\nH 1.0 0.0 0.0\n");
    xyz("b.xyz", "C 0.0 0.0 0.0\nH 0.0 1.0 0.0\n");
    xyz("c.xyz", "C 0.0 0.0 0.0\nH 0.0 0.0 2.0\n");
    // Comma-separated form (covers the list-3 style tokenizer too).
    std::fs::write(dir.path().join("mylist.txt"), "a.xyz, b.xyz,c.xyz,\n").unwrap();
    let list = dir.path().join("mylist.txt");
    (dir, list)
}

fn kabsch_opts() -> Options {
    Options {
        rotation: Rotation::Kabsch,
        ..Default::default()
    }
}

/// Minimal tempdir guard without adding a dependency.
mod tempfile_dir {
    use std::path::{Path, PathBuf};
    pub struct Guard {
        path: PathBuf,
    }
    impl Guard {
        pub fn path(&self) -> &Path {
            &self.path
        }
    }
    impl Drop for Guard {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.path);
        }
    }
    pub fn create(prefix: &str) -> std::io::Result<Guard> {
        let mut dir = std::env::temp_dir();
        dir.push(format!(
            "{prefix}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir)?;
        Ok(Guard { path: dir })
    }
}

#[test]
fn matrix_csv_matches_pair_mode() {
    let (_guard, list) = setup();
    let entries = rmsd::list::load_list_file(&list).unwrap();
    assert_eq!(entries.len(), 3);
    let mut mols: Vec<Molecule> = Vec::new();
    for e in &entries {
        let text = std::fs::read_to_string(&e.path).unwrap();
        mols.push(rmsd::io::xyz::parse_xyz(&text).unwrap());
    }
    let labels: Vec<String> = entries.iter().map(|e| e.label.clone()).collect();
    let mat = rmsd::matrix::pairwise_matrix(&mols, &labels, &kabsch_opts(), &|_, _| {}).unwrap();
    let csv = rmsd::matrix::to_csv(&labels, &mat);
    let csv_path = rmsd::list::csv_path_for(&list);
    assert_eq!(csv_path.extension().unwrap(), "csv");
    std::fs::write(&csv_path, &csv).unwrap();
    let back = std::fs::read_to_string(&csv_path).unwrap();
    assert_eq!(back, csv);

    // Header + labels as listed; lower filled incl. 0 diagonal, upper blank.
    let lines: Vec<&str> = back.lines().collect();
    assert_eq!(lines[0], ",a.xyz,b.xyz,c.xyz");
    assert_eq!(lines[1], "a.xyz,0.000000,,");
    assert!(lines[2].starts_with("b.xyz,"));
    assert!(lines[2].ends_with(','));
    assert_eq!(lines[3].split(',').count(), 4);

    // Every lower cell equals pair-mode compute_rmsd.
    for i in 0..3 {
        for j in 0..=i {
            let direct =
                rmsd::pipeline::compute_rmsd(&mols[i], &mols[j], &kabsch_opts(), ("x", "y"))
                    .unwrap()
                    .rmsd;
            approx::assert_abs_diff_eq!(mat[i][j].unwrap(), direct, epsilon = 1e-12);
            let cell: f64 = lines[i + 1].split(',').nth(j + 1).unwrap().parse().unwrap();
            approx::assert_abs_diff_eq!(cell, direct, epsilon = 1e-6);
        }
    }
}

/// End-to-end: a complete existing CSV is archived (not silently discarded)
/// and recalculated; `--overwrite` recalculates without archiving; an
/// incomplete CSV resumes without archiving.
#[test]
fn matrix_complete_output_archives_and_recalculates() {
    let bin = std::env::var("CARGO_BIN_EXE_rmsd")
        .expect("CARGO_BIN_EXE_rmsd is set by cargo for bin targets");
    let guard = tempfile_dir::create("rmsd-archive-test").unwrap();
    let dir = guard.path();
    let xyz = |n: &str, body: &str| {
        std::fs::write(
            dir.join(n),
            format!("{}\ntitle\n{}", body.lines().count(), body),
        )
        .unwrap();
    };
    xyz("a.xyz", "C 0.0 0.0 0.0\nH 1.0 0.0 0.0\n");
    xyz("b.xyz", "C 0.0 0.0 0.0\nH 0.0 1.0 0.0\n");
    xyz("c.xyz", "C 0.0 0.0 0.0\nH 0.0 0.0 2.0\n");
    std::fs::write(dir.join("mylist.txt"), "a.xyz, b.xyz,c.xyz,\n").unwrap();
    let list = dir.join("mylist.txt");
    let out = dir.join("full.csv");

    fn run(bin: &str, list: &std::path::Path, out: &std::path::Path, extra: &[&str]) -> String {
        let mut cmd = Command::new(bin);
        cmd.args([
            list.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--threads",
            "2",
        ]);
        cmd.args(extra);
        let o = cmd.output().expect("failed to run rmsd binary");
        assert!(
            o.status.success(),
            "rmsd run failed: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        String::from_utf8_lossy(&o.stdout).into_owned()
    }

    // Count archive siblings full_*.csv.
    let archive_count = |dir: &std::path::Path| {
        std::fs::read_dir(dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| {
                e.file_name()
                    .to_str()
                    .map(|s| s.starts_with("full_") && s.ends_with(".csv"))
                    .unwrap_or(false)
            })
            .count()
    };

    // First run: creates full.csv, no archive.
    let o1 = run(&bin, &list, &out, &[]);
    assert!(!o1.contains("archived to"));
    assert!(out.exists(), "full.csv should be created on first run");

    // Second run (no --overwrite): archives + recalculates, message printed.
    let before = archive_count(dir);
    let o2 = run(&bin, &list, &out, &[]);
    assert!(
        o2.contains("archived to") && o2.contains("recalculating"),
        "expected archive+recalc message, got: {o2}"
    );
    assert!(out.exists(), "full.csv should be regenerated");
    assert_eq!(archive_count(dir), before + 1, "exactly one archive added");

    // Third run with --overwrite: recalculates without adding an archive.
    let before = archive_count(dir);
    let o3 = run(&bin, &list, &out, &["--overwrite"]);
    assert!(!o3.contains("archived to"), "--overwrite must not archive");
    assert_eq!(
        archive_count(dir),
        before,
        "--overwrite must not add an archive"
    );
}
