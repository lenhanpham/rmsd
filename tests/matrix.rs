//! Matrix-mode integration tests (temp-dir based; no repo fixtures harmed).

use rmsd::io::Molecule;
use rmsd::pipeline::{Options, Rotation};
use std::path::PathBuf;

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
