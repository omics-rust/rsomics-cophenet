//! Differential tests against scipy.cluster.hierarchy.cophenet.
//!
//! The committed golden (generated once by scipy) always runs. The live
//! oracle diff loud-skips when no scipy python is on PATH / RSOMICS_SCIPY.

use std::path::{Path, PathBuf};
use std::process::Command;

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn bin() -> PathBuf {
    Path::new(env!("CARGO_BIN_EXE_rsomics-cophenet")).to_path_buf()
}

fn run(args: &[&str]) -> (String, String) {
    let out = Command::new(bin())
        .args(args)
        .output()
        .expect("run rsomics-cophenet");
    assert!(out.status.success(), "binary failed: {args:?}");
    (
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8(out.stderr).unwrap(),
    )
}

fn read(p: &Path) -> String {
    std::fs::read_to_string(p).unwrap()
}

#[test]
fn golden_condensed_byte_exact() {
    let g = golden_dir();
    let (stdout, _) = run(&[g.join("linkage.tsv").to_str().unwrap()]);
    assert_eq!(stdout, read(&g.join("cophenetic.condensed.tsv")));
}

#[test]
fn golden_square_byte_exact() {
    let g = golden_dir();
    let (stdout, _) = run(&[g.join("linkage.tsv").to_str().unwrap(), "--square"]);
    assert_eq!(stdout, read(&g.join("cophenetic.square.tsv")));
}

#[test]
fn golden_correlation_exact() {
    let g = golden_dir();
    let (_, stderr) = run(&[
        g.join("linkage.tsv").to_str().unwrap(),
        "--dist",
        g.join("orig.condensed.tsv").to_str().unwrap(),
    ]);
    let ours: f64 = stderr.rsplit(": ").next().unwrap().trim().parse().unwrap();
    let want: f64 = read(&g.join("correlation.txt")).trim().parse().unwrap();
    assert!((ours - want).abs() < 1e-9, "ours={ours} want={want}");
}

#[test]
fn golden_square_dist_matches_condensed_dist() {
    let g = golden_dir();
    let stderr_c = run(&[
        g.join("linkage.tsv").to_str().unwrap(),
        "--dist",
        g.join("orig.condensed.tsv").to_str().unwrap(),
    ])
    .1;
    let stderr_s = run(&[
        g.join("linkage.tsv").to_str().unwrap(),
        "--dist",
        g.join("orig.square.tsv").to_str().unwrap(),
        "--dist-square",
    ])
    .1;
    assert_eq!(stderr_c, stderr_s);
}

fn scipy_python() -> Option<String> {
    if let Ok(p) = std::env::var("RSOMICS_SCIPY") {
        return Some(p);
    }
    for cand in ["python3", "python"] {
        let ok = Command::new(cand)
            .args(["-c", "import scipy.cluster.hierarchy"])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if ok {
            return Some(cand.to_string());
        }
    }
    None
}

#[test]
fn live_scipy_differential() {
    let Some(py) = scipy_python() else {
        eprintln!(
            "SKIP live_scipy_differential: no scipy python (set RSOMICS_SCIPY=/path/to/python)"
        );
        return;
    };

    let scratch = std::env::temp_dir().join("rsomics-cophenet-compat");
    std::fs::create_dir_all(&scratch).unwrap();
    let z = scratch.join("Z.tsv");
    let y = scratch.join("Y.tsv");
    let want_d = scratch.join("d.tsv");
    let want_c = scratch.join("c.txt");

    let script = format!(
        r#"
import numpy as np
from scipy.cluster.hierarchy import linkage, cophenet
np.random.seed(7)
X = np.random.rand(60, 5)
from scipy.spatial.distance import pdist
Y = pdist(X)
Z = linkage(Y, method="average")
c, zz = cophenet(Z, Y)
with open(r"{z}","w") as f:
    for r in Z: f.write(f"{{int(r[0])}}\t{{int(r[1])}}\t{{r[2]:.17f}}\t{{int(r[3])}}\n")
with open(r"{y}","w") as f:
    for v in Y: f.write(f"{{v:.17f}}\n")
with open(r"{d}","w") as f:
    for v in zz: f.write(f"{{v:.17f}}\n")
with open(r"{cc}","w") as f:
    f.write(f"{{c:.17f}}\n")
"#,
        z = z.display(),
        y = y.display(),
        d = want_d.display(),
        cc = want_c.display(),
    );
    let st = Command::new(&py).args(["-c", &script]).status().unwrap();
    assert!(st.success(), "scipy oracle generation failed");

    let (stdout, stderr) = run(&[z.to_str().unwrap(), "--dist", y.to_str().unwrap()]);

    let ours: Vec<f64> = stdout.lines().map(|l| l.parse().unwrap()).collect();
    let want: Vec<f64> = read(&want_d).lines().map(|l| l.parse().unwrap()).collect();
    assert_eq!(ours.len(), want.len());
    for (a, b) in ours.iter().zip(&want) {
        assert!((a - b).abs() < 1e-9, "cophenetic distance: {a} vs {b}");
    }

    let our_c: f64 = stderr.rsplit(": ").next().unwrap().trim().parse().unwrap();
    let want_c: f64 = read(&want_c).trim().parse().unwrap();
    assert!(
        (our_c - want_c).abs() < 1e-9,
        "correlation: {our_c} vs {want_c}"
    );
}
