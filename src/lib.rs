use std::io::{BufRead, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// One merge in a scipy-style linkage matrix Z: children `left`/`right`
/// (ids < n are original observations, ids >= n refer to merge `id - n`),
/// joined at cophenetic `height`, forming a cluster of `size` leaves.
#[derive(Debug, Clone, Copy)]
pub struct Merge {
    pub left: usize,
    pub right: usize,
    pub height: f64,
    pub size: usize,
}

pub fn read_linkage(path: &Path) -> Result<Vec<Merge>> {
    let file = std::fs::File::open(path)?;
    let mut merges = Vec::new();
    for (lineno, line) in std::io::BufReader::new(file).lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 4 {
            return Err(RsomicsError::InvalidInput(format!(
                "linkage line {} has {} columns, expected 4 (left right height size)",
                lineno + 1,
                cols.len()
            )));
        }
        merges.push(Merge {
            left: parse_id(cols[0])?,
            right: parse_id(cols[1])?,
            height: cols[2].parse()?,
            size: parse_id(cols[3])?,
        });
    }
    if merges.is_empty() {
        return Err(RsomicsError::InvalidInput("empty linkage matrix".into()));
    }
    validate_linkage(&merges)?;
    Ok(merges)
}

/// Reject linkage matrices scipy's `is_valid_linkage` rejects, in scipy's order.
/// Every branch here would otherwise index out of bounds or silently emit a
/// wrong cophenetic vector in `cophenetic_distances`. Child ids are already
/// non-negative integers (`parse_id`), so only the height, count, ordering and
/// uniqueness invariants remain to enforce.
pub fn validate_linkage(merges: &[Merge]) -> Result<()> {
    let num_obs = merges.len() + 1;

    for m in merges {
        if m.height < 0.0 {
            return Err(RsomicsError::InvalidInput(
                "linkage contains negative distances".into(),
            ));
        }
        if m.size > num_obs {
            return Err(RsomicsError::InvalidInput(
                "linkage contains excessive observations in a cluster".into(),
            ));
        }
    }

    for (i, m) in merges.iter().enumerate() {
        if m.left >= num_obs + i || m.right >= num_obs + i {
            return Err(RsomicsError::InvalidInput(
                "linkage uses non-singleton cluster before it is formed".into(),
            ));
        }
    }

    let mut seen = std::collections::HashSet::with_capacity(2 * merges.len());
    for m in merges {
        if !seen.insert(m.left) || !seen.insert(m.right) {
            return Err(RsomicsError::InvalidInput(
                "linkage uses the same cluster more than once".into(),
            ));
        }
    }

    Ok(())
}

fn parse_id(s: &str) -> Result<usize> {
    let v: f64 = s.parse()?;
    if v < 0.0 || v.fract() != 0.0 {
        return Err(RsomicsError::InvalidInput(format!(
            "expected non-negative integer, got {s}"
        )));
    }
    Ok(v as usize)
}

#[inline]
fn condensed_index(n: usize, i: usize, j: usize) -> usize {
    let (i, j) = if i < j { (i, j) } else { (j, i) };
    n * i - (i * (i + 1)) / 2 + (j - i - 1)
}

/// Condensed cophenetic distance vector: entry for pair (i,j) is the merge
/// height at which i and j first share a cluster. Iterative post-order DFS
/// over the dendrogram, matching scipy `_hierarchy.pyx::cophenetic_distances`.
pub fn cophenetic_distances(merges: &[Merge]) -> Vec<f64> {
    let n = merges.len() + 1;
    let mut d = vec![0.0f64; n * (n - 1) / 2];

    let mut curr_node = vec![0usize; n];
    let mut left_start = vec![0usize; n];
    let mut members = vec![0usize; n];
    let mut visited = vec![false; 2 * n - 1];

    let mut k: isize = 0;
    curr_node[0] = 2 * n - 2;
    left_start[0] = 0;

    while k >= 0 {
        let ku = k as usize;
        let root = curr_node[ku] - n;
        let i_lc = merges[root].left;
        let i_rc = merges[root].right;

        let n_lc;
        if i_lc >= n {
            n_lc = merges[i_lc - n].size;
            if !visited[i_lc] {
                visited[i_lc] = true;
                k += 1;
                curr_node[k as usize] = i_lc;
                left_start[k as usize] = left_start[ku];
                continue;
            }
        } else {
            n_lc = 1;
            members[left_start[ku]] = i_lc;
        }

        let n_rc;
        if i_rc >= n {
            n_rc = merges[i_rc - n].size;
            if !visited[i_rc] {
                visited[i_rc] = true;
                k += 1;
                curr_node[k as usize] = i_rc;
                left_start[k as usize] = left_start[ku] + n_lc;
                continue;
            }
        } else {
            n_rc = 1;
            members[left_start[ku] + n_lc] = i_rc;
        }

        let dist = merges[root].height;
        let right_start = left_start[ku] + n_lc;
        for &mi in &members[left_start[ku]..right_start] {
            for &mj in &members[right_start..right_start + n_rc] {
                d[condensed_index(n, mi, mj)] = dist;
            }
        }

        k -= 1;
    }

    d
}

/// Pearson r between the condensed cophenetic vector and the original
/// condensed distances Y. Same arithmetic order as scipy's wrapper.
pub fn cophenetic_correlation(coph: &[f64], y: &[f64]) -> Result<f64> {
    if coph.len() != y.len() {
        return Err(RsomicsError::InvalidInput(format!(
            "distance vector has {} entries, cophenetic has {}",
            y.len(),
            coph.len()
        )));
    }
    let m = coph.len() as f64;
    let z_mean = coph.iter().sum::<f64>() / m;
    let y_mean = y.iter().sum::<f64>() / m;
    let (mut num, mut da, mut db) = (0.0f64, 0.0f64, 0.0f64);
    for (&zi, &yi) in coph.iter().zip(y) {
        let yy = yi - y_mean;
        let zz = zi - z_mean;
        num += yy * zz;
        da += yy * yy;
        db += zz * zz;
    }
    Ok(num / (da * db).sqrt())
}

/// Read original distances Y as either a condensed vector (one value per line,
/// or whitespace/tab separated) or a square matrix, returning the condensed form.
pub fn read_distances(path: &Path, square: bool, n: usize) -> Result<Vec<f64>> {
    let text = std::fs::read_to_string(path)?;
    let vals: Vec<f64> = text
        .split(|c: char| c.is_whitespace())
        .filter(|t| !t.is_empty())
        .map(|t| t.parse::<f64>())
        .collect::<std::result::Result<_, _>>()?;

    if square {
        if vals.len() != n * n {
            return Err(RsomicsError::InvalidInput(format!(
                "square distance matrix has {} values, expected {}x{}={}",
                vals.len(),
                n,
                n,
                n * n
            )));
        }
        let mut y = Vec::with_capacity(n * (n - 1) / 2);
        for i in 0..n {
            for j in (i + 1)..n {
                y.push(vals[i * n + j]);
            }
        }
        Ok(y)
    } else {
        let expected = n * (n - 1) / 2;
        if vals.len() != expected {
            return Err(RsomicsError::InvalidInput(format!(
                "condensed distance vector has {} values, expected {expected} for {n} observations",
                vals.len()
            )));
        }
        Ok(vals)
    }
}

use std::fmt::Write as _;

pub fn write_condensed(d: &[f64], out: &mut dyn Write) -> Result<()> {
    let mut w = std::io::BufWriter::new(out);
    let mut buf = String::with_capacity(1 << 16);
    for v in d {
        let _ = writeln!(buf, "{v:.17}");
        if buf.len() >= (1 << 16) - 32 {
            w.write_all(buf.as_bytes())?;
            buf.clear();
        }
    }
    w.write_all(buf.as_bytes())?;
    w.flush()?;
    Ok(())
}

pub fn write_square(d: &[f64], n: usize, out: &mut dyn Write) -> Result<()> {
    let mut w = std::io::BufWriter::new(out);
    let mut buf = String::with_capacity(1 << 16);
    for i in 0..n {
        for j in 0..n {
            if j > 0 {
                buf.push('\t');
            }
            if i == j {
                buf.push_str("0.00000000000000000");
            } else {
                let _ = write!(buf, "{:.17}", d[condensed_index(n, i, j)]);
            }
        }
        buf.push('\n');
        if buf.len() >= (1 << 16) - 1024 {
            w.write_all(buf.as_bytes())?;
            buf.clear();
        }
    }
    w.write_all(buf.as_bytes())?;
    w.flush()?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    linkage: &Path,
    dist: Option<&Path>,
    dist_is_square: bool,
    square: bool,
    out: &mut dyn Write,
) -> Result<()> {
    let merges = read_linkage(linkage)?;
    let n = merges.len() + 1;
    let coph = cophenetic_distances(&merges);

    if let Some(dpath) = dist {
        let y = read_distances(dpath, dist_is_square, n)?;
        let c = cophenetic_correlation(&coph, &y)?;
        eprintln!("cophenetic correlation coefficient: {c:.17}");
    }

    if square {
        write_square(&coph, n, out)
    } else {
        write_condensed(&coph, out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn single_link_example() -> Vec<Merge> {
        let rows: &[(usize, usize, f64, usize)] = &[
            (0, 1, 1.0, 2),
            (2, 12, 1.0, 3),
            (3, 4, 1.0, 2),
            (5, 14, 1.0, 3),
            (6, 7, 1.0, 2),
            (8, 16, 1.0, 3),
            (9, 10, 1.0, 2),
            (11, 18, 1.0, 3),
            (13, 15, 2.0, 6),
            (17, 20, 2.0, 9),
            (19, 21, 2.0, 12),
        ];
        rows.iter()
            .map(|&(l, r, h, s)| Merge {
                left: l,
                right: r,
                height: h,
                size: s,
            })
            .collect()
    }

    #[test]
    fn matches_scipy_doc_example() {
        let expected = [
            1., 1., 2., 2., 2., 2., 2., 2., 2., 2., 2., 1., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
            2., 2., 2., 2., 2., 2., 2., 2., 1., 1., 2., 2., 2., 2., 2., 2., 1., 2., 2., 2., 2., 2.,
            2., 2., 2., 2., 2., 2., 2., 1., 1., 2., 2., 2., 1., 2., 2., 2., 2., 2., 2., 1., 1., 1.,
        ];
        let d = cophenetic_distances(&single_link_example());
        assert_eq!(d.len(), expected.len());
        for (got, want) in d.iter().zip(expected) {
            assert_eq!(*got, want);
        }
    }

    #[test]
    fn condensed_index_roundtrip() {
        let n = 12;
        let mut k = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                assert_eq!(condensed_index(n, i, j), k);
                assert_eq!(condensed_index(n, j, i), k);
                k += 1;
            }
        }
    }

    #[test]
    fn correlation_perfect_for_self() {
        let d = cophenetic_distances(&single_link_example());
        let c = cophenetic_correlation(&d, &d).unwrap();
        assert!((c - 1.0).abs() < 1e-12);
    }
}
