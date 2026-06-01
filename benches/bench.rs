use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_cophenet::{Merge, cophenetic_distances};

/// A balanced binary dendrogram over `n` leaves, heights increasing by level.
fn balanced_linkage(n: usize) -> Vec<Merge> {
    let mut merges = Vec::with_capacity(n - 1);
    let mut next = n;
    let mut frontier: Vec<(usize, usize)> = (0..n).map(|i| (i, 1)).collect();
    let mut h = 1.0;
    while frontier.len() > 1 {
        let mut nxt = Vec::with_capacity(frontier.len().div_ceil(2));
        let mut it = frontier.into_iter();
        while let Some((a, sa)) = it.next() {
            if let Some((b, sb)) = it.next() {
                merges.push(Merge {
                    left: a,
                    right: b,
                    height: h,
                    size: sa + sb,
                });
                nxt.push((next, sa + sb));
                next += 1;
            } else {
                nxt.push((a, sa));
            }
        }
        frontier = nxt;
        h += 1.0;
    }
    merges
}

fn bench(c: &mut Criterion) {
    let merges = balanced_linkage(3000);
    c.bench_function("cophenetic_distances_3000", |b| {
        b.iter(|| cophenetic_distances(black_box(&merges)))
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
