# rsomics-cophenet

Cophenetic distances and the cophenetic correlation coefficient from a
hierarchical-clustering linkage matrix.

The cophenetic distance between two original observations is the height of the
linkage at which they first fall into the same cluster. Given a linkage matrix
`Z` (the scipy / `rsomics-upgma` format), this emits the condensed cophenetic
distance vector. When the original condensed distance matrix is also supplied,
it additionally reports the cophenetic correlation coefficient — the Pearson
correlation between the condensed cophenetic distances and the original
distances, a standard goodness-of-fit measure for a dendrogram.

```
rsomics-cophenet <linkage.tsv> [--dist orig.tsv [--dist-square]] [--square] [-o cophenetic.tsv]
```

- `linkage.tsv` — a linkage matrix, one merge per line: `left  right  height  size`
  (tab-separated). Ids below `n` are original observations, ids `>= n` reference
  merge `id - n`. This is exactly the matrix `rsomics-upgma --linkage` writes.
- `--dist orig.tsv` — the original distances `Z` was built from. Enables the
  cophenetic correlation coefficient, printed to stderr. Condensed by default
  (the upper triangle, row-major); pass `--dist-square` for a full square matrix.
- `--square` — emit a full square cophenetic matrix instead of the condensed
  vector.
- `-o` — output path (`-` or omitted = stdout).

Companion to [`rsomics-upgma`](https://github.com/omics-rust/rsomics-upgma),
which produces the linkage matrix this tool consumes.

## Origin

This crate is an independent Rust reimplementation of
`scipy.cluster.hierarchy.cophenet`, based on:

- SciPy's `scipy/cluster/hierarchy/_hierarchy.pyx` `cophenetic_distances`
  (the iterative post-order dendrogram traversal and the condensed-index
  formula) and the `cophenet` Python wrapper (the Pearson-correlation
  arithmetic), both BSD-3-Clause, read and cited.

Output is value-exact (to ~1e-9, in practice bit-identical for the distances
and to f64 rounding for the correlation) against scipy `cophenet` on the same
linkage matrix.

License: MIT OR Apache-2.0.
Upstream credit: SciPy (https://scipy.org, BSD-3-Clause).
