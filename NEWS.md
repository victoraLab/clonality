# clonality 0.14

## Per-chain independent Hamming in `tenx()` / `clonality()`

`tenx(clonality_input = list(max_edits_per_chain = 1, stringdist_method = "hamming"))`
now enforces an **independent** raw Hamming-1 cutoff on every chain separately.
Two cells cluster together only when *every* chain's CDR3 is within 1 edit of its
counterpart — not 1 edit on the concatenated string.

Internally, `assigntenx()` intercepts `max_edits_per_chain`, auto-detects the
`cdr3_nt_*` per-chain columns in the data frame, and passes them as
`per_chain_cdr3_cols` to `clonality()`. The new `per_chain_cdr3_cols` parameter
is also available directly in `clonality()` for non-10x workflows.

The distance between two cells is `max(hamming(chain_k_a, chain_k_b) for each k)`;
complete-linkage hclust cut at `h = max_edits` ensures within-cluster consistency.

```r
tenx(
  data            = paths,
  project_names   = project_names,
  method          = "sticky_ends",
  cell            = "B",
  clonality_input = list(max_edits_per_chain = 1, stringdist_method = "hamming")
)
```

# clonality 0.13

## Multi-sample ingestion in `tenx()`

`tenx()` now accepts a list or character vector of multiple inputs in one call.
Each sample is loaded and filtered independently; barcodes are suffixed with the
corresponding `project_names` entry so they remain unique across samples.
Clonotype assignment then runs on the combined pool, so cells sharing a CDR3
across samples naturally receive the same clone ID.

```r
out <- tenx(
  data          = list("sample1.csv", "sample2.csv", "sample3.csv"),
  project_names = c("d7", "d20", "boost1"),
  method        = "sticky_ends",
  cell          = "B"
)
# Barcodes: "AAACCTGA-1_d7", "AAACCTGA-1_d20", ...
```

The single-sample API is unchanged. Internally the load → normalize → filter
logic was extracted into a `.tenx_prepare_sample()` helper called by both paths.

# clonality 0.12

## Bug fix: per-chain CDR3 length binning in paired clones (via `tenx()`)

`tenx()` ultimately calls `clonality()` with a single concatenated CDR3 string
per cell (e.g. `"<heavy_cdr3>_<light_cdr3>"`). The bin key inside `clonality()`
included the **total** concatenated CDR3 length but not the per-chain lengths.
As a result, two paired cells with compensating cross-chain indels — e.g.
`cdr3_IGH = 12, cdr3_IGK = 10` and `cdr3_IGH = 13, cdr3_IGK = 9`, both with
total length 22 — ended up in the same bin and could merge under OSA
sub-clustering at `mismatch > 0`.

Fixed by adding an `extra_bin_cols` argument to `clonality()` and setting the
default in `assigntenx()` to `"cdr3_length"`, the per-chain length string
already produced for the output (e.g. `"12_10"`). This makes paired cells
that disagree on per-chain CDR3 length land in separate bins, regardless of
mismatch tolerance.

Behavioral change: at `mismatch > 0`, paired clones that used to absorb
cross-chain-indel cells now split them into separate clonality IDs. Single-
chain workflows are unaffected.

## Mismatch-only scans in `clonality()` (and through `tenx()`)

`clonality()` now accepts a single `cdr3_col` with a vector `mismatch`, running
one independent clustering pass per threshold and producing one
`clonality_<tag>` column per pass:

```r
clonality(data = df, cdr3_col = "Junction", mismatch = c(0, 10, 20))
# → clonality_mm0, clonality_mm10, clonality_mm20
```

The recycling rule is symmetric: `K = max(length(cdr3_col), length(mismatch))`,
and whichever argument is length 1 is recycled to match. So all three of these
work in a single call:

| Call | K | Output columns |
|---|---|---|
| `cdr3_col = "Junction", mismatch = 0` | 1 | `clonality` (legacy single-pass) |
| `cdr3_col = c("nt","aa"), mismatch = c(0, 10)` | 2 | `clonality_nt`, `clonality_aa` |
| `cdr3_col = "Junction", mismatch = c(0, 10, 20)` | 3 | `clonality_mm0`, `clonality_mm10`, `clonality_mm20` |
| `cdr3_col = c("nt","aa"), mismatch = 0` | 2 | `clonality_nt`, `clonality_aa` (both at mm0) |

Tag auto-naming falls back to `"mm<value>"` when the same `cdr3_col` is used
across multiple passes; pass `cdr3_tags = c(...)` to override.

The same vector `mismatch` flows through `tenx(clonality_input = list(mismatch = c(0, 10, 20)))`
to produce a stringency scan on a 10x dataset in one call.

## Other

- Simplified the `?clonality` docstring; the long-form explanation moved to NEWS.

# clonality 0.11

## Sticky-ends rescue overhaul

The single-chain rescue logic in `add_seurat_metadata(sticky = TRUE)` has been
reworked to encode the biology of D-gene vs non-D-gene chains and to surface
per-cell purity scores. A purity analysis on a low-paired-yield nanopore B-cell
dataset showed that the previous defaults rescued IGH cells correctly (median
purity 1.0) but fabricated clonal expansions for IGK / IGL cells (median
purity 0.40 / 0.21).

### Breaking changes (correctness-positive)

All of these reduce the number of rescues but make the surviving ones
trustworthy. Every change has an opt-out path.

- **`purity` is now a per-chain named list.**
  Default: `list(IGH = 0.9, IGK = NA, IGL = NA, TRB = 0.9, TRA = NA, TRD = 0.9, TRG = NA)`.
  D-gene chains rescue at 0.9; non-D chains skip rescue entirely (`NA` ⇒
  `Res_clonality.correction = "Not tested (chain skipped)"`). A scalar value
  still works for back-compat with a one-time message recommending the
  chain-aware defaults.

- **`match_mode` defaults to `"exact"`.**
  The legacy `grepl(fixed = TRUE)` substring match systematically inflated the
  paired-pool match set with length-mismatched neighbors. Pass
  `match_mode = "substring"` to recover the old behavior (intended only for
  truncated-read pipelines).

- **`require_v_gene = TRUE` by default.**
  Paired matches are additionally filtered to share the V gene with the
  single-chain cell when V-gene columns can be resolved. Falls through
  silently when V-gene columns are absent.

### New features

- **`Res_dominant_purity` column** added to the Seurat metadata after a sticky
  pass. Numeric in `[0, 1]` (or `NA` if not tested) representing the fraction
  of paired-pool match-set rows belonging to the dominant clone. Useful for
  alluvial plots of single → paired clone with ribbon alpha = purity.

- **`min_cdr3_aa = 9` floor** on the single-chain CDR3 length. Cells with
  shorter CDR3s are reported as `"Not tested (cdr3 too short)"` and never
  attempt a rescue.

- **`inspect_purity()` helper** prints a per-chain summary table (singles,
  matched, median purity, rescued at thresholds 0.8 / 0.9 / 1.0). Auto-called
  when `add_seurat_metadata(verbose = TRUE)`.

- **`normalize_airr_to_cellranger()` helper** preprocesses IgBLAST AIRR / nanopore
  contigs into a CellRanger-shape table: collapses comma-separated alleles,
  fills C-genes from an external column, coerces the `productive` flag to
  `"true"`/`"false"` strings.

- **`tenx(airr_normalize = "auto")`** auto-detects AIRR-style input and applies
  the normalization above before column renaming.

### Migration

| Previous call | v0.11 equivalent |
|---|---|
| `add_seurat_metadata(seu, sticky = TRUE)` | Same; non-D chains now skipped automatically. Pass `purity = list(IGH = 0.8, IGK = 0.95, IGL = 0.95)` to rescue κ/λ at very high purity. |
| `add_seurat_metadata(seu, sticky = TRUE, purity = 0.8)` | Works with a back-compat message; consider `purity = NULL` to use the chain-aware defaults. |
| Substring rescue desired | `match_mode = "substring"` |
| V-gene check unwanted | `require_v_gene = FALSE` |

# clonality 0.10

- Restored multi-CDR3 / multi-mismatch support in `clonality()` — `cdr3_col` and
  `mismatch` accept vectors and produce one `clonality_<tag>` column per CDR3.
- Updated the bundled circos viewer (sharing matrix + clonotype modes) to
  match the Observable D3 chord-diagram style: slim arc bands, perimeter tick
  marks + numeric labels, white ribbon strokes.
- Fixed `tenx()` "subscript out of bounds" crash when `add_columns` duplicated a
  base column (`row1.R`) and when `clonality_input$cdr3_col` was a vector.
