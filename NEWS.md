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
