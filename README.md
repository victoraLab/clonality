# clonality

`clonality` is an R package for assigning B cell and T cell receptor clonotypes from repertoire sequencing data.
It grew around practical analysis needs for IMGT HighV-QUEST outputs, 10x Genomics Cell Ranger V(D)J annotations, and plate-based repertoire experiments where clone calling and plate quality assessment need to live in the same workflow.

The package currently helps with four broad tasks:

- clonotype assignment from IMGT-style tables
- clonotype assignment from 10x Genomics V(D)J contig annotations
- plate-level contig inspection and QC
- downstream integration with Seurat and immunarch-friendly formats
- interactive Circos-style clone sharing visualization in a browser

## Installation

Core install (pulls all required `Imports` automatically):

```r
# install.packages("devtools")
devtools::install_github("victoraLab/clonality")
```

Recommended install (also pulls the optional `Suggests` packages used by the
Seurat integration, plate-level analysis, and the alluvial purity plot):

```r
devtools::install_github(
  "victoraLab/clonality",
  dependencies = TRUE
)
```

`Suggests` includes:
`Seurat`, `Matrix`, `ggplot2`, `ggalluvial`, `patchwork`, `ggplate`, `purrr`,
`future`, `MASS`, `testthat`. None of these are required for the basic
`clonality()` / `tenx()` workflow, but `add_seurat_metadata()`,
`process_imgt_projects()`, and the plate-level helpers will fail with a
`requireNamespace()` error if their corresponding `Suggests` packages aren't
present. **None of the dependencies are on Bioconductor** — everything is
on CRAN, so no `BiocManager` setup is required.

## What The Package Covers

`clonality` supports:

- IMGT summary and junction outputs
- 10x Genomics Cell Ranger V(D)J annotations, including v3, v7, and v9-style bundled example data
- BCR, TCR, and gamma-delta TCR use cases
- exact or fuzzy grouping based on V gene, J gene, and CDR3 similarity

At a high level, the package separates into a few complementary layers:

1. Read and standardize repertoire inputs.
2. Assign clonotypes from sequence features.
3. Post-process plate-based outputs for QC and visualization.
4. Export or attach results for downstream single-cell analysis.

## Main Workflows

### 1. IMGT-style clonotype assignment

Use `clonality()` when your starting point is an IMGT-derived table or a compatible data frame containing V, J, and CDR3 information.

```r
library(clonality)

head(tra)

tra_out <- clonality(
  data = tra,
  vgene_col = "V_GENE_and_allele",
  jgene_col = "J_GENE_and_allele"
)

# Optionally save to an Excel file
tra_out <- clonality(
  data = tra,
  vgene_col = "V_GENE_and_allele",
  jgene_col = "J_GENE_and_allele",
  file = "tra_output.xlsx"
)
```

Key parameters:

- `mismatch` controls the percent mismatch tolerance (range: 0--100). Internally CDR3 pairwise distances are normalized by length and the clustering tree is cut at `mismatch / 100`. For example, `mismatch = 10` allows up to ~10 % CDR3 divergence within a V-J-CDR3-length bin, while `mismatch = 0` requires identical CDR3 sequences.
- `cdr3_col` accepts one or more CDR3 column names. When multiple columns are provided (e.g. nucleotide and amino-acid), clustering runs independently for each, producing separate `clonality_nt` / `clonality_aa` columns. `mismatch` and `stringdist_method` are recycled to match.
- `file` optionally writes the result to an `.xlsx` file.
- `output_original = TRUE` joins the clonality assignments back to the full input table.

Typical inputs for this workflow are:

- a data frame already loaded in R
- a parsed IMGT table
- package example datasets such as `tra`, `trb`, `IMGT_BCR_Summary_File`, and `IMGT_TCR_Summary_File`

### 2. 10x Genomics clonotype assignment

Use `tenx()` when your starting point is a Cell Ranger V(D)J contig annotation table.
This wrapper handles filtering, grouping by barcode, and delegating clone definition to the package clonality logic.

```r
head(Cellranger9_BCR1)

tenx_out <- tenx(
  data = Cellranger9_BCR1,
  method = "unique_paired",
  cell = "B",
  only_productive = TRUE,
  only_true_cells = TRUE,
  cellranger_version = 9
)
```

Available clonotype strategies include:

- `unique_paired`
- `sticky_ends`
- `unique_all`

Supported cell modes include:

- `cell = "B"`
- `cell = "T"`
- `cell = "Tgd"`

#### Tuning the internal `clonality()` call via `clonality_input`

`tenx()` filters cells, groups them by chain class, and then runs
[`clonality()`](#1-imgt-style-clonotype-assignment) once per chain class to
assign clone IDs. The `clonality_input` argument is a **named `list()`** of
overrides for that internal call.

> **Use `list()`, not `c()`.** A vector flattens nested vectors, so
> `c(cdr3_col = c("CDR3", "cdr3_col2"))` becomes a single character vector
> named `cdr3_col1` / `cdr3_col2` — not the multi-CDR3 spec you intended.
> Always wrap in `list(...)`.

The most useful slots:

| Slot | Default in `tenx()` | What it controls |
|---|---|---|
| `mismatch` | `0` | Percent CDR3 mismatch tolerance (0--100). `10` ≈ 10% divergence. Vectorize alongside `cdr3_col`. |
| `cdr3_col` | `"CDR3"` (nucleotide) | Which concatenated CDR3 column to cluster on. Pass `c("CDR3", "cdr3_col2")` to cluster nucleotide AND amino-acid in one pass, producing `clonality_nt` and `clonality_aa` columns. |
| `stringdist_method` | `"osa"` | String-distance method (`"lv"`, `"dl"`, `"jw"`, ...). Recycled to match `cdr3_col`. |
| `output_original` | `TRUE` | Join clonality back to the per-cell metadata. Set `FALSE` for a minimal table. |
| `file` | `NULL` | Path to write the per-chain-class result as `.xlsx`. |
| `project` | `"project"` | Label written to every row's `project` column — useful when pooling runs. |

Slots managed internally by `tenx()` (rarely worth overriding): `vgene_col`,
`jgene_col`, `ident_col`, `cell`, `suffix`. They are wired to the column names
that `tenx()` constructs in the per-cell metadata.

```r
# Allow 10% CDR3 mismatch
tenx(data, clonality_input = list(mismatch = 10))

# Cluster on nucleotide and amino-acid CDR3 in a single call
tenx(data,
     clonality_input = list(
       cdr3_col = c("CDR3", "cdr3_col2"),
       mismatch = c(0, 10)
     ))

# Save each chain-class clonality table to disk
tenx(data, clonality_input = list(file = "clones.xlsx"))
```

See `?clonality` for the full list of accepted parameters.

### 3. Plate-based QC and visualization

For plate sequencing projects, the package extends beyond clone calling to help identify real contigs, background signals, collapsed wells, and overall plate performance.

```r
output <- clonality(
  data = IMGT_BCR_Summary_File,
  cdr3_col = "DNA_Junction",
  cell = "B",
  mismatch = 90
)

analyzed_contigs <- analyze_well_contigs(
  output,
  barcode_format = "WELL_PLATE"
)

plate_metadata <- extract_plate_metadata(
  df = analyzed_contigs$real_contigs_input,
  barcode_format = "WELL_PLATE"
)

plate_quality <- evaluate_plate_quality(analyzed_contigs)

plot_plate_contig_depth(plate_metadata, plate_id = "P09")
```

![Contig Depth Plot](man/figures/contig_depth_plot1.png)

The same pattern works for TCR plate projects:

```r
output <- clonality(
  data = IMGT_TCR_Summary_File,
  cdr3_col = "DNA_Junction",
  cell = "T"
)

analyzed_contigs <- analyze_well_contigs(
  output,
  barcode_format = "PLATE_WELL"
)

plate_metadata <- extract_plate_metadata(
  df = analyzed_contigs$real_contigs_input,
  barcode_format = "PLATE_WELL"
)

plate_quality <- evaluate_plate_quality(analyzed_contigs)

plot_plate_contig_depth(plate_metadata, plate_id = "P01")
```

![Contig Depth Plot](man/figures/contig_depth_plot2.png)

### 4. Single-cell downstream integration

When clone annotations need to move into single-cell workflows, `clonality` provides direct support for Seurat metadata integration and immunarch-oriented export.

```r
seu <- add_seurat_metadata(
  seu = my_seurat_object,
  clonal.list = ls(pattern = "^Clonal"),
  cell = "T",
  sticky = TRUE,
  overwrite = TRUE
)
```

#### Sticky-ends rescue and per-chain purity (v0.11+)

> For the full explainer including the DNA-strand metaphor that names the
> mode, the biology behind the chain-aware defaults, and an alluvial-plot
> recipe for visualizing rescue purity, see
> [**STICKY_ENDS.md**](STICKY_ENDS.md).

When `sticky = TRUE`, single-chain cells can inherit the dominant paired
clonotype if their CDR3 maps unambiguously into the paired pool. The default
policy encodes the biology of D-gene chains:

```r
purity = list(
  IGH = 0.9, IGK = NA, IGL = NA,    # B cells: rescue heavy only
  TRB = 0.9, TRA = NA,              # αβ T cells: rescue β only
  TRD = 0.9, TRG = NA               # γδ T cells: rescue δ only
)
```

D-gene chains (IGH, TRB, TRD) have hyper-variable CDR3s that uniquely identify
a clone, so rescue at 0.9 purity is reliable. Non-D chains (IGK, IGL, TRA,
TRG) have constrained CDR3 length and limited V repertoires, so the same CDR3
typically appears across many unrelated paired clones — `NA` skips them.

Override per-chain to opt into κ/λ rescue at high stringency:

```r
seu <- add_seurat_metadata(
  seu, sticky = TRUE,
  purity = list(IGH = 0.9, IGK = 0.95, IGL = 0.95),
  verbose = TRUE   # prints inspect_purity() summary after the rescue pass
)
```

After the rescue, the metadata gains `Res_dominant_purity` (numeric in
[0, 1]) suitable for alluvial visualization:

```r
library(ggplot2); library(ggalluvial)
df <- seu[[]] %>% dplyr::filter(!is.na(Res_dominant_purity))
ggplot(df, aes(axis1 = clonality, axis2 = Res_clonality, alpha = Res_dominant_purity)) +
  geom_alluvium(aes(fill = Res_clonality.correction)) +
  geom_stratum() +
  theme_minimal()
```

To audit the rescue at any time:

```r
inspect_purity(seu)
#> Sticky-ends rescue summary
#>  chain singles matched median_purity rescued at_0.8 at_0.9 at_1
#>    IGH    4477    2801         1.000    2264    2264   1750 1576
#>    IGK    6134    5116         0.395       0       0      0    0
#>    IGL     504     461         0.211       0       0      0    0
```

**Migration from < v0.11**: scalar `purity = 0.8` still works (with a message)
but is no longer recommended. The match logic also defaults to `match_mode =
"exact"`; pass `match_mode = "substring"` to recover the legacy behavior.

```r
imm <- gen_immunarch(
  seu,
  proportion = "rm.na",
  sample_column = "Conditions"
)
```

### 5. Interactive Circos viewer

The package also bundles a lightweight D3 viewer for clone sharing across
conditions or tissues.

Open the viewer with:

```r
clonality::open_circos_viewer()
```

If your data are already in a wide clone-by-group table:

```r
wide_df <- data.frame(
  Clone_ID = c("1.1", "2.1", "6.1"),
  M1LP_PC = c(2, 3, 3),
  M1PP_GC = c(11, 0, 0),
  M2LP_PC = c(0, 0, 1)
)

viewer_df <- clonality::circos_viewer_wide_to_long(
  wide_df,
  clone_col = "Clone_ID"
)
```

If your table already has one row per cell or per contig plus a grouping
column:

```r
results$condition <- sample(
  c("M1LP_PC", "M1PP_GC", "M2LP_PC"),
  nrow(results),
  replace = TRUE
)

viewer_df <- clonality::circos_viewer_input(
  data = results,
  clonotype_col = "clonality",
  group_col = "condition"
)
```

Then paste `viewer_df` into the viewer using `Clonotype table` mode.

## Function Map

One useful way to navigate the package is by grouping functions according to the step of the analysis they support.

### Core clone assignment

| Function | Purpose |
| --- | --- |
| `clonality()` | Main clone assignment engine for IMGT-style or compatible repertoire tables |
| `cluster_duplicates_by_cdr3()` | Subclusters duplicated V-J-CDR3-length bins by CDR3 similarity |
| `assign_unique_ids()` | Assigns unique clone IDs when rows do not belong to duplicated bins |
| `simplify_genenames()` | Simplifies gene naming to make grouping more robust |

### 10x Genomics processing

| Function | Purpose |
| --- | --- |
| `tenx()` | End-to-end 10x V(D)J clonotype annotation wrapper |
| `assigntenx()` | Internal-style lower-level assignment helper for 10x workflows |
| `cellranger_version_parser()` | Detects or validates Cell Ranger annotation schema version |
| `extract_and_unite()` | Pulls chain-specific columns and optionally collapses them row-wise |

### Input parsing and project ingestion

| Function | Purpose |
| --- | --- |
| `parse_general()` | Reads general tabular inputs such as `.csv`, `.txt`, and `.xlsx` |
| `parse_imgt()` | Parses IMGT `1_Summary.txt` and `6_Junction.txt` outputs into a merged table |
| `filter_sequences_imgt()` | Filters IMGT-derived sequences before downstream processing |
| `process_imgt_projects()` | Runs multiple IMGT project folders in parallel and combines outputs |

### Plate-level analysis and QC

| Function | Purpose |
| --- | --- |
| `analyze_well_contigs()` | Classifies well contigs into real and background categories |
| `extract_plate_metadata()` | Parses plate, well, contig, and depth metadata from cell IDs |
| `evaluate_plate_quality()` | Summarizes plate efficiency and collapse metrics |
| `plot_plate_contig_depth()` | Displays well-level contig depth for one plate |
| `plot_multiple_plates()` | Arranges several plate-level depth plots together |

### Downstream integration

| Function | Purpose |
| --- | --- |
| `add_seurat_metadata()` | Writes clonality annotations into a Seurat object, with chain-aware sticky-ends rescue |
| `inspect_purity()` | Per-chain summary of sticky-ends rescue purity (singles, matched, median purity, rescued) |
| `normalize_airr_to_cellranger()` | Preprocesses IgBLAST AIRR / nanopore contigs into a CellRanger-shape table |
| `gen_immunarch()` | Builds an immunarch-compatible structure from Seurat metadata |

## Example Pipelines

### IMGT BCR example

```r
output <- clonality(
  data = IMGT_BCR_Summary_File,
  cdr3_col = "DNA_Junction",
  cell = "B",
  mismatch = 90
)
```

### IMGT TCR example

```r
output <- clonality(
  data = IMGT_TCR_Summary_File,
  cdr3_col = "DNA_Junction",
  cell = "T"
)
```

### Cell Ranger v7 BCR example

```r
head(Cellranger7_BCR1)

tenx_out <- tenx(
  data = Cellranger7_BCR1,
  method = "unique_paired",
  cell = "B",
  only_productive = TRUE,
  only_true_cells = TRUE,
  clonality_input = list(mismatch = 10),
  cellranger_version = 7
)
```

### Cell Ranger v3 TCR example

```r
head(Cellranger3_TCR)

tenx_out <- tenx(
  data = Cellranger3_TCR,
  method = "sticky_ends",
  cell = "T",
  only_productive = TRUE,
  cellranger_version = 3
)
```

### Custom data -- single CDR3 column

```r
result <- clonality(
  data      = formattedData,
  ident_col = "SequenceID",
  vgene_col = "V_gene",
  jgene_col = "J_gene",
  cell      = "B",
  cdr3_col  = "Junction",
  mismatch  = 10,
  output_original = TRUE,
  file      = "clonality_10.xlsx"
)
```

### Vector `cdr3_col` and/or `mismatch` (v0.12+)

`cdr3_col` and `mismatch` accept vectors. Whichever is longer determines the
number of independent clustering passes; the shorter argument recycles. Each
pass produces one `clonality_<tag>` column.

```r
# (a) multi-CDR3: cluster nt and aa CDR3s side by side
clonality(df, cdr3_col = c("Junction_nt", "Junction_aa"), mismatch = c(0, 10))
# → clonality_nt, clonality_aa

# (b) stringency scan: same CDR3, three thresholds
clonality(df, cdr3_col = "Junction", mismatch = c(0, 10, 20))
# → clonality_mm0, clonality_mm10, clonality_mm20

# (c) one threshold across two CDR3 columns
clonality(df, cdr3_col = c("Junction_nt", "Junction_aa"), mismatch = 0)
# → clonality_nt, clonality_aa (both at mm0)
```

Tags are auto-detected as `"nt"`/`"aa"` from sequence content; when the same
CDR3 column is reused, they fall back to `"mm<value>"`. Pass
`cdr3_tags = c("strict", "loose")` to label them yourself.

`stringdist_method` and `suffix` follow the same recycling rule.

The same vector `mismatch` works through `tenx()`:

```r
tenx(data = "...filtered_contig_annotations.csv",
     method = "sticky_ends", cell = "B",
     clonality_input = list(mismatch = c(0, 10, 20)))
# Each Clonal.output.10x* table now has clonality_mm0, clonality_mm10, clonality_mm20
```

### Multi-project IMGT processing

```r
results <- process_imgt_projects(
  project_paths = c("project_1", "project_2", "project_3"),
  cell_type = "B",
  mismatch = 10
)
```

## Included Example Data

The package ships with example objects that make it easier to test workflows and document expected inputs:

| Dataset | Description |
| --- | --- |
| `tra` | Small IMGT TCR alpha dataset |
| `trb` | Small IMGT TCR beta dataset |
| `IMGT_BCR_Summary_File` | Processed IMGT HighV-QUEST BCR summary table |
| `IMGT_TCR_Summary_File` | Processed IMGT HighV-QUEST TCR summary table |
| `Cellranger3_TCR` | Example 10x Cell Ranger v3 TCR contig annotations |
| `Cellranger7_BCR1`, `Cellranger7_BCR2` | Example 10x Cell Ranger v7 BCR contig annotations |
| `Cellranger9_BCR1` | Example 10x Cell Ranger v9 BCR contig annotations |
| `Cellranger7_TCRgd` | Example 10x Cell Ranger v7 gamma-delta TCR contig annotations |
| `chain_classes` | Supported chain class definitions used in 10x workflows |

## Suggested Reading Order

If you are new to the package, a good progression is:

1. Start with `clonality()` if your data already looks IMGT-like.
2. Use `tenx()` if your input comes from Cell Ranger V(D)J output.
3. Add `analyze_well_contigs()`, `extract_plate_metadata()`, and `evaluate_plate_quality()` for plate-based projects.
4. Finish with `add_seurat_metadata()` or `gen_immunarch()` if you need downstream single-cell integration.

## Current Focus Of The Package

The package reflects several years of iterative development around real experimental use cases.
That is a strength, because it contains practical workflows, but it also means some parts of the interface are richer than the README previously suggested.
This document is intended to make the current surface area easier to discover before the next pass on code organization, consistency, and efficiency.
