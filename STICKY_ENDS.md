# Sticky-ends rescue and per-chain purity (v0.11+)

## The basic idea

In a 10x BCR/TCR experiment we recover **paired** cells (both heavy and light
chain captured) and **single-chain** cells (only one chain captured, usually
because the other failed library prep or didn't amplify). Canonical
`tenx(method = "unique_paired")` clones are defined from the paired pool only,
so the single-chain cells are dropped — even though many of them probably
belong to clones that *are* present in the paired pool.

The package borrows a metaphor from molecular biology to describe how to
recover them: think of paired cells as **double-stranded DNA** and single-chain
cells as **single strands** with a **sticky end** looking for a complement.

```
                                clone 1.5
   paired cells                 (IGH_IGK_1.5)
   ────────────────────────────────────────────
   IGH-CARWGLRQEGFD       IGH-CARWGLRQEGFD       <-- heavy strand
       |||||||||||           |||||||||||
   IGK-CFQGSHVPWTF        IGK-CFQGSHVPWTF        <-- light strand


   single-chain cells
   ───────────────────
   IGH-CARWGLRQEGFD                              <-- sticky end (heavy only)
       ?
       (looking for a partner in the paired pool)
```

If the heavy CDR3 of a single-chain cell **uniquely** matches one paired clone
(here `1.5`), the package "anneals" the single strand into that clone:
`single_x → clone 1.5`.

If the same heavy CDR3 matches *several* unrelated paired clones, the sticky
end is **promiscuous** — there is no unambiguous partner — and the rescue is
denied to avoid fabricating clonal expansions.

The metric that captures "how unambiguous is the match" is the
**dominant purity score**:

```
                  count of paired-pool rows mapping to the dominant clone
   purity   =   ─────────────────────────────────────────────────────────
                  count of paired-pool rows that match this single's CDR3
```

- **purity = 1.0** → all matches go to one paired clone (true 1-to-1 anneal)
- **purity = 0.5** → matches split equally between two paired clones (ambiguous)
- **purity = 0.2** → matches scattered across 5+ clones (don't rescue)

## The biology: why some chains anneal cleanly and others don't

Not all chains are equally informative. The heuristic that drives v0.11's
defaults:

```
   D-gene chains          ←  hyper-variable CDR3, uniquely fingerprints a clone
   (IGH, TRB, TRD)             → median purity ≈ 1.0, safe to rescue

   Non-D chains           ←  short CDR3, limited V repertoire, recurs
   (IGK, IGL, TRA, TRG)       across many unrelated clones
                              → median purity ≈ 0.2-0.4, rescue is unreliable
```

In a typical mouse BCR dataset, an IGH CDR3 like `CARWGLRQEGFDYW` has ~10⁹
possible sequences from random VH·D·JH joining + N-additions, so the *one*
paired clone that carries it is almost certainly the right partner. A κ CDR3
like `CFQGSHVPWTF` has no D segment, limited N-additions, and only ~140
functional Vκ — so the same κ CDR3 typically appears in dozens of unrelated
heavy-chain clones. Annealing to the paired pool by κ alone fabricates clonal
expansions.

The default `purity` policy in v0.11 encodes this asymmetry directly:

```r
purity = list(
  IGH = 0.9, IGK = NA, IGL = NA,    # B cells: rescue heavy only
  TRB = 0.9, TRA = NA,              # αβ T cells: rescue β only
  TRD = 0.9, TRG = NA               # γδ T cells: rescue δ only
)
```

`NA` means "skip this chain entirely" — non-D singles stay as single-chain
entries with `Res_clonality.correction = "Not tested (chain skipped)"`.

You can override per chain to opt into κ/λ rescue at very high stringency
(only well-justified for projects where you've confirmed the light-chain
repertoire is permissive enough):

```r
purity = list(IGH = 0.9, IGK = 0.95, IGL = 0.95)
```

## Three additional safety checks

On top of the per-chain purity threshold, v0.11 layers three defensive checks:

### 1. Exact CDR3 match (`match_mode = "exact"`)

Pre-v0.11 the rescue used `grepl(fixed = TRUE)` substring matching — a 9-aa
single CDR3 could match any 11-aa paired CDR3 with that prefix or suffix.
That inflated the match set with length-mismatched neighbors that aren't even
the same clone family. v0.11 defaults to exact equality. Pass
`match_mode = "substring"` to recover the legacy behavior (only useful for
truncated-read pipelines).

### 2. V-gene consistency (`require_v_gene = TRUE`)

Even after CDR3 match, the rescue additionally requires the single and
matched paired rows to share the same V gene. True clonal relatives always do.
This filter is silently skipped when V-gene columns aren't resolvable on both
sides.

### 3. CDR3 length floor (`min_cdr3_aa = 9`)

Singles with CDR3s shorter than 9 aa are reported as
`"Not tested (cdr3 too short)"` and never attempt a rescue. Very short CDR3s
are dominated by germline-encoded residues and substring-match indiscriminately.

## Output columns

After `add_seurat_metadata(sticky = TRUE)` you get four `Res_*` columns on
top of the canonical clonality:

| Column | Meaning |
|---|---|
| `Res_clonality` | The (possibly rescued) clonality ID. Equals `Res_clonality_paired` when rescue succeeded. |
| `Res_clonality.correction` | One of `"Normal"`, `"Corrected"`, `"Not tested"`, `"Not tested (chain skipped)"`, `"Not tested (cdr3 too short)"`, `"Not tested (no match)"` |
| `Res_dominant_purity` | Numeric in [0, 1], or `NA` if not tested. Plot this to audit your rescue. |
| `Res_mapping_source` | Which CDR3 column produced the match (e.g. `"cdr3_col2"`) |

## Example workflow

```r
library(clonality)
library(Seurat)

# 1. Build the canonical paired clonality with tenx()
out <- tenx(
  data            = Cellranger9_BCR1,
  method          = "unique_paired",
  cell            = "B",
  only_productive = TRUE,
  only_true_cells = TRUE
)
# tenx() assigns Clonal.output.10xIGH, Clonal.output.10xIGK,
# Clonal.output.10xIGL, Clonal.output.10xIGH_IGK, Clonal.output.10xIGH_IGL
# to the global environment.

# 2. Embed into Seurat with sticky-ends rescue
seu <- add_seurat_metadata(
  seu         = my_seurat_object,
  clonal.list = ls(pattern = "^Clonal"),
  cell        = "B",
  sticky      = TRUE,
  verbose     = TRUE   # prints a per-chain purity summary after the rescue
)
```

`verbose = TRUE` prints something like:

```
Sticky-ends rescue summary
 chain singles matched median_purity rescued at_0.8 at_0.9 at_1
   IGH    4477    2801         1.000    1750    2264   1750 1576
   IGK    6134       0            NA       0       0      0    0   ← skipped (NA)
   IGL     504       0            NA       0       0      0    0   ← skipped (NA)
```

The `0`s for IGK / IGL come from the chain-skipped policy. To audit at any
later point:

```r
inspect_purity(seu)
```

## Visualizing the purity as an alluvial plot

The DNA-strand metaphor maps naturally onto an alluvial: left axis = the
single-chain cells, right axis = the paired clone they were rescued into,
ribbon thickness = number of cells, ribbon transparency = purity.

```r
library(ggplot2)
library(ggalluvial)
library(dplyr)

df <- seu[[]] %>%
  filter(!is.na(Res_dominant_purity)) %>%
  mutate(
    single_clone = .data$clonality,
    paired_clone = .data$Res_clonality
  )

ggplot(df,
       aes(axis1 = single_clone,
           axis2 = paired_clone,
           y = 1,
           alpha = Res_dominant_purity)) +
  geom_alluvium(aes(fill = Res_clonality.correction)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
  scale_x_discrete(limits = c("Single chain", "Paired clone"), expand = c(0.1, 0.1)) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_fill_manual(values = c("Normal" = "#999999",
                                "Corrected" = "#0c7c59",
                                "Not tested" = "#cccccc",
                                "Not tested (chain skipped)" = "#dddddd")) +
  theme_minimal()
```

A clean rescue produces wide, dark ribbons each going to one paired-clone
target. Promiscuous κ/λ singles (when you opt in) fan out across many targets
with low alpha — the visual that motivated the chain-aware defaults in the
first place.

## Migration from < v0.11

| Previous call | v0.11 equivalent |
|---|---|
| `add_seurat_metadata(seu, sticky = TRUE)` | Same; non-D chains now skipped by default. To rescue κ/λ pass `purity = list(IGH = 0.9, IGK = 0.95, IGL = 0.95)` |
| `add_seurat_metadata(seu, sticky = TRUE, purity = 0.8)` | Works with a back-compat message; consider `purity = NULL` for the chain-aware defaults |
| Want substring matching back | `match_mode = "substring"` |
| Want V-gene check off | `require_v_gene = FALSE` |
| Want every CDR3 length tried | `min_cdr3_aa = 0` |

## See also

- `?add_seurat_metadata` — full argument reference
- `?inspect_purity` — diagnostic helper
- `NEWS.md` — full v0.11 release notes
- The original purity analysis report
  (`05_Reports/sticky_ends_purity_analysis.md` in the ELvaron project)
  that motivated this overhaul.
