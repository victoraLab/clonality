#' A small IMGT TCR Beta sequences dataset
#'
#' Small IMGT TRB test dataset for `clonality()`.
#'
#' @format A data frame with 9 columns.
#' @details
#' Key columns: `Sequence_ID`, `V_GENE_and_allele`, `J_GENE_and_allele`, `JUNCTION`, `Sequence`.
#' Use `names(trb)` to see all columns.
"trb"

#' A small IMGT TCR Alpha sequences dataset
#'
#' Small IMGT TRA test dataset for `clonality()`.
#'
#' @format A data frame with 9 columns.
#' @details
#' Key columns: `Sequence_ID`, `V_GENE_and_allele`, `J_GENE_and_allele`, `JUNCTION`, `Sequence`.
#' Use `names(trb)` to see all columns.
"tra"

#' Chain Class Definitions for B and T Cells in 10x Data
#'
#' A list defining the acceptable chain combinations (classes) for B cells, T cells, and γδ T cells,
#' across two clonality methods: `"unique_paired"` and `"sticky_ends"`.
#'
#' These classes are used internally by the `assigntenx()` and `clonality()` functions to determine
#' how contigs should be grouped.
#'
#' @format A named list of length 2 (`unique_paired` and `sticky_ends`), each containing:
#' \describe{
#'   \item{BigCell.classes}{Character vector of BCR chain combinations    (e.g., `"IGH_IGK"`).}
#'   \item{TabCell.classes}{Character vector of αβ TCR chain combinations (e.g., `"TRA_TRB"`).}
#'   \item{TgdCell.classes}{Character vector of γδ TCR chain combinations (e.g., `"TRD_TRG"`).}
#' }
#'
#' @usage data(chain_classes)
#' @examples
#' data(chain_classes)
#' chain_classes$unique_paired$BCell.classes
"chain_classes"

#' IMGT BCR Summary File
#'
#' Processed IMGT HighV-QUEST summary output for BCR sequences.
#'
#' @format A data frame.
#' @details
#' Key columns: `Sequence_ID`, `VDOMAIN_Functionality`, `VGENE_and_allele`,
#' `JGENE_and_allele`, `DNA_Junction`.
#' Use `names(IMGT_BCR_Summary_File)` for all variables.
#'
#' @source Generated from IMGT `1_Summary.txt` and processed in `data-raw/...`.
"IMGT_BCR_Summary_File"

#' IMGT TCR Summary File
#'
#' A processed summary table of T cell receptor (TCR) sequences analyzed using the IMGT HighV-QUEST pipeline.
#' Each row represents a unique TCR sequence annotated with V(D)J gene usage, junction details, and
#' CDR/FR region lengths.
#'
#' @format A data frame.
#' @details
#' Key columns: `Sequence_ID`, `VDOMAIN_Functionality`, `VGENE_and_allele`,
#' `JGENE_and_allele`, `DNA_Junction`.
#' Use `names(IMGT_BCR_Summary_File)` for all variables.
#'
#' @source Generated from `1_Summary.txt`, the output of IMGT HighV-QUEST, and processed via `data-raw/IMGT_BCR_Summary_File.R`.
#'
#' @examples
#' head(IMGT_BCR_Summary_File)
"IMGT_TCR_Summary_File"

#' Cellranger3 TCR Contig Annotations
#'
#' Annotated TCR contigs from 10x Cell Ranger v3 (`filtered_contig_annotations.csv`-style).
#'
#' @format A tibble.
#' @details Use `names(Cellranger3_TCR)` for the full schema.
#' @examples
#' data("Cellranger3_TCR")
#' head(Cellranger3_TCR)
"Cellranger3_TCR"

#' Cellranger7 BCR Contig Annotations
#'
#' @name Cellranger7_BCR
#' @aliases Cellranger7_BCR1 Cellranger7_BCR2
#' @format A tibble with 31 columns.
#' @details Use `names(Cellranger7_BCR1)` for the full schema.
NULL

"Cellranger7_BCR1"
"Cellranger7_BCR2"

#' Cellranger9 BCR Contig Annotations
#'
#' Annotated B cell receptor (BCR) contigs from single-cell V(D)J sequencing experiments processed with
#' 10x Genomics Cell Ranger (v9+).
#'
#' @format Each dataset is a tibble with 31 columns:
#' @details Use `names(Cellranger9_BCR1)` for the full schema.
#' @source 10x Genomics Cell Ranger V(D)J pipeline, using `all_contig_annotations.csv` from samples SN082621 and SN062424.
#'
#' @examples
#' head(Cellranger9_BCR1)
"Cellranger9_BCR1"

#' Cellranger7 TCR γδ Contig Annotations
#'
#' Annotated T cell receptor (TCR) γδ contigs from single-cell V(D)J sequencing experiments
#' processed with 10x Genomics Cell Ranger (v7+).
#'
#' ⚠️ Note: Cell Ranger v7 does not fully support γδ TCR analysis. Consequently:
#' \itemize{
#'   \item The columns \code{is_cell} and \code{high_confidence} are unreliable and should not be used.
#'   \item Several annotation columns are not populated for γδ data and may be entirely \code{NA}.
#'   \item The \code{d_gene} column is sometimes informative and is retained.
#' }
#'
#' @format A tibble with N rows and 17+ columns (γδ-specific subset of the standard contig annotation schema):
#' @details Use `names(Cellranger7_TCRgd)` for the full schema.
#' @source Generated from 10x Genomics Cell Ranger v7 (`all_contig_annotations.csv`) for a γδ TCR dataset.
#'
#' @examples
#' head(Cellranger7_TCRgd)
"Cellranger7_TCRgd"
