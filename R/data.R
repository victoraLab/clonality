#' A small IMGT TCR Beta sequences dataset.
#'
#' A dataset containing TCR Beta sequences to use as a test for clonality.
#'
#' @format A data frame with 49 rows and 9 variables:
#' \describe{
#'   \item{Sequence_number}{Sequence number, not needed}
#'   \item{Sequence_ID}{Unique sequence ID, needed}
#'   \item{V_DOMAIN_Functionality}{If the V sequence are functional, not needed}
#'   \item{V_GENE_and_allele}{V gene name, needed}
#'   \item{J_GENE_and_allele}{J gene name, needed}
#'   \item{D_GENE_and_allele}{D gene name, not needed}
#'   \item{JUNCTION_frame}{If CDR3 junction is in-frame, not needed}
#'   \item{JUNCTION}{Junction sequence}
#'   \item{sequence}{Full TCR sequence}
#'   ...
#' }
#'
"trb"

#' A small IMGT TCR Alpha sequences dataset.
#'
#' A dataset containing TCR Alpha sequences to use as a test for clonality.
#'
#' @format A data frame with 49 rows and 9 variables:
#' \describe{
#'   \item{Sequence_number}{Sequence number, not needed}
#'   \item{Sequence_ID}{Unique sequence ID, needed}
#'   \item{V_DOMAIN_Functionality}{If the V sequence are functional, not needed}
#'   \item{V_GENE_and_allele}{V gene name, needed}
#'   \item{J_GENE_and_allele}{J gene name, needed}
#'   \item{D_GENE_and_allele}{D gene name, not needed}
#'   \item{JUNCTION_frame}{If CDR3 junction is in-frame, not needed}
#'   \item{JUNCTION}{Junction sequence}
#'   \item{sequence}{Full TCR sequence}
#'   ...
#' }
#'
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
#' A processed summary table of B cell receptor (BCR) sequences analyzed using the IMGT HighV-QUEST pipeline.
#' Each row represents a unique immunoglobulin sequence annotated with V(D)J gene usage, junction details, and
#' CDR/FR region lengths.
#'
#' @format A data frame with N rows and 34 columns:
#' \describe{
#'   \item{Sequence_ID}{Unique identifier for the sequence assigned by IMGT.}
#'   \item{VDOMAIN_Functionality}{Functionality status of the V domain (e.g., productive, unproductive).}
#'   \item{VGENE_and_allele}{Identified V gene and allele.}
#'   \item{JGENE_and_allele}{Identified J gene and allele.}
#'   \item{DNA_Junction}{Nucleotide sequence of the junction region.}
#' }
#'
#' Columns not listed here include detailed scores, region insertions/deletions, trimming counts, comments, and
#' alternative identity metrics, all retained for reference and advanced filtering.
#'
#' @source Generated from `1_Summary.txt`, the output of IMGT HighV-QUEST, and processed via `data-raw/IMGT_BCR_Summary_File.R`.
#'
#' @examples
#' head(IMGT_BCR_Summary_File)
"IMGT_BCR_Summary_File"

#' IMGT TCR Summary File
#'
#' A processed summary table of T cell receptor (TCR) sequences analyzed using the IMGT HighV-QUEST pipeline.
#' Each row represents a unique TCR sequence annotated with V(D)J gene usage, junction details, and
#' CDR/FR region lengths.
#'
#' @format A data frame with N rows and 34 columns:
#' \describe{
#'   \item{Sequence_ID}{Unique identifier for the sequence assigned by IMGT.}
#'   \item{VDOMAIN_Functionality}{Functionality status of the V domain (e.g., productive, unproductive).}
#'   \item{VGENE_and_allele}{Identified V gene and allele.}
#'   \item{JGENE_and_allele}{Identified J gene and allele.}
#'   \item{DNA_Junction}{Nucleotide sequence of the junction region.}
#' }
#'
#' Columns not listed here include detailed scores, region insertions/deletions, trimming counts, comments, and
#' alternative identity metrics, all retained for reference and advanced filtering.
#'
#' @source Generated from `1_Summary.txt`, the output of IMGT HighV-QUEST, and processed via `data-raw/IMGT_BCR_Summary_File.R`.
#'
#' @examples
#' head(IMGT_BCR_Summary_File)
"IMGT_TCR_Summary_File"

#' Cellranger3 TCR Contig Annotations
#'
#' A tibble of annotated T cell receptor (TCR) contigs from single-cell V(D)J sequencing,
#' formatted as output from 10x Genomics Cell Ranger v3- (`filtered_contig_annotations.csv`).
#' Each row corresponds to a single TCR contig (productive or non-productive) from a single cell barcode.
#'
#' @format A tibble with N rows and 18 columns.
#' \describe{
#'   \item{barcode}{Cell barcode identifying the source cell.}
#'   \item{is_cell}{Logical indicating whether the barcode passed cell-calling filters.}
#'   \item{contig_id}{Unique identifier for the contig (cell barcode + contig index).}
#'   \item{high_confidence}{Logical indicating if the contig passed quality filters.}
#'   \item{length}{Length of the contig sequence in base pairs.}
#'   \item{chain}{TCR chain type (e.g., TRA, TRB).}
#'   \item{v_gene}{Assigned V gene segment (e.g., TRAV12-3).}
#'   \item{d_gene}{Assigned D gene segment (only for TRB chains; may be None).}
#'   \item{j_gene}{Assigned J gene segment.}
#'   \item{c_gene}{Assigned C gene segment (constant region).}
#'   \item{full_length}{Logical indicating if the contig represents a full-length transcript.}
#'   \item{productive}{String label indicating if the contig is productive ("True"/"None").}
#'   \item{cdr3}{Amino acid sequence of the CDR3 region, if available.}
#'   \item{cdr3_nt}{Nucleotide sequence of the CDR3 region, if available.}
#'   \item{reads}{Number of reads supporting the contig.}
#'   \item{umis}{Number of UMIs supporting the contig.}
#'   \item{raw_clonotype_id}{Clonotype ID assigned to the contig (groups related TCRs from the same cell).}
#'   \item{raw_consensus_id}{Consensus ID for the contig (or `None` if not called).}
#' }
#'
#' @source Generated from 10x Genomics Cell Ranger v3 filtered contig annotations file.
#'
#' @examples
#' head(Cellranger3_TCR)
"Cellranger3_TCR"

#' Cellranger7 BCR Contig Annotations
#'
#' Annotated B cell receptor (BCR) contigs from single-cell V(D)J sequencing experiments processed with
#' 10x Genomics Cell Ranger (v7+). These data frames contain one row per contig and include V(D)J gene assignments,
#' CDR and FWR regions, read support, and clonotype IDs.
#'
#' The datasets are:
#' \itemize{
#'   \item \strong{Cellranger7_BCR1}: Derived from sample \code{SN082621}.
#'   \item \strong{Cellranger7_BCR2}: Derived from sample \code{SN062424}.
#' }
#'
#' Both datasets have the same structure and column names.
#'
#' @format Each dataset is a tibble with 31 columns:
#' \describe{
#'   \item{barcode}{10x cell barcode identifying the source cell.}
#'   \item{is_cell}{Logical indicating if the barcode was called as a valid cell.}
#'   \item{contig_id}{Unique identifier for the contig (cell barcode + contig index).}
#'   \item{high_confidence}{Logical indicating high-confidence contig status.}
#'   \item{length}{Length of the contig sequence (nt).}
#'   \item{chain}{Immunoglobulin chain (e.g., IGH, IGK, IGL).}
#'   \item{v_gene, d_gene, j_gene, c_gene}{Assigned V, D, J, and constant region genes.}
#'   \item{fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4}{Amino acid sequences of framework and CDR regions.}
#'   \item{fwr1_nt, ..., fwr4_nt, cdr1_nt, ..., cdr3_nt}{Nucleotide sequences of each region.}
#'   \item{productive}{Logical indicating if the contig encodes a productive rearrangement.}
#'   \item{full_length}{Logical indicating full-length contig.}
#'   \item{reads}{Number of reads supporting the contig.}
#'   \item{umis}{Number of UMIs supporting the contig.}
#'   \item{raw_clonotype_id}{Assigned clonotype identifier (raw output from Cell Ranger).}
#'   \item{raw_consensus_id}{Consensus sequence ID if available.}
#'   \item{exact_subclonotype_id}{Exact subclonotype cluster ID.}
#' }
#'
#' @source 10x Genomics Cell Ranger V(D)J pipeline, using `all_contig_annotations.csv` from samples SN082621 and SN062424.
#'
#' @examples
#' head(Cellranger7_BCR1)
#' head(Cellranger7_BCR2)
"Cellranger7_BCR1"
"Cellranger7_BCR2"

#' Cellranger9 BCR Contig Annotations
#'
#' Annotated B cell receptor (BCR) contigs from single-cell V(D)J sequencing experiments processed with
#' 10x Genomics Cell Ranger (v9+).
#'
#' @format Each dataset is a tibble with 31 columns:
#' \describe{
#'   \item{sample}{Cellranger vdj running sample ID}
#'   \item{barcode}{10x cell barcode identifying the source cell.}
#'   \item{is_cell}{Logical indicating if the barcode was called as a valid cell.}
#'   \item{contig_id}{Unique identifier for the contig (cell barcode + contig index).}
#'   \item{high_confidence}{Logical indicating high-confidence contig status.}
#'   \item{length}{Length of the contig sequence (nt).}
#'   \item{chain}{Immunoglobulin chain (e.g., IGH, IGK, IGL).}
#'   \item{v_gene, d_gene, j_gene, c_gene}{Assigned V, D, J, and constant region genes.}
#'   \item{fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4}{Amino acid sequences of framework and CDR regions.}
#'   \item{fwr1_nt, ..., fwr4_nt, cdr1_nt, ..., cdr3_nt}{Nucleotide sequences of each region.}
#'   \item{productive}{Logical indicating if the contig encodes a productive rearrangement.}
#'   \item{full_length}{Logical indicating full-length contig.}
#'   \item{reads}{Number of reads supporting the contig.}
#'   \item{umis}{Number of UMIs supporting the contig.}
#'   \item{raw_clonotype_id}{Assigned clonotype identifier (raw output from Cell Ranger).}
#'   \item{raw_consensus_id}{Consensus sequence ID if available.}
#'   \item{exact_subclonotype_id}{Exact subclonotype cluster ID.}
#' }
#'
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
#' \describe{
#'   \item{barcode}{10x cell barcode identifying the source cell. Use with caution, as cell-calling may be incomplete.}
#'   \item{contig_id}{Unique identifier for the contig (cell barcode + contig index).}
#'   \item{length}{Length of the contig sequence (nt).}
#'   \item{chain}{TCR chain type (TRG or TRD).}
#'   \item{v_gene}{Assigned V gene segment.}
#'   \item{d_gene}{Assigned D gene segment (only relevant for δ chains; may be NA).}
#'   \item{j_gene}{Assigned J gene segment.}
#'   \item{c_gene}{Assigned constant region gene.}
#'   \item{full_length}{Logical indicating if the contig is full-length (provided by Cell Ranger).}
#'   \item{productive}{Logical indicating whether the contig encodes a productive rearrangement.}
#'   \item{cdr3}{Amino acid sequence of the CDR3 region (if detected).}
#'   \item{cdr3_nt}{Nucleotide sequence of the CDR3 region.}
#'   \item{reads}{Number of reads supporting the contig.}
#'   \item{umis}{Number of UMIs supporting the contig.}
#'   \item{raw_clonotype_id}{Clonotype ID assigned by Cell Ranger.}
#'   \item{raw_consensus_id}{Consensus sequence ID if available.}
#'   \item{exact_subclonotype_id}{Exact subclonotype identifier.}
#' }
#'
#' @source Generated from 10x Genomics Cell Ranger v7 (`all_contig_annotations.csv`) for a γδ TCR dataset.
#'
#' @examples
#' head(Cellranger7_TCRgd)
"Cellranger7_TCRgd"
