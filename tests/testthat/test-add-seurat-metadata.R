## ── Helpers ────────────────────────────────────────────────────────────────

make_seurat <- function(barcodes) {
  skip_if_not_installed("Seurat")
  counts <- matrix(
    1,
    nrow = 2,
    ncol = length(barcodes),
    dimnames = list(c("gene1", "gene2"), barcodes)
  )
  Seurat::CreateSeuratObject(
    counts = as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  )
}

cleanup_clonal_globals <- function(names) {
  rm(list = names, envir = .GlobalEnv)
}

# A 12-aa CDR3 long enough to clear the default min_cdr3_aa = 9 floor.
LONG_CDR3 <- "CARWGLRQEGFD"

## ── Original behavior tests (refreshed for v0.11 defaults) ─────────────────

test_that("add_seurat_metadata() errors if seu is NULL", {
  expect_error(add_seurat_metadata(NULL), "Please add Seurat object to merge clonality.")
})

test_that("sticky rescue maps an IGH single to the dominant paired clonotype", {
  seu <- make_seurat(c("pair1", "pair2", "single1"))

  pair_obj <- data.frame(
    barcodes = c("pair1", "pair2"),
    clonality = c("cloneA", "cloneA"),
    cdr3_IGH = c(LONG_CDR3, LONG_CDR3),
    stringsAsFactors = FALSE
  )

  single_obj <- data.frame(
    barcodes = "single1",
    clonality = "singleClone",
    cdr3_IGH = LONG_CDR3,
    stringsAsFactors = FALSE
  )

  assign("Clonal.output.10xIGH_IGK", pair_obj, envir = .GlobalEnv)
  assign("Clonal.output.10xIGH", single_obj, envir = .GlobalEnv)
  on.exit(cleanup_clonal_globals(c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGH")), add = TRUE)

  res <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGH",
    cell = "B",
    sticky = TRUE,
    require_v_gene = FALSE  # no V columns in this fixture
  )

  meta <- res[[]]
  expect_equal(meta["single1", "Res_clonality"], "cloneA")
  expect_equal(meta["single1", "Res_clonality.correction"], "Corrected")
  expect_equal(meta["single1", "Res_mapping_source"], "cdr3_IGH")
  expect_equal(meta["single1", "Res_dominant_purity"], 1)
})

test_that("ambiguous paired-pool maps below threshold leave the single uncorrected", {
  seu <- make_seurat(c("pair1", "pair2", "single1"))

  pair_obj <- data.frame(
    barcodes = c("pair1", "pair2"),
    clonality = c("cloneA", "cloneB"),
    cdr3_IGH = c(LONG_CDR3, LONG_CDR3),
    stringsAsFactors = FALSE
  )

  single_obj <- data.frame(
    barcodes = "single1",
    clonality = "singleClone",
    cdr3_IGH = LONG_CDR3,
    stringsAsFactors = FALSE
  )

  assign("Clonal.output.10xIGH_IGK", pair_obj, envir = .GlobalEnv)
  assign("Clonal.output.10xIGH", single_obj, envir = .GlobalEnv)
  on.exit(cleanup_clonal_globals(c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGH")), add = TRUE)

  expect_message(
    res <- add_seurat_metadata(
      seu = seu,
      clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
      cdr3_map = "cdr3_IGH",
      cell = "B",
      sticky = TRUE,
      purity = 0.75,         # scalar → triggers back-compat message
      require_v_gene = FALSE
    ),
    "scalar"
  )

  meta <- res[[]]
  expect_equal(meta["single1", "Res_clonality"], "singleClone")
  expect_equal(meta["single1", "Res_clonality.correction"], "Not tested")
  expect_equal(meta["single1", "Res_dominant_purity"], 0.5)
})

## ── v0.11 new behavior ─────────────────────────────────────────────────────

test_that("default per-chain purity skips IGK rescue (NA threshold)", {
  seu <- make_seurat(c("pair1", "pair2", "single1"))

  # IGK CDR3 shared across two unrelated paired clones (purity = 0.5)
  pair_obj <- data.frame(
    barcodes = c("pair1", "pair2"),
    clonality = c("cloneA", "cloneB"),
    cdr3_IGK = c(LONG_CDR3, LONG_CDR3),
    stringsAsFactors = FALSE
  )

  single_obj <- data.frame(
    barcodes = "single1",
    clonality = "singleClone",
    cdr3_IGK = LONG_CDR3,
    stringsAsFactors = FALSE
  )

  assign("Clonal.output.10xIGH_IGK", pair_obj, envir = .GlobalEnv)
  assign("Clonal.output.10xIGK", single_obj, envir = .GlobalEnv)
  on.exit(cleanup_clonal_globals(c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGK")), add = TRUE)

  res <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGK", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGK",
    cell = "B",
    sticky = TRUE,
    require_v_gene = FALSE
    # purity = NULL → default policy IGK = NA
  )

  meta <- res[[]]
  expect_equal(meta["single1", "Res_clonality.correction"], "Not tested (chain skipped)")
  expect_true(is.na(meta["single1", "Res_dominant_purity"]))
})

test_that("min_cdr3_aa floor blocks short CDR3s", {
  seu <- make_seurat(c("pair1", "pair2", "single1"))

  short <- "CARWGL"  # 6 aa, below default 9-aa floor

  pair_obj <- data.frame(
    barcodes = c("pair1", "pair2"),
    clonality = c("cloneA", "cloneA"),
    cdr3_IGH = c(short, short),
    stringsAsFactors = FALSE
  )
  single_obj <- data.frame(
    barcodes = "single1",
    clonality = "singleClone",
    cdr3_IGH = short,
    stringsAsFactors = FALSE
  )

  assign("Clonal.output.10xIGH_IGK", pair_obj, envir = .GlobalEnv)
  assign("Clonal.output.10xIGH", single_obj, envir = .GlobalEnv)
  on.exit(cleanup_clonal_globals(c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGH")), add = TRUE)

  res <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGH",
    cell = "B",
    sticky = TRUE,
    require_v_gene = FALSE
  )

  meta <- res[[]]
  expect_equal(meta["single1", "Res_clonality.correction"], "Not tested (cdr3 too short)")
  expect_true(is.na(meta["single1", "Res_dominant_purity"]))
})

test_that("match_mode = 'exact' rejects substring-only matches", {
  seu <- make_seurat(c("pair1", "single1"))

  pair_obj <- data.frame(
    barcodes = "pair1",
    clonality = "cloneA",
    cdr3_IGH = paste0("X", LONG_CDR3, "Y"),  # contains LONG_CDR3 as substring
    stringsAsFactors = FALSE
  )
  single_obj <- data.frame(
    barcodes = "single1",
    clonality = "singleClone",
    cdr3_IGH = LONG_CDR3,
    stringsAsFactors = FALSE
  )

  assign("Clonal.output.10xIGH_IGK", pair_obj, envir = .GlobalEnv)
  assign("Clonal.output.10xIGH", single_obj, envir = .GlobalEnv)
  on.exit(cleanup_clonal_globals(c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGH")), add = TRUE)

  # Exact (default): should NOT match
  res_exact <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGH",
    cell = "B",
    sticky = TRUE,
    require_v_gene = FALSE
  )
  expect_equal(res_exact[[]][["Res_clonality.correction"]][res_exact[[]][["Res_clonality"]] == "singleClone"][1],
               "Not tested (no match)")

  # Substring opt-in: should match
  res_sub <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGH",
    cell = "B",
    sticky = TRUE,
    match_mode = "substring",
    require_v_gene = FALSE
  )
  meta_sub <- res_sub[[]]
  expect_equal(meta_sub["single1", "Res_clonality"], "cloneA")
})

test_that("require_v_gene = TRUE filters paired matches by V gene", {
  seu <- make_seurat(c("pair1", "pair2", "single1"))

  pair_obj <- data.frame(
    barcodes  = c("pair1", "pair2"),
    clonality = c("cloneA", "cloneB"),
    cdr3_IGH  = c(LONG_CDR3, LONG_CDR3),
    v_genes_IGH = c("IGHV1-1", "IGHV2-2"),
    stringsAsFactors = FALSE
  )
  single_obj <- data.frame(
    barcodes  = "single1",
    clonality = "singleClone",
    cdr3_IGH  = LONG_CDR3,
    v_genes_IGH = "IGHV1-1",
    stringsAsFactors = FALSE
  )

  assign("Clonal.output.10xIGH_IGK", pair_obj, envir = .GlobalEnv)
  assign("Clonal.output.10xIGH", single_obj, envir = .GlobalEnv)
  on.exit(cleanup_clonal_globals(c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGH")), add = TRUE)

  # Without V check: ambiguous (purity 0.5), denied at default 0.9
  res_no_v <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGH",
    cell = "B",
    sticky = TRUE,
    require_v_gene = FALSE
  )
  expect_equal(res_no_v[[]]["single1", "Res_dominant_purity"], 0.5)

  # With V check: only pair1 matches → purity 1.0, rescued
  res_v <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGH",
    cell = "B",
    sticky = TRUE,
    require_v_gene = TRUE
  )
  expect_equal(res_v[[]]["single1", "Res_clonality"], "cloneA")
  expect_equal(res_v[[]]["single1", "Res_dominant_purity"], 1)
})

test_that("Res_dominant_purity column is always present after sticky", {
  seu <- make_seurat(c("pair1", "single1"))

  pair_obj <- data.frame(
    barcodes = "pair1", clonality = "cloneA",
    cdr3_IGH = LONG_CDR3, stringsAsFactors = FALSE
  )
  single_obj <- data.frame(
    barcodes = "single1", clonality = "singleClone",
    cdr3_IGH = LONG_CDR3, stringsAsFactors = FALSE
  )

  assign("Clonal.output.10xIGH_IGK", pair_obj, envir = .GlobalEnv)
  assign("Clonal.output.10xIGH", single_obj, envir = .GlobalEnv)
  on.exit(cleanup_clonal_globals(c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGH")), add = TRUE)

  res <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGH",
    cell = "B",
    sticky = TRUE,
    require_v_gene = FALSE
  )

  expect_true("Res_dominant_purity" %in% colnames(res[[]]))
  expect_true(is.numeric(res[[]][["Res_dominant_purity"]]))
})

## ── normalize_airr_to_cellranger ──────────────────────────────────────────

test_that("normalize_airr_to_cellranger collapses comma-separated alleles", {
  airr <- data.frame(
    v_call = c("IGHV1-1*01,IGHV1-1*02", "IGHV2-2*01"),
    j_call = c("IGHJ1*01", "IGHJ2*01,IGHJ2*02"),
    stringsAsFactors = FALSE
  )
  out <- normalize_airr_to_cellranger(airr)
  expect_equal(out$v_call, c("IGHV1-1*01", "IGHV2-2*01"))
  expect_equal(out$j_call, c("IGHJ1*01", "IGHJ2*01"))
})

test_that("normalize_airr_to_cellranger fills c_gene from external column", {
  airr <- data.frame(
    c_call = c("", NA_character_, "IGHM"),
    isotype = c("IGHA", "IGHG1", "IGHE"),
    stringsAsFactors = FALSE
  )
  out <- normalize_airr_to_cellranger(airr, c_gene_col = "isotype")
  expect_equal(out$c_call, c("IGHA", "IGHG1", "IGHM"))
})

test_that("normalize_airr_to_cellranger coerces productive to lowercase strings", {
  airr <- data.frame(
    productive = c(TRUE, FALSE, NA),
    stringsAsFactors = FALSE
  )
  out <- normalize_airr_to_cellranger(airr)
  expect_equal(out$productive, c("true", "false", NA_character_))

  airr2 <- data.frame(productive = c("T", "F", "yes", "weird"), stringsAsFactors = FALSE)
  out2 <- normalize_airr_to_cellranger(airr2)
  expect_equal(out2$productive, c("true", "false", "true", NA_character_))
})

## ── inspect_purity ────────────────────────────────────────────────────────

test_that("inspect_purity returns a per-chain summary from a data frame", {
  meta <- data.frame(
    Res_clonal_source = c(rep("Clonal.output.10xIGH", 4), rep("Clonal.output.10xIGK", 4)),
    Res_dominant_purity = c(1, 1, 0.5, NA, 0.4, 0.3, 0.6, 1),
    Res_clonality.correction = c("Corrected", "Corrected", "Not tested", "Not tested",
                                  "Not tested", "Not tested", "Not tested", "Corrected"),
    stringsAsFactors = FALSE
  )

  out <- suppressMessages(inspect_purity(meta))
  expect_true(is.data.frame(out))
  expect_setequal(out$chain, c("IGH", "IGK"))
  expect_equal(out$singles[out$chain == "IGH"], 4)
  expect_equal(out$matched[out$chain == "IGH"], 3)
  expect_equal(out$rescued[out$chain == "IGH"], 2)
})
