test_that("add_seurat_metadata() errors if seu is NULL", {
  expect_error(add_seurat_metadata(NULL), "Please add Seurat object to merge clonality.")
})

test_that("add_seurat_metadata sticky maps singlets to dominant paired clonotype", {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    1,
    nrow = 2,
    ncol = 3,
    dimnames = list(c("gene1", "gene2"), c("pair1", "pair2", "single1"))
  )
  seu <- Seurat::CreateSeuratObject(counts = as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix"))

  pair_obj <- data.frame(
    barcodes = c("pair1", "pair2"),
    clonality = c("cloneA", "cloneA"),
    cdr3_IGH = c("IGHSEQ", "IGHSEQ"),
    stringsAsFactors = FALSE
  )

  single_obj <- data.frame(
    barcodes = "single1",
    clonality = "singleClone",
    cdr3_IGH = "IGHSEQ",
    stringsAsFactors = FALSE
  )

  assign("Clonal.output.10xIGH_IGK", pair_obj, envir = .GlobalEnv)
  assign("Clonal.output.10xIGH", single_obj, envir = .GlobalEnv)
  on.exit(rm(list = c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGH"), envir = .GlobalEnv), add = TRUE)

  res <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGH",
    cell = "B",
    sticky = TRUE
  )

  meta <- res[[]]
  expect_equal(meta["single1", "Res_clonality"], "cloneA")
  expect_equal(meta["single1", "Res_clonality.correction"], "Corrected")
  expect_equal(meta["single1", "Res_mapping_source"], "cdr3_IGH")
})

test_that("add_seurat_metadata purity blocks ambiguous sticky mappings", {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    1,
    nrow = 2,
    ncol = 3,
    dimnames = list(c("gene1", "gene2"), c("pair1", "pair2", "single1"))
  )
  seu <- Seurat::CreateSeuratObject(counts = as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix"))

  pair_obj <- data.frame(
    barcodes = c("pair1", "pair2"),
    clonality = c("cloneA", "cloneB"),
    cdr3_IGH = c("IGHSEQ", "IGHSEQ"),
    stringsAsFactors = FALSE
  )

  single_obj <- data.frame(
    barcodes = "single1",
    clonality = "singleClone",
    cdr3_IGH = "IGHSEQ",
    stringsAsFactors = FALSE
  )

  assign("Clonal.output.10xIGH_IGK", pair_obj, envir = .GlobalEnv)
  assign("Clonal.output.10xIGH", single_obj, envir = .GlobalEnv)
  on.exit(rm(list = c("Clonal.output.10xIGH_IGK", "Clonal.output.10xIGH"), envir = .GlobalEnv), add = TRUE)

  res <- add_seurat_metadata(
    seu = seu,
    clonal.list = c("Clonal.output.10xIGH", "Clonal.output.10xIGH_IGK"),
    cdr3_map = "cdr3_IGH",
    cell = "B",
    sticky = TRUE,
    purity = 0.75
  )

  meta <- res[[]]
  expect_equal(meta["single1", "Res_clonality"], "singleClone")
  expect_equal(meta["single1", "Res_clonality.correction"], "Not tested")
  expect_true(is.na(meta["single1", "Res_mapping_source"]))
})
