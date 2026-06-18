test_that("cellranger_version_parser can detect the version selected is wrong", {
  # These objects must be available in the test environment
  expect_true(exists("Cellranger3_TCR"))
  expect_true(exists("Cellranger7_BCR1"))

  # Invalid version check
  expect_error(cellranger_version_parser(Cellranger3_TCR, cellranger_version = 7), "expected to have 31 columns")
  expect_error(cellranger_version_parser(Cellranger7_BCR1, cellranger_version = 3), "expected to have 18 columns")
})

test_that("cellranger_version_parser can infer version from bundled data", {
  parsed <- cellranger_version_parser(Cellranger7_BCR1, cellranger_version = NULL)

  expect_s3_class(parsed, "data.frame")
  expect_equal(ncol(parsed), 31L)
})

test_that("tenx rejects invalid clonality method", {
  expect_error(
    tenx(data = NULL, method = "not_a_method")
  )
})

test_that("tenx works when cellranger_version is omitted", {
  testthat::local_mocked_bindings(
    assigntenx = function(list.pairs, ...) list.pairs,
    .env = environment(tenx)
  )

  res <- tenx(
    data = Cellranger7_BCR1,
    only_productive = TRUE,
    only_true_cells = TRUE
  )

  expect_true(length(res) > 0)
})

test_that("tenx drops empty/NA CDR3s before grouping", {
  df <- data.frame(
    barcode = c("A", "A", "B"),
    chain   = c("IGH", "IGK", "IGH"),
    cdr3_nt = c("", NA, "ATG"),
    v_gene  = c("IGHV1", "IGKV1", "IGHV1"),
    j_gene  = c("IGHJ1", "IGKJ1", "IGHJ1"),
    productive = c("True", "True", "True"),
    is_cell    = c("True", "True", "True"),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    assigntenx = function(list.pairs, ...) list.pairs,
    .env = environment(tenx)
  )

  res <- tenx(df, ignore_version = TRUE)

  expect_named(res, "IGH")
  expect_equal(length(res[["IGH"]]), 1L)
})

test_that("only_productive=TRUE filters out non-productive contigs", {
  df <- data.frame(
    barcode   = c("A","A","B","B"),
    chain     = c("IGH","IGK","IGH","IGK"),
    cdr3_nt   = c("ATG","TAA","CCC","GGG"),
    v_gene    = c("IGHV1","IGKV1","IGHV1","IGKV1"),
    j_gene    = c("IGHJ1","IGKJ1","IGHJ1","IGKJ1"),
    productive= c("True","None","True","None"),   # IGK of both A and B are non-productive
    is_cell   = c("True","True","True","True"),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings(
    assigntenx = function(list.pairs, ...) list.pairs,                     # <- named
    .env = environment(tenx)
    , {
    res <- tenx(df, only_productive = TRUE, ignore_version = TRUE)
    # All remaining contigs are IGH; classes should be just "IGH" (no IGK)
    expect_named(res, "IGH")
    # two barcodes (A,B) each with a single IGH entry
    expect_equal(length(res[["IGH"]]), 2L)
  })
})

test_that("only_true_cells=TRUE filters by is_cell", {
  df <- data.frame(
    barcode   = c("A","A","B","B"),
    chain     = c("IGH","IGK","IGH","IGK"),
    cdr3_nt   = c("ATG","TAA","CCC","GGG"),
    v_gene    = c("IGHV1","IGKV1","IGHV1","IGKV1"),
    j_gene    = c("IGHJ1","IGKJ1","IGHJ1","IGKJ1"),
    productive= c("True","True","True","True"),
    is_cell   = c("True","True","False","False"),  # drop barcode B completely
    stringsAsFactors = FALSE
  )

  with_mocked_bindings(
    assigntenx = function(list.pairs, ...) list.pairs,                     # <- named
    .env = environment(tenx)
    , {
    res <- tenx(df, only_true_cells = TRUE, ignore_version = TRUE)
    # Only barcode A remains; its chain class is IGH_IGK (ordered)
    expect_true("IGH_IGK" %in% names(res))
    expect_equal(length(res[["IGH_IGK"]]), 1L)
  })
})

test_that("Tgd special handling removes chain=='None' and CDR3 with '*'", {
  df <- data.frame(
    barcode   = c("A","A","B","B","C"),
    chain     = c("TRG","None","TRD","TRG","TRD"),
    cdr3      = c("CAXYZ", "CAXYZ", "CA*BAD", "CAGGG", "CATT"),  # one with '*'
    cdr3_nt   = c("AAA","TTT","CCC","GGG","CCC"),
    v_gene    = c("TRGV1","TRGV1","TRDV1","TRGV1","TRDV1"),
    j_gene    = c("TRGJ1","TRGJ1","TRDJ1","TRGJ1","TRDJ1"),
    productive= c("True","True","True","True","True"),
    is_cell   = c("True","True","True","True","True"),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings(
    assigntenx = function(list.pairs, ...) list.pairs,                     # <- named
    .env = environment(tenx)
    , {
    res <- tenx(df, cell = "Tgd", ignore_version = TRUE)
    # Drops barcode A's 'None' contig, and drops barcode B's TRD with '*'
    # Remaining classes:
    #   A: TRG
    #   B: TRG   (TRD dropped)
    #   C: TRD
    expect_setequal(names(res), c("TRD","TRG"))
    expect_equal(length(res[["TRG"]]), 2L)  # barcodes A and B
    expect_equal(length(res[["TRD"]]), 1L)  # barcode C
  })
})

## ── multi-sample (project_names) ─────────────────────────────────────────────

make_bcr_df <- function(barcode_prefix = "") {
  data.frame(
    barcode    = paste0(barcode_prefix, c("A", "A", "B", "B")),
    chain      = c("IGH", "IGK", "IGH", "IGK"),
    cdr3_nt    = c("ATGATG", "TAATAA", "CCCGGG", "GGGTTT"),
    v_gene     = c("IGHV1", "IGKV1", "IGHV1", "IGKV1"),
    j_gene     = c("IGHJ1", "IGKJ1", "IGHJ1", "IGKJ1"),
    productive = c("True", "True", "True", "True"),
    is_cell    = c("True", "True", "True", "True"),
    stringsAsFactors = FALSE
  )
}

test_that("multi-sample: barcodes are suffixed with project_names", {
  df1 <- make_bcr_df()
  df2 <- make_bcr_df()

  with_mocked_bindings(
    assigntenx = function(list.pairs, ...) list.pairs,
    .env = environment(tenx),
    {
      res <- tenx(
        data          = list(df1, df2),
        project_names = c("S1", "S2"),
        ignore_version = TRUE
      )

      all_barcodes <- unlist(lapply(res, function(cls) {
        vapply(cls, function(cell_df) cell_df$barcode[[1]], character(1))
      }))

      expect_true(all(grepl("_S1$|_S2$", all_barcodes)))
      expect_true(any(grepl("_S1$", all_barcodes)))
      expect_true(any(grepl("_S2$", all_barcodes)))
    }
  )
})

test_that("multi-sample: barcodes from different samples stay unique", {
  # Both samples have identical raw barcodes; after suffix they must differ.
  df1 <- make_bcr_df()
  df2 <- make_bcr_df()

  with_mocked_bindings(
    assigntenx = function(list.pairs, ...) list.pairs,
    .env = environment(tenx),
    {
      res <- tenx(
        data          = list(df1, df2),
        project_names = c("S1", "S2"),
        ignore_version = TRUE
      )

      all_barcodes <- unlist(lapply(res, function(cls) {
        vapply(cls, function(cell_df) cell_df$barcode[[1]], character(1))
      }))

      # No duplicated barcodes in the combined result
      expect_false(any(duplicated(all_barcodes)))
    }
  )
})

test_that("multi-sample: error when project_names is missing", {
  df1 <- make_bcr_df()
  df2 <- make_bcr_df()

  expect_error(
    tenx(data = list(df1, df2), ignore_version = TRUE),
    "`project_names` must be provided"
  )
})

test_that("multi-sample: error when project_names length mismatches data", {
  df1 <- make_bcr_df()
  df2 <- make_bcr_df()

  expect_error(
    tenx(data = list(df1, df2), project_names = "S1", ignore_version = TRUE),
    "same length"
  )
})

test_that("multi-sample: error when project_names are not unique", {
  df1 <- make_bcr_df()
  df2 <- make_bcr_df()

  expect_error(
    tenx(data = list(df1, df2), project_names = c("S1", "S1"), ignore_version = TRUE),
    "unique"
  )
})

test_that("tenx normalizes common AIRR-style compiled column names", {
  df <- data.frame(
    cell_id = c("A", "A", "B", "B"),
    locus = c("IGH", "IGK", "IGH", "IGK"),
    junction = c("ATG", "TAA", "CCC", "GGG"),
    junction_aa = c("CAR", "CTR", "CQQ", "CGG"),
    v_call_genotyped = c("IGHV1", "IGKV1", "IGHV3", "IGKV3"),
    j_call = c("IGHJ1", "IGKJ1", "IGHJ3", "IGKJ3"),
    c_call = c("IGHM", "IGKC", "IGHG", "IGKC"),
    clone_id = c("cl1", "cl1", "cl2", "cl2"),
    productive = c(TRUE, TRUE, TRUE, TRUE),
    is_cell = c(TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings(
    assigntenx = function(list.pairs, ...) list.pairs,
    .env = environment(tenx)
    , {
      res <- tenx(df, ignore_version = TRUE)
      expect_true("IGH_IGK" %in% names(res))
      expect_equal(length(res[["IGH_IGK"]]), 2L)
    })
})
