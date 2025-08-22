test_that("cellranger_version_parser can detect the version selected is wrong", {
  # These objects must be available in the test environment
  expect_true(exists("Cellranger3_TCR"))
  expect_true(exists("Cellranger7_BCR1"))

  # Invalid version check
  expect_error(cellranger_version_parser(Cellranger3_TCR, cellranger_version = 7), "expected to have 31 columns")
  expect_error(cellranger_version_parser(Cellranger7_BCR1, cellranger_version = 3), "expected to have 18 columns")
})

test_that("tenx rejects invalid clonality method", {
  expect_error(
    tenx(data = NULL, method = "not_a_method")
  )
})

test_that("tenx drops empty/NA CDR3s before grouping", {
  df <- data.frame(
    barcode = c("A", "A", "B"),
    chain   = c("IGH", "IGK", "IGH"),
    cdr3_nt = c("", NA, "ATG"),          # first two should be removed
    productive = c("True", "True", "True"),
    is_cell    = c("True", "True", "True"),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings(
    cellranger_version_parser = function(data, cellranger_version) data,  # <- named
    assigntenx = function(list.pairs, ...) list.pairs,                     # <- named
    .env = environment(tenx)
    , {
      res <- tenx(df)

      # Only barcode B remains (A's contigs dropped for empty/NA CDR3)
      expect_named(res, "IGH")               # one chain class
      expect_equal(length(res[["IGH"]]), 1L) # one cell (barcode B)
    })
})


test_that("only_productive=TRUE filters out non-productive contigs", {
  df <- data.frame(
    barcode   = c("A","A","B","B"),
    chain     = c("IGH","IGK","IGH","IGK"),
    cdr3_nt   = c("ATG","TAA","CCC","GGG"),
    productive= c("True","None","True","None"),   # IGK of both A and B are non-productive
    is_cell   = c("True","True","True","True"),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings(
    cellranger_version_parser = function(data, cellranger_version) data,  # <- named
    assigntenx = function(list.pairs, ...) list.pairs,                     # <- named
    .env = environment(tenx)
    , {
    res <- tenx(df, only_productive = TRUE)
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
    productive= c("True","True","True","True"),
    is_cell   = c("True","True","False","False"),  # drop barcode B completely
    stringsAsFactors = FALSE
  )

  with_mocked_bindings(
    cellranger_version_parser = function(data, cellranger_version) data,  # <- named
    assigntenx = function(list.pairs, ...) list.pairs,                     # <- named
    .env = environment(tenx)
    , {
    res <- tenx(df, only_true_cells = TRUE)
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
    productive= c("True","True","True","True","True"),
    is_cell   = c("True","True","True","True","True"),
    stringsAsFactors = FALSE
  )

  with_mocked_bindings(
    cellranger_version_parser = function(data, cellranger_version) data,  # <- named
    assigntenx = function(list.pairs, ...) list.pairs,                     # <- named
    .env = environment(tenx)
    , {
    res <- tenx(df, cell = "Tgd")
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

