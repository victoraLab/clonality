test_that("clonality works on default test dataset", {
  res <- clonality(
    data = IMGT_BCR_Summary_File,
    ident_col = "Sequence_ID",
    vgene_col = "VGENE_and_allele",
    jgene_col = "JGENE_and_allele",
    cdr3_col = "DNA_Junction",
    cell = "B",
    output_original = FALSE,
    mismatch = 100
  )

  expect_s3_class(res, "data.frame")
  expect_true("clonality" %in% names(res))
  expect_false(any(is.na(res$clonality)))
})

test_that("clonality removes empty and NA CDR3s", {
  df <- data.frame(
    Sequence_ID = c("X", "X1", "X2"),
    VGENE_and_allele = c("TRAV1*01 F", "TRAV1*01 F", "TRAV1*01 F"),
    JGENE_and_allele = c("TRAJ1*01 F", "TRAJ1*01 F", "TRAJ1*01 F"),
    JUNCTION = c("", "CXXXXXXRC", NA),
    stringsAsFactors = FALSE
  )

  result <- clonality(df)
  # expected: (rows without empty JUNCTION) - (total rows)
  expected <- nrow(df) - sum(is.na(df$JUNCTION) | df$JUNCTION %in% "")
  expect_equal(nrow(result), expected)
})

test_that("simplify_genenames simplifies TRAV/TRAJ correctly", {
  input <- data.frame(
    v_genes = "TRAV12-2*01 F",
    j_genes = "TRAJ33*01 F",
    stringsAsFactors = FALSE
  )
  result <- simplify_genenames(input, cell = "T")
  expect_equal(result$clonal.df$v_genes, "TRAV12-2")
  expect_equal(result$clonal.df$j_genes, "TRAJ33")
})

test_that("original columns and [vj]_all columns are generated", {
  input <- data.frame(
    v_genes = "TRAV12-2*01 F, or TRAV13-1*01",
    j_genes = "TRAJ33*01 F",
    stringsAsFactors = FALSE
  )
  result <- simplify_genenames(input, "T")$clonal.df
  expect_true("v_genes_original" %in% names(result))
  expect_true(grepl("TRAV12-2;TRAV13-1", result$v_genes_all))
})


test_that("simplify_genenames throws for invalid cell type", {
  expect_error(simplify_genenames(data.frame(v_genes = "", j_genes = ""), "X"))

  })

test_that("clonality does not write an xlsx file unless file is requested", {
  skip_if_not_installed("openxlsx")

  input_path <- tempfile(fileext = ".xlsx")
  auto_output_path <- paste0(input_path, ".xlsx")

  openxlsx::write.xlsx(IMGT_BCR_Summary_File[1:5, ], input_path, rowNames = FALSE)

  res <- clonality(
    data = input_path,
    ident_col = "Sequence_ID",
    vgene_col = "VGENE_and_allele",
    jgene_col = "JGENE_and_allele",
    cdr3_col = "DNA_Junction",
    cell = "B",
    output_original = FALSE,
    mismatch = 100
  )

  expect_s3_class(res, "data.frame")
  expect_false(file.exists(auto_output_path))
})

test_that("mismatch controls normalized CDR3 distance tolerance, not similarity", {
  df <- data.frame(
    Sequence_ID = c("s1", "s2"),
    VGENE_and_allele = c("IGHV1-1*01", "IGHV1-1*01"),
    JGENE_and_allele = c("IGHJ1*01", "IGHJ1*01"),
    JUNCTION = c("AAAA", "AAAT"),
    stringsAsFactors = FALSE
  )

  strict <- clonality(df, cell = "B", mismatch = 0)
  tolerant <- clonality(df, cell = "B", mismatch = 25)

  expect_equal(length(unique(strict$clonality)), 2L)
  expect_equal(length(unique(tolerant$clonality)), 1L)
})

## ── v0.12: single CDR3 column with vector mismatch (stringency scan) ──────

test_that("single cdr3_col + vector mismatch yields one clonality_mm<v> column per threshold", {
  df <- data.frame(
    Sequence_ID = c("s1", "s2", "s3"),
    VGENE_and_allele = c("IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*01"),
    JGENE_and_allele = c("IGHJ1*01", "IGHJ1*01", "IGHJ1*01"),
    JUNCTION = c("CARWGLRQEGFD", "CARWGLRQEGFE", "CARWGLRQEGFF"),
    stringsAsFactors = FALSE
  )

  res <- clonality(df, cell = "B",
                   cdr3_col = "JUNCTION",
                   mismatch = c(0, 100))

  # Tag auto-fallback to "mm<value>" when the same CDR3 col is reused
  expect_true("clonality_mm0" %in% names(res))
  expect_true("clonality_mm100" %in% names(res))

  # mm0: each unique CDR3 distinct
  expect_equal(length(unique(res$clonality_mm0)), 3L)
  # mm100: all collapse into one
  expect_equal(length(unique(res$clonality_mm100)), 1L)
})

test_that("recycling is symmetric: scalar cdr3_col + vector mismatch and vice versa", {
  df <- data.frame(
    Sequence_ID = c("s1", "s2"),
    VGENE_and_allele = c("IGHV1-1*01", "IGHV1-1*01"),
    JGENE_and_allele = c("IGHJ1*01", "IGHJ1*01"),
    JUNCTION = c("CARWGLRQEGFD", "CARWGLRQEGFE"),
    stringsAsFactors = FALSE
  )

  # 1 col, 2 thresholds → K=2
  r1 <- clonality(df, cell = "B", cdr3_col = "JUNCTION", mismatch = c(0, 100))
  expect_equal(sum(grepl("^clonality", names(r1))), 2L)

  # 2 cols, 1 threshold → K=2 (mismatch recycled)
  r2 <- clonality(df, cell = "B", cdr3_col = c("JUNCTION", "JUNCTION"), mismatch = 0)
  expect_equal(sum(grepl("^clonality", names(r2))), 2L)
})

test_that("mismatched vector lengths (>1 vs >1, unequal) error clearly", {
  df <- data.frame(
    Sequence_ID = c("s1", "s2"),
    VGENE_and_allele = c("IGHV1-1*01", "IGHV1-1*01"),
    JGENE_and_allele = c("IGHJ1*01", "IGHJ1*01"),
    JUNCTION = c("CARWGLRQEGFD", "CARWGLRQEGFE"),
    stringsAsFactors = FALSE
  )

  expect_error(
    clonality(df, cell = "B",
              cdr3_col = c("JUNCTION", "JUNCTION"),
              mismatch = c(0, 10, 20)),
    "must be length 1 or 3"
  )
})

test_that("extra_bin_cols prevents compensating-indel cells from sharing a bin", {
  # Two synthetic 'paired cells' (heavy_light concatenated) where total CDR3
  # length is identical but per-chain lengths differ. Without extra_bin_cols
  # these would share a bin and (with sufficient mismatch) cluster together.
  df <- data.frame(
    Sequence_ID = c("a", "b"),
    VGENE_and_allele = c("IGHV1-1*01_IGKV1-1*01", "IGHV1-1*01_IGKV1-1*01"),
    JGENE_and_allele = c("IGHJ1*01_IGKJ1*01", "IGHJ1*01_IGKJ1*01"),
    # Cell a: heavy length 12, light length 10 → concatenated "12_10"
    # Cell b: heavy length 13, light length 9  → concatenated "13_9"
    # Both have same total chars after the _ separator (12+1+10 == 13+1+9 == 23)
    JUNCTION        = c("CARAAAAAAAAA_CQQAAAAAAA",
                        "CARAAAAAAAAAB_CQQAAAAAA"),
    cdr3_length     = c("12_10", "13_9"),
    stringsAsFactors = FALSE
  )

  # Without extra_bin_cols: bin key is V_J_total-length → same bin → can merge
  # at any tolerant mismatch (here mm=20 normalised CDR3 distance ≈ 2/23 ≈ 8.7%).
  permissive <- clonality(df, cell = "B", cdr3_col = "JUNCTION", mismatch = 20)
  expect_equal(length(unique(permissive$clonality)), 1L)

  # With extra_bin_cols = "cdr3_length": the per-chain length string is part
  # of the bin key → different bins → cannot share a clonality ID.
  strict <- clonality(df, cell = "B", cdr3_col = "JUNCTION",
                      mismatch = 20, extra_bin_cols = "cdr3_length")
  expect_equal(length(unique(strict$clonality)), 2L)
})

test_that("cdr3_tags override controls the output column names", {
  df <- data.frame(
    Sequence_ID = c("s1", "s2"),
    VGENE_and_allele = c("IGHV1-1*01", "IGHV1-1*01"),
    JGENE_and_allele = c("IGHJ1*01", "IGHJ1*01"),
    JUNCTION = c("CARWGLRQEGFD", "CARWGLRQEGFE"),
    stringsAsFactors = FALSE
  )

  res <- clonality(df, cell = "B",
                   cdr3_col = "JUNCTION",
                   mismatch = c(0, 25),
                   cdr3_tags = c("strict", "loose"))

  expect_true(all(c("clonality_strict", "clonality_loose") %in% names(res)))
})
