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




