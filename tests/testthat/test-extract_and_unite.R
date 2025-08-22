test_that("extract_and_unite returns concatenated strings with prefix match", {
  df <- tibble::tibble(
    v_gene_1 = c("V1", "V2"),
    v_gene_2 = c("V3", "V4"),
    unrelated = c("x", "y")
  )

  result <- extract_and_unite(df, "v_gene", unite = TRUE)
  expect_equal(result, c("V1_V3", "V2_V4"))
})

test_that("extract_and_unite returns original columns when unite = FALSE", {
  df <- tibble::tibble(
    j_gene_1 = c("J1", "J2"),
    j_gene_2 = c("J3", "J4")
  )

  result <- extract_and_unite(df, "j_gene", unite = FALSE)
  expect_s3_class(result, "tbl_df")
  expect_equal(ncol(result), 2)
  expect_equal(result$j_gene_1, c("J1", "J2"))
})

test_that("extract_and_unite works with regex matching", {
  df <- tibble::tibble(
    cdr3nt_IGH = c("A", "B"),
    cdr3nt_IGK = c("C", "D")
  )

  result <- extract_and_unite(df, "cdr3nt", unite = TRUE, use_regex = TRUE)
  expect_equal(result, c("A_C", "B_D"))
})

test_that("extract_and_unite handles no matching columns", {
  df <- tibble::tibble(a = 1:2)

  result <- suppressWarnings(extract_and_unite(df, "zzz", unite = TRUE))
  expect_equal(result, character(0))

  result <- suppressWarnings(extract_and_unite(df, "zzz", unite = FALSE))
  expect_s3_class(result, "tbl_df")
  expect_equal(ncol(result), 0)
})

