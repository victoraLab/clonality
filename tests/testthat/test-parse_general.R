test_that("parse_general returns data frames unchanged", {
  input <- data.frame(a = 1:2, b = c("x", "y"))

  result <- parse_general(input)

  expect_equal(result, input)
})

test_that("parse_general reads csv and txt files", {
  csv_path <- tempfile(fileext = ".csv")
  txt_path <- tempfile(fileext = ".txt")

  writeLines(c("a,b", "1,x", "2,y"), csv_path)
  writeLines(c("a\tb", "1\tx", "2\ty"), txt_path)

  csv_res <- parse_general(csv_path)
  txt_res <- parse_general(txt_path)

  expect_s3_class(csv_res, "tbl_df")
  expect_s3_class(txt_res, "tbl_df")
  expect_equal(csv_res$a, c(1, 2))
  expect_equal(txt_res$b, c("x", "y"))
})

test_that("parse_general rejects unsupported extensions and invalid inputs", {
  bad_path <- tempfile(fileext = ".tsv")
  writeLines("a\tb", bad_path)

  expect_error(parse_general(bad_path), "Unsupported file format")
  expect_error(parse_general(123), "Input must be a data frame or a file path")
})
