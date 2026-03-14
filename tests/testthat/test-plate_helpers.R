test_that("extract_plate_metadata parses both barcode formats and auto-detects", {
  df_wp <- data.frame(cell_id = "SampleA01P09_2-150", stringsAsFactors = FALSE)
  df_pw <- data.frame(cell_id = "SampleP09A01_2-150", stringsAsFactors = FALSE)

  res_wp <- extract_plate_metadata(df_wp, barcode_format = "WELL_PLATE")
  res_pw <- extract_plate_metadata(df_pw, barcode_format = "auto")

  expect_equal(res_wp$PLATE, "P09")
  expect_equal(res_wp$WELL, "A01")
  expect_equal(res_wp$CONTIG_NUMBER, 2L)
  expect_equal(res_wp$CONTIG_DEPTH, 150L)
  expect_equal(res_wp$ggplate_well, "A1")

  expect_equal(res_pw$PLATE_WELL, "P09A01")
  expect_equal(res_pw$PLATE, "P09")
  expect_equal(res_pw$WELL, "A01")
})

test_that("extract_plate_metadata errors for malformed IDs", {
  df_bad <- data.frame(cell_id = "not_a_valid_id", stringsAsFactors = FALSE)

  expect_error(
    extract_plate_metadata(df_bad, barcode_format = "WELL_PLATE"),
    "Failed to parse some cell IDs"
  )
})

test_that("evaluate_plate_quality summarizes per-plate metrics", {
  skip_if_not_installed("MASS")

  analysis_output <- list(
    collapsed = tibble::tibble(
      PLATE_WELL = c("A01P01", "A02P01", "A03P01", "A04P01", "A05P01"),
      PLATE = c("P01", "P01", "P01", "P01", "P01"),
      WELL = c("A01", "A02", "A03", "A04", "A05"),
      CONTIG_NUMBER = c(1L, 1L, 1L, 1L, 1L),
      CONTIG_DEPTH = c(100L, 80L, 60L, 40L, 20L),
      TREE = c("1", "2", "3", "4", "5"),
      collapsed = c(TRUE, FALSE, FALSE, TRUE, FALSE),
      background_tag = c("REAL", "REAL", "REAL", "REAL", "BG1")
    ),
    real_contigs_input = tibble::tibble(dummy = 1)
  )

  result <- evaluate_plate_quality(analysis_output, plate_size = 4)

  expect_s3_class(result, "tbl_df")
  expect_equal(result$PLATE, "P01")
  expect_equal(result$n_wells_detected, 4)
  expect_equal(result$n_successful_contigs, 4)
  expect_equal(result$n_fully_collapsed, 2)
  expect_equal(result$efficiency_score, 100)
  expect_equal(result$collapse_quality, 0.5)
  expect_false(is.na(result$mu))
  expect_false(is.na(result$size))
})

test_that("evaluate_plate_quality validates expected input structure", {
  expect_error(
    evaluate_plate_quality(list(collapsed = tibble::tibble())),
    "Input must be a list returned by `analyze_well_contigs\\(\\)`"
  )
})
