test_that("analyze_well_contigs works for WELL_PLATE and PLATE_WELL formats", {
  # Mock dataset (minimal example)
  df_mock <- data.frame(
    cell_id = c("SampleA01P01_1-100", "SampleA01P01_2-50", "SampleA01P02_1-200"),
    clonality = c("10.1", "10.1", "20.1"),
    stringsAsFactors = FALSE
  )

  df_mock2 <- data.frame(
    cell_id = c("SampleP01A01_1-100", "SampleP01A01_2-50", "SampleP02A01_1-200"),
    clonality = c("10.1", "10.1", "20.1"),
    stringsAsFactors = FALSE
  )

  # Test WELL_PLATE format
  result_wp <- analyze_well_contigs(df_mock, cell_id_col = "cell_id", barcode_format = "WELL_PLATE")
  expect_type(result_wp, "list")
  expect_named(result_wp, c("collapsed", "real_contigs_input"))
  expect_true(all(c("PLATE", "WELL", "CONTIG_NUMBER", "CONTIG_DEPTH", "background_tag") %in% colnames(result_wp$collapsed)))

  # Test PLATE_WELL format
  result_pw <- analyze_well_contigs(df_mock2, cell_id_col = "cell_id", barcode_format = "PLATE_WELL")
  expect_type(result_pw, "list")
  expect_named(result_pw, c("collapsed", "real_contigs_input"))
  expect_true(all(c("PLATE", "WELL", "CONTIG_NUMBER", "CONTIG_DEPTH", "background_tag") %in% colnames(result_pw$collapsed)))
})

test_that("analyze_well_contigs classifies REAL, BG1, and BG2 correctly", {
  df_mock <- data.frame(
    cell_id = c(
      "SampleA01P01_1-100",
      "SampleA01P01_2-40",
      "SampleA01P01_3-20",
      "SampleB01P01_1-60"
    ),
    clonality = c("10.1", "10.2", "20.1", "30.1"),
    stringsAsFactors = FALSE
  )

  result <- analyze_well_contigs(df_mock, barcode_format = "WELL_PLATE")
  tagged <- result$collapsed[result$collapsed$PLATE_WELL == "A01P01", ]

  expect_equal(tagged$background_tag[match("10.1", tagged$clonality)], "REAL")
  expect_equal(tagged$background_tag[match("10.2", tagged$clonality)], "BG1")
  expect_equal(tagged$background_tag[match("20.1", tagged$clonality)], "BG2")

  expect_equal(nrow(result$real_contigs_input), 2L)
  expect_true(all(c("SampleA01P01_1-100", "SampleB01P01_1-60") %in% result$real_contigs_input$cell_id))
})
