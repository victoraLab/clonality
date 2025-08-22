test_that("coalescing raw_clonotype_id and barcode columns works correctly", {
  # Simulate wide-format row1() output for 2 cells
  mock_res <- tibble::tibble(
    raw_clonotype_id_IGH = c("cloneA", NA),
    raw_clonotype_id_IGK = c(NA, "cloneB"),
    barcode_IGH = c("cellA", NA),
    barcode_IGK = c(NA, "cellB")
  )

  # Step 1: Coerce all matching columns to character
  mock_res <- mock_res %>%
    dplyr::mutate_at(dplyr::vars(dplyr::matches("raw")), as.character) %>%
    dplyr::mutate(sc.raw_clonotypes = dplyr::coalesce(!!!dplyr::select(., dplyr::matches("raw"))))

  mock_res <- mock_res %>%
    dplyr::mutate_at(dplyr::vars(dplyr::matches("bar")), as.character) %>%
    dplyr::mutate(sc.barcodes = dplyr::coalesce(!!!dplyr::select(., dplyr::matches("bar"))))

  # Expectations
  expect_equal(mock_res$sc.raw_clonotypes, c("cloneA", "cloneB"))
  expect_equal(mock_res$sc.barcodes, c("cellA", "cellB"))

  # Check that output columns exist
  expect_true(all(c("sc.raw_clonotypes", "sc.barcodes") %in% colnames(mock_res)))
})




