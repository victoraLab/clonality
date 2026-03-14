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

test_that("assigntenx does not auto-populate file output paths", {
  mock_list_pairs <- list(
    IGH_IGK = list(
      cell1 = data.frame(
        v_gene = c("IGHV1-1", "IGKV1-1"),
        j_gene = c("IGHJ1", "IGKJ1"),
        c_gene = c("IGHM", "IGKC"),
        cdr3 = c("CARDR", "CQQY"),
        cdr3_nt = c("ATGCGA", "CCAAAT"),
        chain = c("IGH", "IGK"),
        raw_clonotype_id = c("raw1", "raw1"),
        barcode = c("cell1", "cell1"),
        stringsAsFactors = FALSE
      ),
      cell2 = data.frame(
        v_gene = c("IGHV1-1", "IGKV1-1"),
        j_gene = c("IGHJ1", "IGKJ1"),
        c_gene = c("IGHM", "IGKC"),
        cdr3 = c("CARER", "CQQF"),
        cdr3_nt = c("ATGCAA", "CCAATT"),
        chain = c("IGH", "IGK"),
        raw_clonotype_id = c("raw2", "raw2"),
        barcode = c("cell2", "cell2"),
        stringsAsFactors = FALSE
      )
    )
  )

  old_wd <- setwd(tempdir())
  on.exit(setwd(old_wd), add = TRUE)

  out <- assigntenx(
    list.pairs = mock_list_pairs,
    method = "unique_paired",
    experiment = "experiment1",
    clonality_input = NULL,
    cell = "B",
    col_res = "reduced",
    add_columns = NULL
  )

  expect_true("IGH_IGK" %in% names(out))
  expect_equal(length(list.files(pattern = "\\.xlsx$")), 0L)
})

