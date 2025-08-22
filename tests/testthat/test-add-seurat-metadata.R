test_that("add_seurat_metadata() errors if seu is NULL", {
  expect_error(add_seurat_metadata(NULL), "Please add Seurat object to merge clonality.")
})

