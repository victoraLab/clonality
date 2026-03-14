test_that("parse_imgt merges summary and junction files and cleans names", {
  path <- tempfile()
  dir.create(path)

  writeLines(
    c(
      "Sequence ID\tVGENE and allele\tExtra Column",
      "seq1\tIGHV1-1*01\tfoo",
      "seq2\tIGHV1-2*01\tbar"
    ),
    file.path(path, "1_Summary.txt")
  )

  writeLines(
    c(
      "Sequence ID\tJUNCTION",
      "seq1\tATGC",
      "seq2\tTTAA"
    ),
    file.path(path, "6_Junction.txt")
  )

  result <- parse_imgt(path)

  expect_true(all(c("Sequence_ID", "VGENE_and_allele", "Extra_Column", "DNA_Junction") %in% names(result)))
  expect_equal(result$DNA_Junction, c("ATGC", "TTAA"))
})

test_that("parse_imgt errors when Sequence_ID is missing", {
  path <- tempfile()
  dir.create(path)

  writeLines(
    c(
      "Bad Column\tValue",
      "seq1\tfoo"
    ),
    file.path(path, "1_Summary.txt")
  )

  writeLines(
    c(
      "Sequence ID\tJUNCTION",
      "seq1\tATGC"
    ),
    file.path(path, "6_Junction.txt")
  )

  expect_error(parse_imgt(path), "required 'Sequence_ID' column is missing")
})

test_that("parse_imgt coalesces existing DNA_Junction without duplicate names", {
  path <- tempfile()
  dir.create(path)

  writeLines(
    c(
      "Sequence ID\tDNA Junction\tVGENE and allele",
      "seq1\tATGC\tIGHV1-1*01",
      "seq2\t\tIGHV1-2*01"
    ),
    file.path(path, "1_Summary.txt")
  )

  writeLines(
    c(
      "Sequence ID\tJUNCTION",
      "seq1\tATGC",
      "seq2\tTTAA"
    ),
    file.path(path, "6_Junction.txt")
  )

  result <- parse_imgt(path)

  expect_equal(sum(names(result) == "DNA_Junction"), 1L)
  expect_equal(result$DNA_Junction, c("ATGC", "TTAA"))
})
