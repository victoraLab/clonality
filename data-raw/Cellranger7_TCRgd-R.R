## code to prepare `Cellranger7_TCRgd.R` dataset goes here
Cellranger7_TCRgd <- readr::read_csv("data-raw/BR010622.all_contig_annotations.csv")
usethis::use_data(Cellranger7_TCRgd, overwrite = TRUE)
