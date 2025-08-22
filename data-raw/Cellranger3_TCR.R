## code to prepare `Cellranger3_TCR` dataset goes here
Cellranger3_TCR <- readr::read_csv("data-raw/AB070119.all_contig_annotations.csv")

usethis::use_data(Cellranger3_TCR, overwrite = TRUE)
