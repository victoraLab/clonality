## code to prepare `Cellranger7_BCR` dataset goes here
Cellranger7_BCR1 <- readr::read_csv("data-raw/SN082621.all_contig_annotations.csv")
Cellranger7_BCR2 <- readr::read_csv("data-raw/SN062424.all_contig_annotations.csv")
usethis::use_data(Cellranger7_BCR1, overwrite = TRUE)
usethis::use_data(Cellranger7_BCR2, overwrite = TRUE)
