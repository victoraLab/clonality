## code to prepare `Cellranger9_BCR1` dataset goes here
Cellranger9_BCR1 <- read_delim("data-raw/TC073025.all_contig_annotations.csv")
usethis::use_data(Cellranger9_BCR1, overwrite = TRUE)
