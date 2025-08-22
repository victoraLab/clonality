library(readr)

df <- parse_imgt("data-raw/Carla082024/")

# Save cleaned data to package as IMGT_BCR_Summary_File
IMGT_BCR_Summary_File <- df
usethis::use_data(IMGT_BCR_Summary_File, overwrite = TRUE)
