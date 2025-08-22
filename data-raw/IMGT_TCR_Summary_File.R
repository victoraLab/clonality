library(readr)

df <- parse_imgt("data-raw/RP230221//")

# Save cleaned data to package as IMGT_BCR_Summary_File
IMGT_TCR_Summary_File <- df
usethis::use_data(IMGT_TCR_Summary_File, overwrite = TRUE)
