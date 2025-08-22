## code to prepare `TCR_Example` dataset goes here
tra <- readr::read_csv("data-raw/tra.csv")
trb <- readr::read_csv("data-raw/trb.csv")
usethis::use_data(tra, overwrite = TRUE)
usethis::use_data(trb, overwrite = TRUE)
