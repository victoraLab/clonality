# data-raw/chain_classes.R

chain_classes <- list(
  unique_paired = list(
    BigCell.classes = c(
      "IGH_IGK", "IGH_IGL", "IGH_IGK_IGK", "IGH_IGL_IGL",
      "IGH_IGH_IGK", "IGH_IGH_IGL", "IGH_IGH_IGK_IGK"
    ),
    TabCell.classes = c(
      "TRA_TRB", "TRA_TRA_TRB", "TRA_TRB_TRB", "TRA_TRA_TRB_TRB"
    ),
    TgdCell.classes = c(
      "TRD_TRG", "TRD_TRD_TRG", "TRD_TRG_TRG", "TRD_TRD_TRG_TRG"
    )
  ),
  sticky_ends = list(
    BigCell.classes = c(
      "IGH", "IGK", "IGL", "IGH_IGH", "IGH_IGK", "IGH_IGL",
      "IGK_IGK", "IGL_IGL", "IGH_IGK_IGK", "IGH_IGL_IGL",
      "IGH_IGK_IGL", "IGH_IGH_IGK", "IGH_IGH_IGL", "IGH_IGH_IGK_IGK"
    ),
    TabCell.classes = c(
      "TRA", "TRB", "TRA_TRB", "TRA_TRA", "TRB_TRB",
      "TRA_TRA_TRB", "TRA_TRB_TRB", "TRA_TRA_TRB_TRB"
    ),
    TgdCell.classes = c(
      "TRD", "TRG", "TRD_TRG", "TRD_TRD", "TRG_TRG",
      "TRD_TRD_TRG", "TRD_TRG_TRG", "TRD_TRD_TRG_TRG"
    )
  )
)

usethis::use_data(chain_classes, overwrite = TRUE)
