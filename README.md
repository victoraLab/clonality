# Clonality

The Clonality package function helps to define B or T cell repertoire clonality automatically.

# Instalation
```R
devtools::install_github("victoraLab/clonality")
```

# How to use:

```R
Clonality(File = "test.bcell.xlsx", ID_Column = "Sequence_ID", Vgene_Column = "V_GENE_and_allele", Jgene_Column = "J_GENE_and_allele", CDR3_Column = "AA_JUNCTION", cell = "B")
```

