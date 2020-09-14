# Clonality

The Clonality package function helps to define B or T cell repertoire clonality automatically.

# Instalation:
```R
devtools::install_github("victoraLab/clonality")
```

# How to use:

The tra and trb objects are small sequencing datasets that can be used as a test

```R
head(tra)

clonality(File = tra, NewF = "tra.output")
```


On B cells:
Running the function on a xlsx file containing BCR sequences and getting clonal subclusters (Mismatch parameter) that have up to 20% CDR3 differences based on the Optimal string aligment, (restricted Damerau-Levenshtein distance from stringdist package).

```R
clonality(File = "example.xlsx",
ID_Column = "Sequence_ID",
Vgene_Column = "V_GENE_and_allele",
Jgene_Column = "J_GENE_and_allele",
CDR3_Column = "AA_JUNCTION",
cell = "B",
Mismatch = 0.2)
```

On T cells:
Running the function on a xlsx file containing TCR sequences and without subclustering, only identical sequences based on V gene, J gene, and CDR3 sequence will cluster.

```R
clonality(File = "example.xlsx",
ID_Column = "Sequence_ID",
Vgene_Column = "V_GENE_and_allele",
Jgene_Column = "J_GENE_and_allele",
CDR3_Column = "AA_JUNCTION",
cell = "T",
Mismatch = 0)
```

The function can also be used on a data.frame like:

```R
clonality(File = my_tcr, NewF = "output")
```


