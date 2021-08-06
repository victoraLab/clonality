#clonality

The clonality package function helps to define B or T cell repertoire clonality automatically.

#Instalation:

```R
devtools::install_github("victoraLab/clonality")
```

# How to use:

The tra and trb objects are small sequencing datasets that can be used as a test

```R
head(tra)
clonality(data =  tra, output =  "tra.output")
```

On B cells:
Running the function on a xlsx file containing BCR sequences and getting clonal subclusters (mm (Mismatch) parameter) that have up to 20% CDR3 differences based on the string aligment, (restricted Damerau-Levenshtein distance from stringdist package).

```R
clonality(data = "example.xlsx",
id_col = "Sequence_ID",
vgene_col = "V_GENE_and_allele",
jgene_col = "J_GENE_and_allele",
cdr3_col = "AA_JUNCTION",
cell = "B",
mm = 0.2)
```

On T cells:
Running the function on a xlsx file containing TCR sequences and without subclustering, only identical sequences based on V gene, J gene, and CDR3 sequence will cluster.

```R
clonality(data = "example.xlsx",
id_col = "Sequence_ID",
vgene_col = "V_GENE_and_allele",
jgene_col = "J_GENE_and_allele",
cdr3_col = "AA_JUNCTION",
cell = "T",
mm = 0)
```

The function can also be used on a data.frame like:

```R
clonality(data = my_tcr, output = "output")
```

This package also has a function to call and import clonality on the `filtered_contig_annotations` file from 10x genomics experiments.

Use the `sticky` parameter = `TRUE` to map cells with only one chain to cells with the paired sequence. 

```R
tenx(data = filtered_contig_annotations, method = "sticky_ends", cell = "T", only_productive = T)
seu <- add_seurat_metadata(seu = singlets,  clonal.list = ls(pattern = "^Clonal"), sticky = T, purity = 0.8)
```


