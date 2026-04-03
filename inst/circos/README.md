# Clonality Circos Viewer

This folder contains the standalone HTML viewer bundled with the `clonality`
package.

## Open it from R

```r
clonality::open_circos_viewer()
```

## Expected input

The viewer accepts a long table with three columns:

```text
clonotype,size,group
1.1,2,M1LP_PC
1.1,11,M1PP_GC
2.1,3,M1LP_PC
6.1,3,M1LP_PC
6.1,1,M2LP_PC
```

## Prepare input from a wide table

```r
viewer_df <- clonality::circos_viewer_wide_to_long(wide_df, clone_col = "Clone_ID")
```

## Prepare input from clonality output

If you already have a `clonality` column and a sample or condition column:

```r
viewer_df <- clonality::circos_viewer_input(
  data = results,
  clonotype_col = "clonality",
  group_col = "condition"
)
```

If you have a numeric abundance column instead of one row per cell:

```r
viewer_df <- clonality::circos_viewer_input(
  data = results,
  clonotype_col = "clonality",
  group_col = "condition",
  size_col = "n"
)
```

Paste the resulting table into the viewer in `Clonotype table` mode.
