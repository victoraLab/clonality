# Clonality Sharing Circos Viewer

Standalone HTML tool for the lab to visualize BCR/TCR clonality sharing between groups with D3.

## Use

Open [index.html](/path/to/project) in a browser, choose an input mode, paste data, and click `Render Plot`.

## Input modes

### Clonotype table

This mode is closer to the old Circos workflow.

- Paste `id,size,group` style data
- Clonotypes are sorted from larger to smaller within each group
- Shared clonotypes get ribbons
- Unique clonotypes remain visible as ideogram-like bands with no ribbons

Example:

```text
clonotype,size,group
clA,40,Blood
clB,28,Blood
clA,36,Tumor
clC,20,Tumor
clB,24,Spleen
clC,16,Spleen
```

### Sharing matrix

- Labeled CSV / TSV with column labels in the first row and row labels in the first column
- Numeric-only square matrix with auto-generated group names
- Optional label override entered separately as one label per line or comma-separated

## Example matrix

```text
Group,Blood,Tumor,Spleen,LN
Blood,18,7,3,5
Tumor,6,22,8,9
Spleen,2,7,15,4
LN,4,10,5,19
```

## Export

The page can download the current figure as:

- `SVG`
- `PNG`
