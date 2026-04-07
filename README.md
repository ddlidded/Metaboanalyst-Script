# Metabolomics Analysis Shiny GUI

This repository contains a Shiny application (`app.R`) for interactive metabolomics analysis inspired by common MetaboAnalyst workflows.

## Features

- Upload metabolomics input as:
  - CSV (`.csv`, `.txt`)
  - Excel (`.xls`, `.xlsx`)
- Select metadata columns:
  - Sample ID column (optional)
  - Group/Class column (optional)
- Choose metabolite columns to analyze
- Optional preprocessing:
  - Missing value handling (none, half-minimum imputation, remove rows with missing values)
  - Data transformation (none, `log2`, `log10`, `sqrt`)
  - Scaling/normalization (none, auto, pareto, range 0-1)
- Select which analyses/visualizations to generate:
  - PCA score plot
  - Heatmap
  - Correlation matrix
  - Volcano plot (requires exactly 2 groups)
  - Top-feature boxplots
  - t-test table (requires exactly 2 groups)
  - ANOVA table (requires 3+ groups)
- Export outputs:
  - ZIP of selected analysis plots as PNG
  - ZIP of per-metabolite bar plots (one PNG per metabolite)
  - ZIP of generated statistical result tables as CSV

## Required R packages

Install dependencies in R:

```r
install.packages(c("shiny", "ggplot2", "dplyr", "tidyr", "readxl"))
```

`readxl` is only needed when uploading Excel files.

## Run the app

From the repository root:

```bash
Rscript -e "shiny::runApp('app.R')"
```

Or inside an interactive R session:

```r
shiny::runApp("app.R")
```

## Sample file list (import templates)

Use the files in `sample_data/` as formatting references before importing your own dataset:

- `sample_data/metabolomics_two_groups.csv`
  - Includes `SampleID`, `Group` (2 classes: Control/Treatment), and metabolite columns.
  - Good for: PCA, heatmap, correlation, volcano plot, t-test, top-feature boxplots, per-metabolite bar plots.

- `sample_data/metabolomics_three_groups.csv`
  - Includes `SampleID`, `Group` (3 classes), and metabolite columns.
  - Good for: PCA, heatmap, correlation, ANOVA, top-feature boxplots, per-metabolite bar plots.

- `sample_data/metabolomics_no_group.csv`
  - Includes `SampleID` and metabolite columns only (no group/class column).
  - Good for: PCA, heatmap, correlation, per-metabolite bar plots.
  - Group-dependent analyses (volcano, t-test, ANOVA) will be skipped.

### Required format rules

- **One row = one sample**
- **One column = one field**
  - Metadata columns (e.g., `SampleID`, `Group`) are optional but recommended.
  - Metabolite columns should be numeric (or numeric-like strings).
- Keep the first row as column headers.
- Avoid merged cells and free-text notes inside the data table.

### Using Excel files

If you use Excel (`.xls/.xlsx`), mirror the same column layout as the sample CSV files:

- same header names
- same row-wise sample layout
- one sample per row

## Expected data shape

- Rows should represent samples.
- Columns should include:
  - one or more metabolite abundance columns (numeric or numeric-like)
  - optional sample and/or group metadata columns

When no group column is selected, analyses that require grouping (volcano, t-test, ANOVA by group) are skipped with a message.
