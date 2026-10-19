# CP mitochondrial health – Olink inflammation analysis

This repository contains R code used for the analysis of circulating inflammatory proteins measured using the Olink Target 96 Inflammation panel in children with cerebral palsy (CP) and typically developing (TD) controls.

The analyses form part of a study investigating mitochondrial health and associated systemic alterations in skeletal muscle from children with cerebral palsy.

## Analysis workflow

The scripts are intended to be run sequentially:

### `00_load_and_clean_olink.R`
Imports and prepares the Olink dataset for downstream analyses.

The script:
- imports the raw Olink data,
- corrects sample identifiers,
- removes internal Olink control samples,
- assigns participants to CP and TD groups, and
- generates a sample × protein matrix containing normalized protein expression (NPX) values.

### `01_PCA_olink_inflammation.R`
Performs principal component analysis (PCA) of the Olink Target 96 Inflammation dataset.

Protein NPX values are standardized prior to PCA. The first two principal components are visualized together with 95% confidence ellipses for the CP and TD groups.

### `02_volcano_CP_vs_TD_olink.R`
Performs protein-wise comparisons between CP and TD groups and generates a volcano plot.

For each protein:
- CP and TD groups are compared using an independent-samples t-test,
- log2 fold change is calculated as mean NPX(CP) − mean NPX(TD), and
- p-values are adjusted for multiple comparisons using the Benjamini–Hochberg procedure.

Proteins are classified as differentially expressed when:

- adjusted p-value < 0.05, and
- |log2 fold change| ≥ 0.58 (approximately 1.5-fold).

Positive fold changes indicate higher protein levels in CP.

### `03_boxplots_selected_proteins.R`
Generates boxplots with individual participant data for selected proteins identified from the Olink analysis.

### `04_export_olink_summary_table.xlsx.R`
Exports a protein-level summary table containing:
- protein name,
- log2 fold change,
- Benjamini–Hochberg adjusted p-value, and
- individual participant NPX values.

## Software and R packages

Analyses were performed in R.

The following R packages are used:

- `dplyr`
- `tidyr`
- `readr`
- `ggplot2`
- `ggrepel`
- `purrr`
- `writexl`

## Data availability

The repository contains analysis code but does not contain participant-level source data.

Access to the underlying human research data is subject to applicable ethical approvals and data-protection requirements.

## Authors

Sebastian Edman  
Karolinska Institutet, Stockholm, Sweden

## Citation

If you use this code, please cite the associated publication:

[Publication details to be added upon publication.]
