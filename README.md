# Integrative Inherited Metabolic Disease Gene and Model Organism Explorer

An interactive **R Shiny** application for exploring inherited metabolic diseases (IMDs), associated human genes, orthology relationships, published model-organism evidence, and pathway-level conservation across commonly used experimental organisms.

## Web application

The public web interface is available at:

**[Integrative IMD Model Organism Explorer](https://www.johnson-flylab.com/tools/integrative-imd-model-organism-explorer)**

## Overview

Inherited metabolic diseases are individually rare and many have limited clinical and mechanistic evidence. Selecting an appropriate experimental model can therefore be difficult, particularly when orthology relationships are complex.

This application integrates disease annotation, gene-level orthology, published model-organism studies, and KEGG pathway information to help researchers:

- search for IMDs by gene or disease name;
- compare orthology support across model organisms;
- interpret one-to-one and complex orthology relationships;
- browse IMDs by disease classification or pathway;
- identify published model-organism studies;
- explore pathway-level gene coverage; and
- prioritise candidate model organisms for further investigation.

Orthology information is intended to support model selection rather than replace biological judgement. Gene function, phenotype, disease relevance, experimental feasibility, and tissue or pathway context should also be considered.

## Supported model organisms

The application currently includes:

| Common name | Scientific name |
|---|---|
| Mouse | *Mus musculus* |
| Zebrafish | *Danio rerio* |
| Fruit fly | *Drosophila melanogaster* |
| Nematode worm | *Caenorhabditis elegans* |
| Budding yeast | *Saccharomyces cerevisiae* |
| Fission yeast | *Schizosaccharomyces pombe* |

## Main features

### Disease search

Search by:

- human gene symbol;
- ICIEM disease name;
- alternative disease name; or
- disease abbreviation.

Multiple terms can be entered, one per line. The search includes partial and approximate matching, so results should be reviewed before an entry is selected.
<p align="center">
  <img src="README_Figure/Disease%20search.png" alt="Disease search interface" width="1200">
</p>
<p align="center"><em>Example of the disease-search interface.</em></p>


### Adjustable DIOPT threshold

The advanced search panel provides a global DIOPT score threshold. The default threshold is **7**.

Increasing the threshold retains more strongly supported orthology predictions. Decreasing it includes broader, potentially weaker orthology support. The selected threshold is applied throughout the current session, including disease searches, model browsing, disease-classification browsing, pathway browsing, and disease-page outputs.
<p align="center">
  <img src="README_Figure/Adjustable%20DIOPT%20threshold.png" alt="Adjustable DIOPT score threshold" width="1000">
</p>
<p align="center"><em>The global DIOPT threshold control used throughout the application.</em></p>


### Browse by model organism

Users can browse IMD-associated genes according to orthology support in each model organism and distinguish:

- one-to-one orthology;
- complex orthology; and
- no detected ortholog support.
<p align="center">
  <img src="README_Figure/Browse%20by%20model%20organism.png" alt="Browse by model organism page" width="1200">
</p>
<p align="center"><em>Model-organism overview showing one-to-one and total ortholog support.</em></p>


### ICIMD disease classification

IMD records can be explored using the ICIMD classification hierarchy. Each disease category displays model-organism orthology coverage and can be expanded to show subcategories and individual disease records.
<p align="center">
  <img src="README_Figure/ICIMD%20disease%20classification_instruction.png" alt="ICIMD classification instructions" width="700">
</p>
<p align="center"><em>In-app instructions for browsing the ICIMD disease classification.</em></p>

<p align="center">
  <img src="README_Figure/ICIMD%20disease%20classification.png" alt="ICIMD disease classification interface" width="1100">
</p>
<p align="center"><em>Example of an expanded ICIMD disease category.</em></p>


### IMD–KEGG pathway classification

IMD-associated genes can be explored by KEGG pathway category, subcategory, and pathway. Interactive pathway maps allow users to inspect mapped genes and open linked disease records.
<p align="center">
  <img src="README_Figure/IMD%E2%80%93KEGG%20pathway%20classification_instruction%20step%200.png" alt="IMD-KEGG classification instructions" width="700">
</p>
<p align="center"><em>In-app instructions for using the IMD–KEGG pathway classification.</em></p>

<p align="center">
  <img src="README_Figure/IMD%E2%80%93KEGG%20pathway%20classification.png" alt="IMD-KEGG pathway classification interface" width="1100">
</p>
<p align="center"><em>Browse IMD-associated genes by KEGG pathway category and pathway.</em></p>

<p align="center">
  <img src="README_Figure/IMD%E2%80%93KEGG%20pathway%20classification_Pathway%20map%20step1.png" alt="Model selection for a KEGG pathway" width="1200">
</p>
<p align="center"><em>Select a model organism to assess pathway-level orthology coverage.</em></p>

<p align="center">
  <img src="README_Figure/IMD%E2%80%93KEGG%20pathway%20classification_Pathway%20map%20step2.png" alt="Interactive IMD-KEGG pathway map" width="1100">
</p>
<p align="center"><em>Example pathway map showing IMD genes with and without ortholog support.</em></p>


### IMD disease pages

Each disease page can display:

- ICIEM disease name;
- associated human gene;
- disease category and subcategory;
- alternative names and abbreviations;
- mode of inheritance;
- treatability;
- prevalence;
- OMIM, Orphanet, IEMbase, and GeneReviews links;
- published IMD studies in model organisms;
- orthology relationships across all supported organisms;
- DIOPT-supported human and model-gene components; and
- links to organism-specific model databases.
<p align="center">
  <img src="README_Figure/IMD%20disease%20pages.png" alt="IMD disease page" width="1200">
</p>
<p align="center"><em>Example disease page containing disease annotation and published model-organism evidence.</em></p>

<p align="center">
  <img src="README_Figure/IMD%20disease%20pages-pathway%20map.png" alt="Pathways associated with an IMD gene" width="1000">
</p>
<p align="center"><em>Pathway links associated with the selected IMD-causing gene.</em></p>


### Orthology-based model prioritisation

The application uses the following interpretation:

| Orthology relationship | Application label | Interpretation |
|---|---|---|
| 1:1 | Prioritized | Clearest gene-level correspondence |
| 1:2 | Retained | Potentially useful, with possible paralogue redundancy |
| Other complex relationships | Deprioritized | Requires additional biological interpretation |
| No supported relationship | No ortholog | No ortholog support at the selected threshold |

When multiple models have one-to-one orthology, a model-selection matrix is used to rank them according to the configured model weights.
<p align="center">
  <img src="README_Figure/Orthology-based%20model%20prioritisation.png" alt="Orthology-based model prioritisation output" width="1000">
</p>
<p align="center"><em>Example of the orthology-based model-selection output for an IMD-causing gene.</em></p>


## Repository structure

The application expects the following general structure:

```text
Model-Organism-IMD-Explorer/
├── app.R
├── prepare_data.R
├── app_data.rds
├── pathway_gene_clickmap.rds
│
├── KEGG_hsa_pathway_to_genes_with_category.xlsx
├── IMDs_Info_With_Category.xlsx
├── IMDs_Others.xlsx
├── Gene_reference_table.xlsx
├── Childhood_Dementia_Gene_List.xlsx
├── Disease_Associations_With_Paper_Titles.docx
├── model_selection_matrix.xlsx
├── pathway_gene_clickmap_with_image_size.csv
│
├── Human_to_Mouse_Orthology_State_All_DIOPT_Thresholds.xlsx
├── Human_to_Zebrafish_Orthology_State_All_DIOPT_Thresholds.xlsx
├── Human_to_Fly_Orthology_State_All_DIOPT_Thresholds.xlsx
├── Human_to_Worm_Orthology_State_All_DIOPT_Thresholds.xlsx
├── Human_to_BuddingYeast_Orthology_State_All_DIOPT_Thresholds.xlsx
├── Human_to_FissionYeast_Orthology_State_All_DIOPT_Thresholds.xlsx
│
├── Model_Link/
│   ├── Mouse.xlsx
│   ├── Zebrafish.xlsx
│   ├── Fly.xlsx
│   ├── Worm.xlsx
│   ├── BuddingYeast.xlsx
│   └── FissionYeast.xlsx
│
├── README_Figure/
│   └── screenshots used in README.md
│
├── Pathway_View/
│   ├── Mus musculus/
│   ├── Danio rerio/
│   ├── Drosophila melanogaster/
│   ├── Caenorhabditis elegans/
│   ├── Saccharomyces cerevisiae/
│   └── Schizosaccharomyces pombe/
│
└── www/
    └── static images and interface assets
```

The exact contents of `Pathway_View/` and `www/` depend on the pathway images and interface assets included in the deployed version.

## Requirements

The application uses the following R packages:

```r
required_packages <- c(
  "shiny",
  "shinyjs",
  "DT",
  "readxl",
  "dplyr",
  "stringr",
  "tidyr",
  "purrr",
  "htmltools",
  "scales",
  "officer",
  "tibble"
)

new_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}
```

## Data preparation

`prepare_data.R` reads and validates the source files, standardises gene symbols, combines disease annotations, processes orthology results at each DIOPT threshold, imports model-database links, calculates model-selection information, and prepares pathway click-map coordinates.

It creates two deployment-ready files:

```text
app_data.rds
pathway_gene_clickmap.rds
```

### Important path setting

Before running the preparation script, make sure `base_dir` points to the repository directory.

For a portable GitHub version, use:

```r
base_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
```

Then open R in the repository root and run:

```r
source("prepare_data.R")
```

Alternatively, from a terminal:

```bash
Rscript prepare_data.R
```

The script should finish with messages confirming that `app_data.rds` and `pathway_gene_clickmap.rds` were created.

## Running the application locally

Clone the repository:

```bash
git clone https://github.com/johnsonflygroup/Model-Organism-IMD-Explorer.git
cd Model-Organism-IMD-Explorer
```

Start the application from R:

```r
shiny::runApp()
```

The app requires the following files in the repository root before startup:

```text
app_data.rds
pathway_gene_clickmap.rds
```

It also requires the pathway images and static interface assets used by the deployed version.

## Updating the data

When an input table or source file is changed:

1. replace or update the relevant source file;
2. confirm that its filename and required column names remain compatible with `prepare_data.R`;
3. rerun `prepare_data.R`;
4. check that both RDS files are recreated successfully;
5. launch the app locally and inspect several representative disease, model, and pathway records; and
6. commit the updated source files, scripts, and generated RDS files as appropriate.

## Minimum files needed for deployment

A deployment that does not rebuild the source data requires at least:

```text
app.R
app_data.rds
pathway_gene_clickmap.rds
Pathway_View/
www/
```

For full reproducibility, also retain `prepare_data.R` and the source data tables, subject to the reuse and redistribution conditions of the original data providers.

## Data resources

The application integrates or links information derived from resources including:

- IEMbase and the ICIMD disease classification;
- DIOPT orthology predictions;
- KEGG pathway annotations;
- Alliance of Genome Resources;
- OMIM;
- Orphanet;
- GeneReviews;
- Mouse Genome Informatics;
- ZFIN;
- FlyBase;
- WormBase;
- Saccharomyces Genome Database; and
- PomBase.

Users should consult and cite the original databases and publications when using information derived from these resources.

## Reproducibility notes

- Human gene symbols are standardised using `Gene_reference_table.xlsx`.
- Orthology results are organised by DIOPT threshold.
- The default application threshold is 7 when that value is available.
- The same selected threshold is applied across the application during a user session.
- The application reads preprocessed RDS objects at startup to reduce repeated parsing of Excel, CSV, and Word files.
- Pathway click-map coordinates are stored separately and loaded when required.

## Citation

A formal citation will be added following publication of the associated manuscript.

## Disclaimer

This resource is intended for research and model-selection support. It is not a diagnostic or clinical decision-making tool. Orthology support does not guarantee conservation of phenotype, biochemical function, tissue context, or therapeutic response.

## Contributing

Suggestions, bug reports, and reproducibility issues can be submitted through the repository's **Issues** page. When reporting an issue, please include:

- the affected gene, disease, model organism, or pathway;
- the selected DIOPT threshold;
- the expected and observed result; and
- steps needed to reproduce the problem.

## License
