# prepare_data.R
# Run this locally before deployment
# This creates app_data.rds so app.R does not need to read all Excel/CSV/Word files at startup.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(officer)
  library(tibble)
})

# -------------------------
# 1) PATHS
# -------------------------

base_dir <- "C:/Users/a2618/Desktop/Shiny"
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)

file_kegg <- file.path(base_dir, "KEGG_hsa_pathway_to_genes_with_category.xlsx")
file_info <- file.path(base_dir, "IMDs_Info_With_Category.xlsx")
file_other <- file.path(base_dir, "IMDs_Others.xlsx")
file_gene_reference <- file.path(base_dir, "Gene_reference_table.xlsx")

file_childhood_dementia <- file.path(base_dir, "Childhood_Dementia_Gene_List.xlsx")
file_research_doc <- file.path(base_dir, "Disease_Associations_With_Paper_Titles.docx")

file_mouse <- file.path(base_dir, "Human_to_Mouse_Orthology_State_All_DIOPT_Thresholds.xlsx")
file_zfish <- file.path(base_dir, "Human_to_Zebrafish_Orthology_State_All_DIOPT_Thresholds.xlsx")
file_fly <- file.path(base_dir, "Human_to_Fly_Orthology_State_All_DIOPT_Thresholds.xlsx")
file_worm <- file.path(base_dir, "Human_to_Worm_Orthology_State_All_DIOPT_Thresholds.xlsx")
file_yeast <- file.path(base_dir, "Human_to_FissionYeast_Orthology_State_All_DIOPT_Thresholds.xlsx")
file_budding_yeast <- file.path(base_dir, "Human_to_BuddingYeast_Orthology_State_All_DIOPT_Thresholds.xlsx")

file_gene_clickmap <- file.path(base_dir, "Pathway_gene_clickmap_with_image_size.csv")

file_matrix <- file.path(base_dir, "model_selection_matrix.xlsx")

model_link_dir <- file.path(base_dir, "Model_Link")

file_model_mouse <- file.path(model_link_dir, "Mouse.xlsx")
file_model_zfish <- file.path(model_link_dir, "Zebrafish.xlsx")
file_model_fly <- file.path(model_link_dir, "Fly.xlsx")
file_model_worm <- file.path(model_link_dir, "Worm.xlsx")
file_model_yeast <- file.path(model_link_dir, "FissionYeast.xlsx")
file_model_budding_yeast <- file.path(model_link_dir, "BuddingYeast.xlsx")


# -------------------------
# 2) HELPER FUNCTIONS
# -------------------------

clean_text <- function(x) {
  x <- as.character(x)
  x <- ifelse(is.na(x), "", x)
  stringr::str_squish(x)
}

pick_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_col <- function(df, col, default = "") {
  if (is.na(col) || !(col %in% names(df))) {
    return(rep(default, nrow(df)))
  }
  df[[col]]
}

is_has_ortholog <- function(rel, model_gene_n = NULL) {
  rel2 <- tolower(clean_text(rel))
  rel2 <- gsub("\\s+", "", rel2)
  
  out <- rel2 != "" &
    !rel2 %in% c("1:0", "0:1", "0:0", "many:0", "0:many") &
    !stringr::str_detect(rel2, "no.*ortholog")
  
  if (!is.null(model_gene_n)) {
    model_gene_n <- suppressWarnings(as.numeric(model_gene_n))
    out <- out & !is.na(model_gene_n) & model_gene_n > 0
  }
  
  out
}

priority_from_relationship <- function(rel) {
  rel2 <- tolower(clean_text(rel))
  rel2 <- gsub("\\s+", "", rel2)
  
  dplyr::case_when(
    rel2 %in% c("", "1:0", "0:1", "0:0", "many:0", "0:many") |
      stringr::str_detect(rel2, "no.*ortholog") ~ "No ortholog",
    rel2 == "1:1" ~ "Prioritized",
    rel2 == "1:2" ~ "Retained",
    TRUE ~ "Deprioritized"
  )
}

priority_score <- function(priority_label) {
  dplyr::case_when(
    priority_label == "Prioritized" ~ 2,
    priority_label == "Retained" ~ 1,
    priority_label == "Deprioritized" ~ 0,
    TRUE ~ -1
  )
}

parse_research_docx <- function(docx_path) {
  if (!file.exists(docx_path)) {
    return(tibble(
      Gene = character(),
      Species = character(),
      Title = character(),
      Reference = character(),
      Status = character()
    ))
  }
  
  doc <- officer::read_docx(docx_path)
  txt <- officer::docx_summary(doc)$text
  
  txt <- clean_text(txt)
  txt <- txt[txt != ""]
  
  out <- list()
  
  current_gene <- NA_character_
  current_species <- NA_character_
  pending_title <- NA_character_
  
  for (line in txt) {
    
    if (grepl("^\\d+\\.\\s+[A-Za-z0-9._-]+$", line)) {
      current_gene <- sub("^\\d+\\.\\s+", "", line)
      current_species <- NA_character_
      pending_title <- NA_character_
      next
    }
    
    if (grepl("^No disease associations found with paper-like references for selected species\\.$", line)) {
      if (!is.na(current_gene)) {
        out[[length(out) + 1]] <- tibble(
          Gene = current_gene,
          Species = "",
          Title = "",
          Reference = "",
          Status = "unavailable"
        )
      }
      current_species <- NA_character_
      pending_title <- NA_character_
      next
    }
    
    if (grepl("^Species\\s+\\d+\\s*:\\s*", line)) {
      current_species <- sub("^Species\\s+\\d+\\s*:\\s*", "", line)
      pending_title <- NA_character_
      next
    }
    
    if (grepl("^Title\\s*:\\s*", line)) {
      pending_title <- sub("^Title\\s*:\\s*", "", line)
      next
    }
    
    if (grepl("^References\\s*:\\s*", line)) {
      ref_text <- sub("^References\\s*:\\s*", "", line)
      
      if (!is.na(current_gene) && !is.na(current_species) && !is.na(pending_title)) {
        out[[length(out) + 1]] <- tibble(
          Gene = current_gene,
          Species = current_species,
          Title = pending_title,
          Reference = ref_text,
          Status = "available"
        )
      }
      
      pending_title <- NA_character_
      next
    }
  }
  
  if (length(out) == 0) {
    return(tibble(
      Gene = character(),
      Species = character(),
      Title = character(),
      Reference = character(),
      Status = character()
    ))
  }
  
  bind_rows(out)
}

read_childhood_dementia_genes <- function(path) {
  if (!file.exists(path)) {
    return(character())
  }
  
  df <- readxl::read_excel(path)
  
  gene_col <- pick_first_existing_col(df, c(
    "Gene", "Human Symbol", "Human_Symbol", "Gene Symbol", "Symbol"
  ))
  
  if (is.na(gene_col)) {
    gene_col <- names(df)[1]
  }
  
  genes <- clean_text(df[[gene_col]])
  genes <- genes[genes != ""]
  unique(genes)
}

read_orthology <- function(path, model_name) {
  sheet_to_read <- if ("All_threshold_results" %in% readxl::excel_sheets(path)) {
    "All_threshold_results"
  } else {
    1
  }
  
  df <- readxl::read_excel(path, sheet = sheet_to_read)
  
  
  human_symbol_col <- pick_first_existing_col(df, c(
    "Human Symbol",
    "Human gene",
    "Human Gene",
    "Human_Symbol",
    "Gene"
  ))
  
  if (is.na(human_symbol_col)) {
    stop(
      "No human-gene column was found in ",
      basename(path),
      ". Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  threshold_col <- pick_first_existing_col(df, c(
    "DIOPT threshold",
    "DIOPT Threshold",
    "Threshold"
  ))
  
  rel_col <- pick_first_existing_col(df, c(
    "Orthology relationship",
    "Fission yeast orthology state",
    "Budding yeast orthology state",
    "BuddingYeast orthology state",
    "Saccharomyces cerevisiae orthology state",
    "Orthology state"
  ))
  
  diopt_col <- pick_first_existing_col(df, c("DIOPT SCORE", "DIOPT Score"))
  orth_score_col <- pick_first_existing_col(df, c("Orthology Score"))
  
  model_n_col <- pick_first_existing_col(df, c(
    "n_mice",
    "n_zebrafish",
    "n_flies",
    "n_worms",
    "n_yeasts",
    "n_budding_yeasts",
    "n_buddingyeasts",
    "n_buddingyeast",
    "n_saccharomyces",
    "n_model",
    "n_models",
    "n_model_genes"
  ))
  
  human_component_col <- pick_first_existing_col(df, c(
    "Human(s) in component",
    "Humans in component",
    "Human in component"
  ))
  
  model_component_col <- pick_first_existing_col(df, c(
    paste0(model_name, "(s) in component"),
    paste0(model_name, "s in component"),
    paste0(model_name, " in component"),
    "Mice in component",
    "Zebrafish in component",
    "Zebrafish(s) in component",
    "Flies in component",
    "Worms in component",
    "Worm(s) in component",
    "Fission yeast(s) in component",
    "Fission yeasts in component",
    "Budding yeast(s) in component",
    "Budding yeasts in component",
    "Budding Yeasts in component",
    "BuddingYeast genes in component",
    "BuddingYeasts in component",
    "Saccharomyces cerevisiae in component",
    "Saccharomyces cerevisiae gene(s) in component",
    "Model(s) in component"
  ))
  
  message(
    "  ", model_name,
    ": human column = ", human_symbol_col,
    "; relationship column = ", ifelse(is.na(rel_col), "<not found>", rel_col),
    "; model component column = ", ifelse(is.na(model_component_col), "<not found>", model_component_col),
    "; model count column = ", ifelse(
      is.na(model_n_col),
      "<not found; relationship fallback used>",
      model_n_col
    )
  )
  
  threshold_values <- if (is.na(threshold_col)) {
    rep(7, nrow(df))
  } else {
    suppressWarnings(as.numeric(df[[threshold_col]]))
  }
  
  out <- tibble(
    DIOPT_Threshold = threshold_values,
    Human_Symbol = clean_text(safe_col(df, human_symbol_col)),
    Orthology_Relationship = clean_text(safe_col(df, rel_col)),
    DIOPT_Score = clean_text(safe_col(df, diopt_col)),
    Orthology_Score = suppressWarnings(as.numeric(safe_col(df, orth_score_col, NA))),
    Human_Component = clean_text(safe_col(df, human_component_col)),
    Model_Component = clean_text(safe_col(df, model_component_col)),
    Model_Gene_N = suppressWarnings(as.numeric(safe_col(df, model_n_col, NA))),
    Model = model_name
  ) %>%
    filter(!is.na(DIOPT_Threshold)) %>%
    mutate(
      Has_Ortholog = if (!is.na(model_n_col)) {
        is_has_ortholog(Orthology_Relationship, Model_Gene_N)
      } else {
        is_has_ortholog(Orthology_Relationship)
      },
      Priority = ifelse(
        Has_Ortholog,
        priority_from_relationship(Orthology_Relationship),
        "No ortholog"
      ),
      PriorityScore = priority_score(Priority)
    )
  
  out
}

read_model_links <- function(path, model_name) {
  if (!file.exists(path)) {
    return(tibble(
      Human_Symbol = character(),
      Model_DB_Link = character(),
      Model = character()
    ))
  }
  
  df <- readxl::read_excel(path)
  
  human_col <- pick_first_existing_col(df, c("Human Symbol"))
  link_col <- pick_first_existing_col(df, c(
    "FlyBase Link",
    "MGI Link",
    "WormBase Link",
    "ZFIN Link",
    "PomBase Link",
    "SGD Link",
    "Saccharomyces Genome Database Link",
    "Budding Yeast Link"
  ))
  
  tibble(
    Human_Symbol = clean_text(safe_col(df, human_col)),
    Model_DB_Link = clean_text(safe_col(df, link_col)),
    Model = model_name
  ) %>%
    filter(Human_Symbol != "")
}

parse_diopt_max <- function(x) {
  x <- clean_text(x)
  
  denominator <- suppressWarnings(as.numeric(stringr::str_extract(x, "(?<=/)\\d+")))
  plain_number <- suppressWarnings(as.numeric(stringr::str_extract(x, "^\\d+$")))
  
  dplyr::coalesce(denominator, plain_number)
}

split_component_genes <- function(x) {
  x <- clean_text(x)
  if (x == "") return(character(0))
  
  genes <- unlist(strsplit(x, "\\s*;\\s*|\\s*,\\s*"))
  genes <- clean_text(genes)
  genes[genes != ""]
}


read_gene_reference <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Gene_reference_table.xlsx was not found: ",
      path
    )
  }
  
  available_sheets <- readxl::excel_sheets(path)
  
  sheet_to_read <- if ("Gene_Reference" %in% available_sheets) {
    "Gene_Reference"
  } else {
    available_sheets[1]
  }
  
  df <- readxl::read_excel(path, sheet = sheet_to_read)
  
  input_gene_col <- pick_first_existing_col(df, c(
    "Gene Name",
    "Input Gene Symbol",
    "Input_Gene_Symbol",
    "Original Gene",
    "Original_Gene"
  ))
  
  human_symbol_col <- pick_first_existing_col(df, c(
    "Human Symbol",
    "Human_Symbol",
    "Search Term",
    "Updated Human Symbol",
    "Current Human Symbol"
  ))
  
  if (is.na(input_gene_col)) {
    stop(
      "No Gene Name column was found in ",
      basename(path),
      ". Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  if (is.na(human_symbol_col)) {
    stop(
      "No Human Symbol column was found in ",
      basename(path),
      ". Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  out <- df %>%
    transmute(
      Input_Gene_Symbol = clean_text(.data[[input_gene_col]]),
      Human_Symbol_Reference = clean_text(.data[[human_symbol_col]]),
      Input_Gene_Key = stringr::str_to_upper(Input_Gene_Symbol),
      Human_Symbol_Key = stringr::str_to_upper(Human_Symbol_Reference)
    ) %>%
    filter(
      Input_Gene_Symbol != "",
      Human_Symbol_Reference != ""
    ) %>%
    distinct()
  
  input_gene_conflicts <- out %>%
    distinct(Input_Gene_Key, Human_Symbol_Key) %>%
    count(Input_Gene_Key, name = "n_human_symbols") %>%
    filter(n_human_symbols > 1)
  
  human_symbol_conflicts <- out %>%
    distinct(Human_Symbol_Key, Input_Gene_Key) %>%
    count(Human_Symbol_Key, name = "n_input_symbols") %>%
    filter(n_input_symbols > 1)
  
  if (nrow(input_gene_conflicts) > 0) {
    stop(
      "Gene_reference_table.xlsx contains one or more Gene Name ",
      "values that map to multiple Human Symbol values: ",
      paste(input_gene_conflicts$Input_Gene_Key, collapse = ", ")
    )
  }
  
  if (nrow(human_symbol_conflicts) > 0) {
    stop(
      "Gene_reference_table.xlsx contains one or more Human Symbol ",
      "values that map to multiple Gene Name values: ",
      paste(human_symbol_conflicts$Human_Symbol_Key, collapse = ", ")
    )
  }
  
  out %>%
    distinct(Input_Gene_Key, .keep_all = TRUE) %>%
    arrange(Input_Gene_Symbol)
}


map_to_reference_human_symbol <- function(input_symbol, gene_reference) {
  input_symbol2 <- clean_text(input_symbol)
  
  lookup <- stats::setNames(
    gene_reference$Human_Symbol_Reference,
    gene_reference$Input_Gene_Key
  )
  
  mapped <- unname(
    lookup[stringr::str_to_upper(input_symbol2)]
  )
  
  dplyr::if_else(
    !is.na(mapped) & mapped != "",
    mapped,
    input_symbol2
  )
}


map_component_to_reference_human_symbol <- function(
    component,
    gene_reference
) {
  component2 <- clean_text(component)
  
  vapply(
    component2,
    function(value) {
      if (value == "") {
        return("")
      }
      
      genes <- split_component_genes(value)
      
      if (length(genes) == 0) {
        return(value)
      }
      
      mapped_genes <- map_to_reference_human_symbol(
        genes,
        gene_reference
      )
      
      paste(mapped_genes, collapse = "; ")
    },
    character(1)
  )
}


# -------------------------
# 3) READ GENE REFERENCE TABLE
# -------------------------

message("Reading gene reference table...")

gene_reference <- read_gene_reference(file_gene_reference)

message(
  "Gene reference mappings loaded: ",
  nrow(gene_reference)
)


# -------------------------
# 4) READ DISEASE FILES
# -------------------------

message("Reading disease files...")

imd_info <- readxl::read_excel(file_info)
imd_other <- readxl::read_excel(file_other)

imd_info2 <- imd_info %>%
  mutate(
    Disease_Category = clean_text(`Disease Category`),
    Disease_Subcategory = clean_text(`Disease Subcategory`),
    ICIEM_Name = clean_text(`ICIEM Name`),
    Gene = clean_text(Gene),
    Alternative_Names = clean_text(`Alternative Names`),
    Disease_Abbreviation = clean_text(`Disease Abbreviation`),
    OMIM = clean_text(OMIM),
    OMIM_Link = clean_text(`OMIM Link`),
    OrphaCode = clean_text(OrphaCode),
    OrphaCode_Link = clean_text(`OrphaCode Link`),
    IEMBASE_Link = clean_text(`IEMBASE Link (URL)`),
    GeneReviews = clean_text(GeneReviews),
    GeneReviews_Link = clean_text(`GeneReviews Link`),
    Mode_of_Inheritance = clean_text(`Mode of Inheritance`),
    Treatable = clean_text(`Treatable?`),
    Prevalence = clean_text(Prevalence),
    Source_Table = "IMDs_Info_With_Category"
  ) %>%
  dplyr::select(
    Disease_Category, Disease_Subcategory, ICIEM_Name, Gene,
    Alternative_Names, Disease_Abbreviation,
    OMIM, OMIM_Link, OrphaCode, OrphaCode_Link,
    IEMBASE_Link, GeneReviews, GeneReviews_Link,
    Mode_of_Inheritance, Treatable, Prevalence,
    everything()
  )

imd_other2 <- imd_other %>%
  mutate(
    Disease_Category = "Other",
    Disease_Subcategory = "Other",
    ICIEM_Name = clean_text(`ICIEM Name`),
    Gene = clean_text(Gene),
    Alternative_Names = clean_text(`Alternative Names`),
    Disease_Abbreviation = clean_text(`Disease Abbreviation`),
    OMIM = clean_text(OMIM),
    OMIM_Link = clean_text(`OMIM Link`),
    OrphaCode = clean_text(OrphaCode),
    OrphaCode_Link = clean_text(`OrphaCode Link`),
    IEMBASE_Link = clean_text(`IEMBASE Link`),
    GeneReviews = "",
    GeneReviews_Link = clean_text(`GeneReviews Link`),
    Mode_of_Inheritance = clean_text(`Mode of Inheritance`),
    Treatable = clean_text(`Treatable?`),
    Prevalence = clean_text(Prevalence),
    Source_Table = "IMDs_Others"
  ) %>%
  dplyr::select(
    Disease_Category, Disease_Subcategory, ICIEM_Name, Gene,
    Alternative_Names, Disease_Abbreviation,
    OMIM, OMIM_Link, OrphaCode, OrphaCode_Link,
    IEMBASE_Link, GeneReviews, GeneReviews_Link,
    Mode_of_Inheritance, Treatable, Prevalence,
    everything()
  )

imd_all <- bind_rows(imd_info2, imd_other2) %>%
  mutate(
    Gene = clean_text(Gene),
    Gene_Key = stringr::str_to_upper(Gene)
  ) %>%
  left_join(
    gene_reference %>%
      select(
        Human_Symbol_Key,
        Input_Gene_Symbol
      ),
    by = c("Gene_Key" = "Human_Symbol_Key")
  ) %>%
  mutate(
    Input_Gene_Symbol = dplyr::coalesce(
      Input_Gene_Symbol,
      ""
    ),
    Gene_Reference_Applied = Input_Gene_Symbol != "",
    record_id = row_number(),
    search_blob = str_to_lower(paste(
      Gene,
      Input_Gene_Symbol,
      ICIEM_Name,
      Alternative_Names,
      Disease_Abbreviation,
      sep = " | "
    ))
  ) %>%
  select(-Gene_Key)

message(
  "IMD records linked to an alternative input symbol: ",
  sum(imd_all$Gene_Reference_Applied, na.rm = TRUE)
)

# -------------------------
# 5) READ KEGG PATHWAY FILE
# -------------------------

message("Reading KEGG pathway file...")

kegg_raw <- readxl::read_excel(file_kegg)

kegg_tbl <- kegg_raw %>%
  mutate(
    Category = clean_text(Category),
    Subcategory = clean_text(Subcategory),
    KEGG_Pathway = clean_text(`KEGG Pathway`),
    Pathway_KEGG_ID = clean_text(Pathway_KEGG_ID),
    Human_Gene_Symbol = clean_text(Human_Gene_Symbol)
  ) %>%
  filter(
    Category != "",
    Subcategory != "",
    KEGG_Pathway != "",
    Pathway_KEGG_ID != "",
    Human_Gene_Symbol != ""
  )

imd_gene_lookup <- imd_all %>%
  filter(Gene != "") %>%
  distinct(Gene, record_id, ICIEM_Name, Disease_Category, Disease_Subcategory)

kegg_imd_tbl <- kegg_tbl %>%
  inner_join(
    imd_gene_lookup,
    by = c("Human_Gene_Symbol" = "Gene"),
    relationship = "many-to-many"
  )


# -------------------------
# 6) READ PATHWAY CLICKMAP
# -------------------------

message("Reading pathway clickmap...")

pathway_gene_clickmap <- read.csv(
  file_gene_clickmap,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    Model = clean_text(Model),
    Model = dplyr::case_when(
      stringr::str_to_lower(Model) %in% c(
        "budding yeast",
        "buddingyeast",
        "saccharomyces cerevisiae",
        "s. cerevisiae"
      ) ~ "Saccharomyces cerevisiae",
      TRUE ~ Model
    ),
    Pathway_KEGG_ID = clean_text(Pathway_KEGG_ID),
    Gene = clean_text(Gene),
    x1 = suppressWarnings(as.numeric(x1)),
    y1 = suppressWarnings(as.numeric(y1)),
    x2 = suppressWarnings(as.numeric(x2)),
    y2 = suppressWarnings(as.numeric(y2)),
    image_width = suppressWarnings(as.numeric(image_width)),
    image_height = suppressWarnings(as.numeric(image_height))
  ) %>%
  filter(
    Model != "",
    Pathway_KEGG_ID != "",
    Gene != "",
    !is.na(x1), !is.na(y1), !is.na(x2), !is.na(y2),
    !is.na(image_width), !is.na(image_height)
  ) %>%
  mutate(
    Box_ID = paste(Pathway_KEGG_ID, x1, y1, x2, y2, sep = "__")
  ) %>%
  distinct()

# The coordinate CSV uses Mus musculus as the common reference click-map.
# The same box coordinates are reused for all model-organism pathway figures
# in app.R whenever species-specific click-map rows are unavailable.
clickmap_reference_model <- "Mus musculus"

if (!clickmap_reference_model %in% pathway_gene_clickmap$Model) {
  stop(
    paste0(
      "The reference click-map rows ('", clickmap_reference_model,
      "') were not found in ", basename(file_gene_clickmap), ". ",
      "At least one complete reference model is required."
    )
  )
}

message(
  "Pathway click-map loaded using ", clickmap_reference_model,
  " as the shared coordinate reference for all model organisms, ",
  "including Saccharomyces cerevisiae."
)

saveRDS(
  pathway_gene_clickmap,
  file = file.path(base_dir, "pathway_gene_clickmap.rds"),
  compress = "gzip"
)

# -------------------------
# 7) READ RESEARCH DOCX
# -------------------------

message("Reading research Word document...")

research_models_tbl <- parse_research_docx(file_research_doc)


# -------------------------
# 8) READ CHILDHOOD DEMENTIA GENE LIST
# -------------------------

message("Reading childhood dementia gene list...")

childhood_dementia_genes <- read_childhood_dementia_genes(file_childhood_dementia)


# -------------------------
# 9) READ ORTHOLOGY TABLES
# -------------------------

message("Reading orthology files...")

orth_mouse <- read_orthology(file_mouse, "Mouse")
orth_zfish <- read_orthology(file_zfish, "Zebrafish")
orth_fly <- read_orthology(file_fly, "Fly")
orth_worm <- read_orthology(file_worm, "Worm")
orth_yeast <- read_orthology(file_yeast, "Fission Yeast")
orth_budding_yeast <- read_orthology(file_budding_yeast, "Budding Yeast")

orth_all_raw <- bind_rows(
  orth_mouse,
  orth_zfish,
  orth_fly,
  orth_worm,
  orth_yeast,
  orth_budding_yeast
) %>%
  mutate(
    Input_Human_Symbol = clean_text(Human_Symbol),
    Human_Component_Input = clean_text(Human_Component),
    Human_Symbol = map_to_reference_human_symbol(
      Input_Human_Symbol,
      gene_reference
    ),
    Human_Component = map_component_to_reference_human_symbol(
      Human_Component_Input,
      gene_reference
    ),
    Gene_Reference_Applied =
      Human_Symbol != Input_Human_Symbol
  )

orthology_symbol_updates <- orth_all_raw %>%
  filter(Gene_Reference_Applied) %>%
  distinct(
    Model,
    Input_Human_Symbol,
    Human_Symbol
  ) %>%
  arrange(
    Model,
    Human_Symbol
  )

message(
  "Orthology records updated using Gene_reference_table.xlsx: ",
  sum(orth_all_raw$Gene_Reference_Applied, na.rm = TRUE)
)

message(
  "Unique model and gene-symbol mappings updated: ",
  nrow(orthology_symbol_updates)
)

orth_by_threshold <- split(
  orth_all_raw,
  as.character(orth_all_raw$DIOPT_Threshold)
)

available_diopt_thresholds <- sort(unique(orth_all_raw$DIOPT_Threshold))

diopt_default <- if (7 %in% available_diopt_thresholds) {
  7
} else {
  available_diopt_thresholds[which.min(abs(available_diopt_thresholds - 7))]
}


# -------------------------
# 10) DIOPT MAX MARKERS AND GENE-LEVEL DIOPT SCORES
# -------------------------

message("Calculating DIOPT score summaries...")

diopt_species_max <- orth_all_raw %>%
  mutate(
    DIOPT_Max_From_Score = parse_diopt_max(DIOPT_Score)
  ) %>%
  group_by(Model) %>%
  summarise(
    Max_DIOPT = max(DIOPT_Max_From_Score, DIOPT_Threshold, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Display_Model = recode(
      Model,
      "Mouse" = "Mouse",
      "Zebrafish" = "Zebrafish",
      "Fly" = "Fly",
      "Worm" = "Worm",
      "Fission Yeast" = "Fission yeast",
      "Budding Yeast" = "Budding yeast"
    )
  ) %>%
  group_by(Max_DIOPT) %>%
  summarise(
    Display_Model = paste(Display_Model, collapse = " / "),
    .groups = "drop"
  ) %>%
  mutate(
    Position_Pct = 100 * (Max_DIOPT - min(available_diopt_thresholds)) /
      (max(available_diopt_thresholds) - min(available_diopt_thresholds))
  )

model_gene_diopt_scores <- orth_all_raw %>%
  mutate(
    Human_Symbol = clean_text(Human_Symbol),
    Model = clean_text(Model),
    Model_Component = clean_text(Model_Component)
  ) %>%
  filter(
    Human_Symbol != "",
    Model != "",
    Model_Component != "",
    Has_Ortholog
  ) %>%
  mutate(
    Model_Gene = strsplit(Model_Component, "\\s*;\\s*|\\s*,\\s*")
  ) %>%
  tidyr::unnest(Model_Gene) %>%
  mutate(
    Model_Gene = clean_text(Model_Gene),
    Model_Gene_Key = stringr::str_to_lower(Model_Gene)
  ) %>%
  filter(Model_Gene != "") %>%
  group_by(Human_Symbol, Model, Model_Gene_Key) %>%
  summarise(
    Model_Gene_Display = first(Model_Gene),
    Max_DIOPT_Threshold = max(DIOPT_Threshold, na.rm = TRUE),
    .groups = "drop"
  )

human_gene_diopt_scores <- orth_all_raw %>%
  mutate(
    Human_Symbol = clean_text(Human_Symbol),
    Model = clean_text(Model),
    Human_Component = clean_text(Human_Component)
  ) %>%
  filter(
    Human_Symbol != "",
    Model != "",
    Human_Component != "",
    Has_Ortholog
  ) %>%
  mutate(
    Human_Component_Gene = strsplit(Human_Component, "\\s*;\\s*|\\s*,\\s*")
  ) %>%
  tidyr::unnest(Human_Component_Gene) %>%
  mutate(
    Human_Component_Gene = clean_text(Human_Component_Gene),
    Human_Component_Gene_Key = stringr::str_to_upper(Human_Component_Gene)
  ) %>%
  filter(Human_Component_Gene != "") %>%
  group_by(Human_Symbol, Model, Human_Component_Gene_Key) %>%
  summarise(
    Human_Component_Gene_Display = first(Human_Component_Gene),
    Max_DIOPT_Threshold = max(DIOPT_Threshold, na.rm = TRUE),
    .groups = "drop"
  )


# -------------------------
# 11) READ MODEL DATABASE LINKS
# -------------------------

message("Reading model database links...")

model_link_mouse <- read_model_links(file_model_mouse, "Mouse")
model_link_zfish <- read_model_links(file_model_zfish, "Zebrafish")
model_link_fly <- read_model_links(file_model_fly, "Fly")
model_link_worm <- read_model_links(file_model_worm, "Worm")
model_link_yeast <- read_model_links(file_model_yeast, "Fission Yeast")
model_link_budding_yeast <- read_model_links(
  file_model_budding_yeast,
  "Budding Yeast"
)

model_links_all <- bind_rows(
  model_link_mouse,
  model_link_zfish,
  model_link_fly,
  model_link_worm,
  model_link_yeast,
  model_link_budding_yeast
) %>%
  mutate(
    Model_Link_Input_Human_Symbol = clean_text(Human_Symbol),
    Human_Symbol = map_to_reference_human_symbol(
      Model_Link_Input_Human_Symbol,
      gene_reference
    )
  ) %>%
  distinct(Model, Human_Symbol, .keep_all = TRUE)


# -------------------------
# 12) MODEL SELECTION MATRIX
# -------------------------

message("Reading model selection matrix...")

default_model_weights <- tibble(
  Model = c(
    "Budding Yeast",
    "Fission Yeast",
    "Worm",
    "Fly",
    "Zebrafish",
    "Mouse"
  ),
  Ethic = c(3, 3, 3, 3, 2, 1),
  Generation_time = c(3, 3, 2, 2, 1, 1),
  Throughput = c(3, 3, 3, 3, 2, 1),
  IMDs_Coverage = c(1, 1, 2, 2, 3, 3),
  Frozen_Revived = c(1, 1, 1, 0, 0, 0)
) %>%
  mutate(
    Weight = Ethic + Generation_time + Throughput + IMDs_Coverage + Frozen_Revived
  )

model_weights <- if (file.exists(file_matrix)) {
  x <- readxl::read_excel(file_matrix)
  
  required_cols <- c(
    "Model", "Ethic", "Generation_time", "Throughput",
    "IMDs_Coverage", "Frozen_Revived"
  )
  
  if (all(required_cols %in% names(x))) {
    x %>%
      mutate(
        Model = clean_text(Model),
        Ethic = suppressWarnings(as.numeric(Ethic)),
        Generation_time = suppressWarnings(as.numeric(Generation_time)),
        Throughput = suppressWarnings(as.numeric(Throughput)),
        IMDs_Coverage = suppressWarnings(as.numeric(IMDs_Coverage)),
        Frozen_Revived = suppressWarnings(as.numeric(Frozen_Revived))
      ) %>%
      mutate(
        Weight = Ethic + Generation_time + Throughput + IMDs_Coverage + Frozen_Revived
      ) %>%
      filter(Model %in% default_model_weights$Model)
  } else if (all(c("Model", "Weight") %in% names(x))) {
    x %>%
      mutate(
        Model = clean_text(Model),
        Weight = suppressWarnings(as.numeric(Weight))
      ) %>%
      filter(Model %in% default_model_weights$Model)
  } else {
    default_model_weights
  }
} else {
  default_model_weights
}

# If model_selection_matrix.xlsx has not yet been updated with Budding Yeast,
# append any missing model(s) from the defaults so all six species are retained.
missing_weight_models <- setdiff(
  default_model_weights$Model,
  model_weights$Model
)

if (length(missing_weight_models) > 0) {
  model_weights <- bind_rows(
    model_weights,
    default_model_weights %>%
      filter(Model %in% missing_weight_models)
  )
}

model_weights <- model_weights %>%
  mutate(
    model_order = match(Model, default_model_weights$Model)
  ) %>%
  arrange(model_order) %>%
  dplyr::select(-model_order)


# -------------------------
# 13) CATEGORY COVERAGE BASE TABLES
# -------------------------

message("Preparing category tables...")

all_models <- c(
  "Mouse",
  "Zebrafish",
  "Fly",
  "Worm",
  "Budding Yeast",
  "Fission Yeast"
)

category_gene_table <- imd_all %>%
  mutate(
    Disease_Category = clean_text(Disease_Category),
    Gene = clean_text(Gene)
  ) %>%
  filter(Gene != "") %>%
  distinct(Disease_Category, Gene)

all_categories <- category_gene_table %>%
  distinct(Disease_Category) %>%
  mutate(
    category_order = suppressWarnings(as.numeric(stringr::str_extract(Disease_Category, "^\\d+"))),
    category_order = ifelse(is.na(category_order), 9999, category_order)
  ) %>%
  arrange(category_order, Disease_Category) %>%
  pull(Disease_Category)

total_imd_gene_n <- imd_all %>%
  filter(Gene != "") %>%
  distinct(Gene) %>%
  nrow()

model_display_order <- c(
  "Mus musculus",
  "Danio rerio",
  "Drosophila melanogaster",
  "Caenorhabditis elegans",
  "Saccharomyces cerevisiae",
  "Schizosaccharomyces pombe"
)

model_display_map <- c(
  "Mouse" = "Mus musculus",
  "Zebrafish" = "Danio rerio",
  "Fly" = "Drosophila melanogaster",
  "Worm" = "Caenorhabditis elegans",
  "Fission Yeast" = "Schizosaccharomyces pombe",
  "Budding Yeast" = "Saccharomyces cerevisiae"
)


# -------------------------
# 14) SAVE EVERYTHING
# -------------------------

message("Saving app_data.rds...")

app_data <- list(
  kegg_raw = kegg_raw,
  kegg_tbl = kegg_tbl,
  imd_all = imd_all,
  imd_gene_lookup = imd_gene_lookup,
  kegg_imd_tbl = kegg_imd_tbl,
  
  gene_reference = gene_reference,
  orthology_symbol_updates = orthology_symbol_updates,
  
  childhood_dementia_genes = childhood_dementia_genes,
  research_models_tbl = research_models_tbl,
  
  orth_all_raw = orth_all_raw,
  orth_by_threshold = orth_by_threshold,
  available_diopt_thresholds = available_diopt_thresholds,
  diopt_default = diopt_default,
  diopt_species_max = diopt_species_max,
  model_gene_diopt_scores = model_gene_diopt_scores,
  human_gene_diopt_scores = human_gene_diopt_scores,
  
  model_links_all = model_links_all,
  model_weights = model_weights,
  default_model_weights = default_model_weights,
  
  all_models = all_models,
  category_gene_table = category_gene_table,
  all_categories = all_categories,
  total_imd_gene_n = total_imd_gene_n,
  model_display_order = model_display_order,
  model_display_map = model_display_map
)

saveRDS(
  app_data,
  file = file.path(base_dir, "app_data.rds"),
  compress = "gzip"
)

message("Done. app_data.rds created at: ", file.path(base_dir, "app_data.rds"))