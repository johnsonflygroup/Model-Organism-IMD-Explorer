global_start_time <- Sys.time()
message("[TIMER] app.R sourcing started")
#############################
suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(DT)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(htmltools)
  library(scales)
  library(officer)
  library(tibble)
})
#################################
message(
  "[TIMER] packages loaded in ",
  round(as.numeric(difftime(Sys.time(), global_start_time, units = "secs")), 3),
  " sec"
)

# -------------------------
# 1) PATHS
# -------------------------

base_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

addResourcePath("assets", file.path(base_dir, "www"))
# Expose the app folder only for static images that are kept outside www.
# For deployment, the cleaner option is still to put these image files in www/.
addResourcePath("root_assets", base_dir)

find_static_asset <- function(stem) {
  exts <- c("png", "jpg", "jpeg", "webp", "gif")
  candidates <- unlist(lapply(exts, function(ext) {
    c(
      file.path(base_dir, "www", paste0(stem, ".", ext)),
      file.path(base_dir, paste0(stem, ".", ext))
    )
  }))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return("")
  hit <- normalizePath(hit[1], winslash = "/", mustWork = FALSE)
  fname <- basename(hit)
  if (grepl("/www/", hit, fixed = TRUE)) {
    paste0("assets/", utils::URLencode(fname, reserved = TRUE))
  } else {
    paste0("root_assets/", utils::URLencode(fname, reserved = TRUE))
  }
}

tool_background_src <- find_static_asset("Tool_Background")
ortholog_classification_src <- find_static_asset("Ortholog_classification")

pathway_view_dir <- file.path(base_dir, "Pathway_View")
addResourcePath("pathway_view", pathway_view_dir)

# -------------------------
# 1.2) LOAD PREPROCESSED DATA
# -------------------------

app_data_file <- file.path(base_dir, "app_data.rds")

if (!file.exists(app_data_file)) {
  stop(
    "app_data.rds was not found. ",
    "Please run prepare_data.R locally before launching or deploying the app."
  )
}

t0 <- Sys.time()
message("[TIMER] reading app_data.rds...")

app_data <- readRDS(app_data_file)

message(
  "[TIMER] app_data.rds loaded in ",
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3),
  " sec"
)

kegg_raw <- app_data$kegg_raw
kegg_tbl <- app_data$kegg_tbl
imd_all <- app_data$imd_all
imd_gene_lookup <- app_data$imd_gene_lookup
kegg_imd_tbl <- app_data$kegg_imd_tbl

childhood_dementia_genes <- app_data$childhood_dementia_genes
research_models_tbl <- app_data$research_models_tbl


orth_all_raw <- app_data$orth_all_raw
orth_by_threshold <- app_data$orth_by_threshold
available_diopt_thresholds <- app_data$available_diopt_thresholds
diopt_default <- app_data$diopt_default
diopt_species_max <- app_data$diopt_species_max
model_gene_diopt_scores <- app_data$model_gene_diopt_scores
human_gene_diopt_scores <- app_data$human_gene_diopt_scores

model_links_all <- app_data$model_links_all
model_weights <- app_data$model_weights
default_model_weights <- app_data$default_model_weights

all_models <- c(
  "Mouse",
  "Zebrafish",
  "Fly",
  "Worm",
  "Budding Yeast",
  "Fission Yeast"
)
category_gene_table <- app_data$category_gene_table
all_categories <- app_data$all_categories
total_imd_gene_n <- app_data$total_imd_gene_n
model_display_order <- app_data$model_display_order
model_display_map <- app_data$model_display_map

rm(app_data)
message(
  "[TIMER] global startup finished in ",
  round(as.numeric(difftime(Sys.time(), global_start_time, units = "secs")), 3),
  " sec"
)

pathway_gene_clickmap_cache <- NULL

get_pathway_gene_clickmap <- function() {
  if (is.null(pathway_gene_clickmap_cache)) {
    pathway_clickmap_file <- file.path(base_dir, "pathway_gene_clickmap.rds")
    
    if (!file.exists(pathway_clickmap_file)) {
      stop(
        "pathway_gene_clickmap.rds was not found. ",
        "Please run prepare_data.R before launching or deploying the app."
      )
    }
    
    t0 <- Sys.time()
    message("[TIMER] loading pathway_gene_clickmap.rds...")
    
    pathway_gene_clickmap_cache <<- readRDS(pathway_clickmap_file)
    
    message(
      "[TIMER] pathway_gene_clickmap.rds loaded in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3),
      " sec"
    )
  }
  
  pathway_gene_clickmap_cache
}
# -------------------------
# 2) HELPERS
# -------------------------
clean_text <- function(x) {
  x <- as.character(x)
  x <- ifelse(is.na(x), "", x)
  str_squish(x)
}

pick_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_col <- function(df, col, default = "") {
  if (is.na(col) || !(col %in% names(df))) return(rep(default, nrow(df)))
  df[[col]]
}

coverage_colour <- function(pct) {
  # red -> white -> blue
  pal <- col_numeric(
    palette = c("#d73027", "#fee8c8", "#91bfdb", "#4575b4"),
    domain = c(0, 100),
    na.color = "#d9d9d9"
  )
  pal(pct)
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
  
  case_when(
    rel2 %in% c("", "1:0", "0:1", "0:0", "many:0", "0:many") |
      stringr::str_detect(rel2, "no.*ortholog") ~ "No ortholog",
    rel2 == "1:1" ~ "Prioritized",
    rel2 == "1:2" ~ "Retained",
    TRUE ~ "Deprioritized"
  )
}

priority_score <- function(priority_label) {
  case_when(
    priority_label == "Prioritized" ~ 2,
    priority_label == "Retained" ~ 1,
    priority_label == "Deprioritized" ~ 0,
    TRUE ~ -1
  )
}

get_box_genes <- function(model_name, pathway_id, box_id) {
  clickmap_tbl <- get_pathway_gene_clickmap()
  
  x <- clickmap_tbl %>%
    filter(
      Model == model_name,
      Pathway_KEGG_ID == pathway_id,
      Box_ID == box_id
    )
  
  if (nrow(x) == 0) {
    x <- clickmap_tbl %>%
      filter(
        Model == "Mus musculus",
        Pathway_KEGG_ID == pathway_id,
        Box_ID == box_id
      )
  }
  
  x %>%
    distinct(Gene) %>%
    arrange(Gene) %>%
    pull(Gene)
}

get_box_disease_records <- function(model_name, pathway_id, box_id) {
  genes <- get_box_genes(model_name, pathway_id, box_id)
  
  if (length(genes) == 0) {
    return(imd_all[0, ])
  }
  
  imd_all %>%
    filter(Gene %in% genes) %>%
    distinct(record_id, .keep_all = TRUE) %>%
    arrange(Disease_Category, Disease_Subcategory, ICIEM_Name)
}

get_gene_pathways <- function(gene_symbol) {
  gene_symbol <- clean_text(gene_symbol)
  if (gene_symbol == "") return(tibble())
  
  kegg_tbl %>%
    filter(Human_Gene_Symbol == gene_symbol) %>%
    distinct(KEGG_Pathway, Pathway_KEGG_ID) %>%
    arrange(KEGG_Pathway, Pathway_KEGG_ID)
}

is_childhood_dementia_gene <- function(gene_symbol) {
  gene_symbol <- clean_text(gene_symbol)
  if (gene_symbol == "") return(FALSE)
  gene_symbol %in% childhood_dementia_genes
}

get_childhood_dementia_records <- function() {
  imd_all %>%
    dplyr::filter(Gene %in% childhood_dementia_genes) %>%
    dplyr::distinct(record_id, .keep_all = TRUE) %>%
    dplyr::arrange(Disease_Category, Disease_Subcategory, ICIEM_Name)
}

make_clickable_link <- function(url, label = NULL) {
  url <- clean_text(url)
  if (url == "") return("")
  if (is.null(label) || is.na(label) || clean_text(label) == "") label <- url
  as.character(tags$a(href = url, target = "_blank", label))
}

split_component_genes <- function(x) {
  x <- clean_text(x)
  if (x == "") return(character(0))
  
  genes <- unlist(strsplit(x, "\\s*;\\s*|\\s*,\\s*"))
  genes <- clean_text(genes)
  genes[genes != ""]
}


format_model_component_with_diopt <- function(human_symbol, model_name, model_component) {
  human_symbol <- clean_text(human_symbol)
  model_name <- clean_text(model_name)
  
  genes <- split_component_genes(model_component)
  if (length(genes) == 0) return(clean_text(model_component))
  
  out <- vapply(genes, function(g) {
    hit <- model_gene_diopt_scores %>%
      filter(
        Human_Symbol == human_symbol,
        Model == model_name,
        Model_Gene_Key == stringr::str_to_lower(clean_text(g))
      )
    
    if (nrow(hit) == 0 || is.na(hit$Max_DIOPT_Threshold[1])) {
      return(g)
    }
    
    paste0(g, " (DIOPT ", hit$Max_DIOPT_Threshold[1], ")")
  }, character(1))
  
  paste(out, collapse = "; ")
}


format_human_component_with_diopt <- function(human_symbol, model_name, human_component) {
  human_symbol <- clean_text(human_symbol)
  model_name <- clean_text(model_name)
  
  genes <- split_component_genes(human_component)
  if (length(genes) == 0) return(clean_text(human_component))
  
  out <- vapply(genes, function(g) {
    hit <- human_gene_diopt_scores %>%
      filter(
        Human_Symbol == human_symbol,
        Model == model_name,
        Human_Component_Gene_Key == stringr::str_to_upper(clean_text(g))
      )
    
    if (nrow(hit) == 0 || is.na(hit$Max_DIOPT_Threshold[1])) {
      return(g)
    }
    
    paste0(g, " (DIOPT ", hit$Max_DIOPT_Threshold[1], ")")
  }, character(1))
  
  paste(out, collapse = "; ")
}

model_display_order <- c(
  "Mus musculus",
  "Danio rerio",
  "Drosophila melanogaster",
  "Caenorhabditis elegans",
  "Saccharomyces cerevisiae",
  "Schizosaccharomyces pombe"
)

safe_id <- function(...) {
  make.names(paste(..., collapse = "__"))
}

js_escape <- function(x) {
  x <- clean_text(x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("'", "\\\\'", x)
  x <- gsub("\r|\n", " ", x)
  x
}

get_model_figure_rel <- function(model_name) {
  image_exts <- c("png", "jpg", "jpeg", "webp")
  
  # First look in Pathway_View itself, which is how the existing
  # model overview figures are normally stored.
  root_files <- list.files(pathway_view_dir, full.names = FALSE)
  root_images <- root_files[
    tolower(tools::file_ext(root_files)) %in% image_exts
  ]
  root_hit <- root_images[
    tolower(tools::file_path_sans_ext(root_images)) == tolower(model_name)
  ]
  
  if (length(root_hit) > 0) {
    return(
      paste0(
        "pathway_view/",
        utils::URLencode(root_hit[1], reserved = TRUE)
      )
    )
  }
  
  # Also allow the model overview figure to be stored inside its
  # species folder, for example:
  # Pathway_View/Saccharomyces cerevisiae/Saccharomyces cerevisiae.png
  model_folder <- file.path(pathway_view_dir, model_name)
  
  if (!dir.exists(model_folder)) {
    return(NULL)
  }
  
  folder_files <- list.files(model_folder, full.names = FALSE)
  folder_images <- folder_files[
    tolower(tools::file_ext(folder_files)) %in% image_exts
  ]
  folder_hit <- folder_images[
    tolower(tools::file_path_sans_ext(folder_images)) == tolower(model_name)
  ]
  
  if (length(folder_hit) == 0) {
    return(NULL)
  }
  
  paste0(
    "pathway_view/",
    utils::URLencode(model_name, reserved = TRUE),
    "/",
    utils::URLencode(folder_hit[1], reserved = TRUE)
  )
}

get_model_pathway_map_rel <- function(model_name, pathway_id) {
  folder <- file.path(pathway_view_dir, model_name)
  if (!dir.exists(folder)) return(NULL)
  
  files <- list.files(folder, full.names = FALSE)
  img_files <- files[tolower(tools::file_ext(files)) %in% c("png", "jpg", "jpeg", "webp")]
  hit <- img_files[startsWith(tolower(img_files), tolower(pathway_id))]
  if (length(hit) == 0) return(NULL)
  
  paste0(
    "pathway_view/",
    utils::URLencode(model_name, reserved = TRUE),
    "/",
    utils::URLencode(hit[1], reserved = TRUE)
  )
}

model_name_map_for_clickmap <- c(
  "Mus musculus" = "Mus musculus",
  "Danio rerio" = "Danio rerio",
  "Drosophila melanogaster" = "Drosophila melanogaster",
  "Caenorhabditis elegans" = "Caenorhabditis elegans",
  "Schizosaccharomyces pombe" = "Schizosaccharomyces pombe",
  "Saccharomyces cerevisiae" = "Saccharomyces cerevisiae"
)

get_clickmap_model_name <- function(display_model_name) {
  out <- unname(model_name_map_for_clickmap[display_model_name])
  if (length(out) == 0 || is.na(out)) return(display_model_name)
  out
}

get_click_genes <- function(model_name, pathway_id, click_label) {
  imd_ec_tbl %>%
    filter(
      Model == model_name,
      Pathway_KEGG_ID == pathway_id,
      Click_Label == click_label
    ) %>%
    distinct(Human_Symbol) %>%
    arrange(Human_Symbol) %>%
    pull(Human_Symbol)
}

get_click_disease_records <- function(model_name, pathway_id, click_label) {
  genes <- imd_ec_tbl %>%
    filter(
      Model == model_name,
      Pathway_KEGG_ID == pathway_id,
      Click_Label == click_label
    ) %>%
    distinct(Human_Symbol) %>%
    pull(Human_Symbol)
  
  if (length(genes) == 0) {
    return(imd_all[0, ])
  }
  
  imd_all %>%
    filter(Gene %in% genes) %>%
    distinct(record_id, .keep_all = TRUE) %>%
    arrange(Disease_Category, Disease_Subcategory, ICIEM_Name)
}

make_pubmed_link <- function(ref_text) {
  ref_text <- clean_text(ref_text)
  if (ref_text == "") return("")
  
  pmid <- stringr::str_extract(ref_text, "\\d+")
  if (is.na(pmid) || pmid == "") return(ref_text)
  
  url <- paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/")
  as.character(tags$a(href = url, target = "_blank", ref_text))
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

get_research_in_models <- function(gene_symbol) {
  gene_symbol <- clean_text(gene_symbol)
  if (gene_symbol == "") return(NULL)
  
  x <- research_models_tbl %>%
    dplyr::filter(Gene == gene_symbol) %>%
    dplyr::distinct()
  
  if (nrow(x) == 0) return(NULL)
  x
}


# -------------------------
# 7) SEARCH FUNCTION
# -------------------------
do_search <- function(query_text, data) {
  # preserve line breaks first
  query_text <- as.character(query_text)
  if (is.na(query_text) || trimws(query_text) == "") return(data[0, ])
  
  # split by line first, then clean each line
  queries <- unlist(strsplit(query_text, "\\r?\\n"))
  queries <- clean_text(queries)
  queries <- queries[queries != ""]
  
  if (length(queries) == 0) return(data[0, ])
  
  res <- purrr::map_dfr(queries, function(q) {
    q2 <- str_to_lower(q)
    
    gene_lower <- str_to_lower(clean_text(data$Gene))
    name_lower <- str_to_lower(clean_text(data$ICIEM_Name))
    alt_lower  <- str_to_lower(clean_text(data$Alternative_Names))
    abbr_lower <- str_to_lower(clean_text(data$Disease_Abbreviation))
    
    gene_dist <- as.vector(adist(q2, gene_lower))
    name_dist <- as.vector(adist(q2, name_lower))
    
    tmp <- data %>%
      mutate(
        score_gene = ifelse(str_detect(gene_lower, fixed(q2)), 4, 0),
        score_name = ifelse(str_detect(name_lower, fixed(q2)), 3, 0),
        score_alt  = ifelse(str_detect(alt_lower, fixed(q2)), 2, 0),
        score_abbr = ifelse(str_detect(abbr_lower, fixed(q2)), 2, 0),
        score_fuzzy = ifelse(gene_dist <= 1 | name_dist <= 2, 1, 0),
        score_total = score_gene + score_name + score_alt + score_abbr + score_fuzzy
      ) %>%
      filter(score_total > 0) %>%
      arrange(desc(score_total), Disease_Category, Disease_Subcategory, ICIEM_Name)
    
    tmp
  })
  
  res %>%
    distinct(record_id, .keep_all = TRUE) %>%
    arrange(desc(score_total), Disease_Category, Disease_Subcategory, ICIEM_Name)
}
# -------------------------
# 8) DISEASE DETAIL + MODEL RANKING
# -------------------------
get_disease_models <- function(disease_row, orth_tbl) {
  gene <- clean_text(disease_row$Gene[[1]])
  if (gene == "") {
    return(tibble(
      Model = all_models,
      Human_Symbol = gene,
      Orthology_Relationship = "",
      DIOPT_Score = "",
      Orthology_Score = NA_real_,
      Human_Component = "",
      Model_Component = "",
      Model_DB_Link = "",
      Has_Ortholog = FALSE,
      Priority = "No ortholog",
      PriorityScore = -1,
      Weight = model_weights$Weight[match(all_models, model_weights$Model)],
      FinalScore = NA_real_
    ))
  }
  
  if (!"Human_Component" %in% names(orth_tbl)) {
    orth_tbl$Human_Component <- ""
  }
  if (!"Model_Component" %in% names(orth_tbl)) {
    orth_tbl$Model_Component <- ""
  }
  
  x <- orth_tbl %>%
    filter(Human_Symbol == gene) %>%
    select(Model, Human_Symbol, Orthology_Relationship, DIOPT_Score,
           Orthology_Score, Human_Component, Model_Component,
           Has_Ortholog, Priority, PriorityScore) %>%
    left_join(
      model_links_all,
      by = c("Model", "Human_Symbol")
    )
  
  if (nrow(x) == 0) {
    x <- tibble(
      Model = all_models,
      Human_Symbol = gene,
      Orthology_Relationship = "",
      DIOPT_Score = "",
      Orthology_Score = NA_real_,
      Human_Component = "",
      Model_Component = "",
      Model_DB_Link = "",
      Has_Ortholog = FALSE,
      Priority = "No ortholog",
      PriorityScore = -1
    )
  } else {
    x <- tibble(Model = all_models) %>%
      left_join(x, by = "Model") %>%
      mutate(
        Human_Symbol = replace_na(Human_Symbol, gene),
        Orthology_Relationship = replace_na(Orthology_Relationship, ""),
        DIOPT_Score = replace_na(DIOPT_Score, ""),
        Human_Component = replace_na(Human_Component, ""),
        Model_Component = replace_na(Model_Component, ""),
        Model_DB_Link = replace_na(Model_DB_Link, ""),
        Has_Ortholog = replace_na(Has_Ortholog, FALSE),
        Priority = replace_na(Priority, "No ortholog"),
        PriorityScore = replace_na(PriorityScore, -1)
      )
  }
  
  n_prioritized <- sum(x$Priority == "Prioritized", na.rm = TRUE)
  
  x <- x %>%
    left_join(model_weights, by = "Model") %>%
    mutate(
      Weight = replace_na(Weight, 0),
      FinalScore = ifelse(
        n_prioritized >= 2 & Priority == "Prioritized",
        Weight,
        ifelse(Priority == "Prioritized", 100, NA_real_)
      ),
      Human_Component_With_DIOPT = purrr::pmap_chr(
        list(Human_Symbol, Model, Human_Component),
        format_human_component_with_diopt
      ),
      Model_Component_With_DIOPT = purrr::pmap_chr(
        list(Human_Symbol, Model, Model_Component),
        format_model_component_with_diopt
      )
    ) %>%
    arrange(desc(PriorityScore), desc(FinalScore), desc(Orthology_Score))
  
  x
}

make_diopt_max_markers <- function(marker_df) {
  tags$div(
    id = "diopt_max_marker_layer",
    class = "diopt-max-marker-layer",
    
    lapply(seq_len(nrow(marker_df)), function(i) {
      
      tooltip_text <- paste0(
        marker_df$Display_Model[i],
        " max DIOPT score = ",
        marker_df$Max_DIOPT[i]
      )
      
      tags$div(
        class = "diopt-max-marker",
        `data-value` = marker_df$Max_DIOPT[i],
        title = tooltip_text,
        
        tags$div(class = "diopt-max-marker-line")
      )
    })
  )
}


# -------------------------
# 9) UI
# -------------------------
ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML(paste0("
      html, body { min-height: 100%; }
      body {
        font-family: Arial, sans-serif;
        background-color: #f7fbff;
        background-image: linear-gradient(rgba(255,255,255,0.45), rgba(255,255,255,0.45)), url('", tool_background_src, "');
        background-size: cover;
        background-position: center center;
        background-attachment: fixed;
        background-repeat: no-repeat;
      }

      .container-fluid {
        max-width: 1820px;
        margin-left: auto;
        margin-right: auto;
      }
      
      /* Make result-page tabs more dominant */
      .results-page-wrap .nav-tabs {
        border-bottom: 1px solid rgba(120, 140, 160, 0.35);
        margin-top: 18px;
        margin-bottom: 22px;
        padding-left: 10px;
      }

      /* Inactive tab buttons */
      .results-page-wrap .nav-tabs > li > a {
        font-size: 20px;
        font-weight: 600;
        padding: 14px 28px;
        margin-right: 10px;

        color: #1f5f9c;
        background: rgba(255, 255, 255, 0.58);

        border: 1px solid rgba(180, 200, 220, 0.75);
        border-radius: 10px 10px 0 0;

        box-shadow: 0 2px 6px rgba(60, 90, 130, 0.08);
        transition: all 0.15s ease;
      }

      /* Hover effect */
      .results-page-wrap .nav-tabs > li > a:hover {
        color: #123b62;
        background: rgba(235, 246, 255, 0.92);
        border-color: rgba(90, 145, 190, 0.75);
        box-shadow: 0 4px 10px rgba(60, 90, 130, 0.14);
      }

      /* Active tab button */
      .results-page-wrap .nav-tabs > li.active > a,
      .results-page-wrap .nav-tabs > li.active > a:hover,
      .results-page-wrap .nav-tabs > li.active > a:focus {
        font-size: 21px;
        font-weight: 700;

        color: #1f2d3a;
        background: rgba(255, 255, 255, 0.96);

        border: 1px solid rgba(120, 145, 165, 0.75);
        border-bottom-color: rgba(255, 255, 255, 0.96);

        box-shadow: 0 -1px 0 rgba(255, 255, 255, 0.9),
                    0 4px 12px rgba(60, 90, 130, 0.12);
      }

      .app-header-title {
        font-size: 32px;
        font-weight: 700;
        color: #2f2f2f;
        margin: 18px 0 12px 0;
        line-height: 1.18;
        letter-spacing: -0.2px;
      }
      
      .landing-search-card .form-group,
      .landing-search-card .shiny-input-container {
        width: 100% !important;
        max-width: none !important;
      }
      
      .landing-search-card .section-helper-text {
        font-size: 17px;
        line-height: 1.5;
        margin-top: 14px;
      }

      .landing-search-card .help-btn {
        width: 36px;
        height: 36px;
        line-height: 31px;
        font-size: 22px;
      }
      

      .landing-search-card .irs {
        width: 100% !important;
      }
      
      .top-nav-wrap {
        width: 100%;
        border-bottom: 1px solid rgba(80, 90, 105, 0.22);
        padding: 0 0 11px 0;
        margin-bottom: 24px;
      }

      .landing-nav-bar {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        justify-content: flex-start;
        gap: 14px;
      }

      /* Smaller and more polished top navigation buttons */
      .btn.landing-nav-btn {
        min-width: 148px;
        height: 38px;

        font-size: 16px !important;
        font-weight: 500 !important;
        line-height: 1.15 !important;

        padding: 6px 18px !important;

        border-radius: 9px !important;
        border: 1px solid rgba(95, 110, 125, 0.42) !important;

        background: linear-gradient(
          180deg,
          rgba(255, 255, 255, 0.88) 0%,
          rgba(238, 244, 250, 0.76) 100%
        ) !important;

        color: #1f2d3a !important;

        box-shadow:
          0 1px 3px rgba(40, 70, 100, 0.10),
          inset 0 1px 0 rgba(255, 255, 255, 0.75);

        transition:
          background 0.15s ease,
          border-color 0.15s ease,
          box-shadow 0.15s ease,
          transform 0.15s ease;
      }

      /* Hover effect */
      .btn.landing-nav-btn:hover {
        background: linear-gradient(
          180deg,
          rgba(255, 255, 255, 0.98) 0%,
          rgba(229, 240, 250, 0.90) 100%
        ) !important;

        border-color: rgba(60, 110, 160, 0.58) !important;

        box-shadow:
          0 3px 8px rgba(40, 80, 120, 0.16),
          inset 0 1px 0 rgba(255, 255, 255, 0.85);

        transform: translateY(-1px);
      }

      /* Currently selected top mode */
      .btn.landing-nav-btn.active-mode {
       background: linear-gradient(
          180deg,
          rgba(225, 241, 255, 0.96) 0%,
          rgba(202, 224, 246, 0.88) 100%
        ) !important;

        border-color: rgba(45, 105, 165, 0.72) !important;
        color: #123b62 !important;

        box-shadow:
          0 3px 9px rgba(50, 100, 150, 0.18),
          inset 0 1px 0 rgba(255, 255, 255, 0.85);
      }

      /* Pressed/focus state */
      .btn.landing-nav-btn:focus,
      .btn.landing-nav-btn:active {
        outline: none !important;
        box-shadow:
          0 1px 4px rgba(40, 80, 120, 0.18),
          0 0 0 3px rgba(80, 140, 200, 0.14) !important;
      }

      /* The KEGG label is longer, so give only this button more width */
      .btn.landing-nav-btn-wide {
        min-width: 210px;
      }

      /* More Info can be shorter */
      .btn.landing-nav-btn-small {
        min-width: 120px;
      }

      .page-panel {
        background: rgba(255,255,255,0.70);
        border-radius: 2px;
        padding: 12px;
      }

      .landing-search-layout {
        width: min(96vw, 1680px);
        margin: 32px auto 30px auto;
        display: grid;
        grid-template-columns: minmax(850px, 1080px) minmax(350px, 430px);
        column-gap: 110px;
        align-items: start;
        justify-content: center;
      }
      
      .landing-search-card {
        width: 100%;
        max-width: none;
      }

      .landing-search-card .form-group,
      .landing-search-card .shiny-input-container {
        width: 100% !important;
        max-width: none !important;
      }

      .landing-search-card textarea.form-control {
        width: 100% !important;
        max-width: none !important;
        min-height: 190px;
        font-size: 14px;
        line-height: 1.45;
      }
      
      .landing-search-card label {
        font-size: 18px;
        line-height: 1.35;
        font-weight: 700;
        margin-bottom: 12px;
      }

      .landing-search-card textarea.form-control {
        width: 100% !important;
        max-width: none !important;
        min-height: 310px;
        font-size: 21px;
        line-height: 1.45;
        padding: 16px 18px;
      }

      .landing-search-card textarea.form-control::placeholder {
        font-size: 21px;
        line-height: 1.45;
        color: #8a8f96;
      }

      .landing-search-card .irs {
        width: 100% !important;
      }

      #landing_one_to_one_summary {
        width: 100%;
      }

      .search-panel-box {
        background: rgba(247,249,252,0.88);
        border: 1px solid #dde6f0;
        border-radius: 18px;
        padding: 34px 42px 32px 42px;
        box-shadow: 0 4px 16px rgba(60,90,130,0.12);
      }

      .search-panel-box h4 {
        font-size: 30px;
        margin-top: 0;
        margin-bottom: 24px;
        font-weight: 500;
      }

      .landing-summary-box {
        width: 100%;
        max-width: none;
        background: rgba(255,255,255,0.68);
        border: 3px solid rgba(55,55,55,0.55);
        border-radius: 34px;
        padding: 22px 28px 22px 28px;
        box-shadow: 0 3px 12px rgba(60,90,130,0.08);
      }

      .landing-summary-title {
        font-size: 18px;
        font-weight: 700;
        margin-bottom: 2px;
        color: #111;
      }

      .landing-summary-threshold {
        font-size: 13px;
        margin-bottom: 18px;
        color: #333;
      }

      .landing-summary-subtitle {
        font-size: 17px;
        font-weight: 700;
        margin-bottom: 10px;
      }

      .landing-summary-row {
        display: grid;
        grid-template-columns: 1fr auto;
        gap: 16px;
        align-items: baseline;
        font-size: 14px;
        line-height: 1.38;
      }

      .landing-summary-value {
        font-weight: 700;
        text-align: right;
        white-space: nowrap;
      }
      
      #landing_one_to_one_summary {
        width: 100%;
        margin-top: 26px;
      }

      .landing-summary-species {
        font-style: italic;
      }

      .landing-summary-value {
        font-weight: 700;
        text-align: right;
        white-space: nowrap;
      }

      .classification-page {
        max-width: 980px;
        margin: 10px auto 50px auto;
        background: rgba(255,255,255,0.68);
        padding: 10px 14px 20px 14px;
      }

      .browse-model-page {
        max-width: 1680px;
        margin: 0 auto 50px auto;
        background: rgba(255,255,255,0.42);
        padding: 18px 26px 34px 26px;
      }

      .more-info-page {
        max-width: 1600px;
        margin: 10px 0 50px 18px;
      }

      .summary-records-box {
        background: rgba(247,248,251,0.88);
        border: 1px solid #dfe5ef;
        border-radius: 14px;
        padding: 24px 30px 28px 30px;
        max-width: 780px;
        margin-bottom: 46px;
      }
      .summary-records-btn-row {
        display: flex;
        flex-direction: row;
        flex-wrap: wrap;
        gap: 34px;
        align-items: center;
        justify-content: center;
        margin-top: 16px;
        width: 100%;
      }

      .summary-records-btn {
        min-width: 230px;
        height: 52px;
        font-size: 14px;
        font-weight: 700;
        padding: 6px 16px;
        border-radius: 10px;
        border: 1px solid #cfd8e3;
        background: linear-gradient(180deg, #ffffff 0%, #eef3f8 100%);
        color: #33485c !important;
        box-shadow: 0 1px 3px rgba(60, 90, 130, 0.08);
      }
      
      .advanced-search-box {
        margin-top: 18px;
        margin-bottom: 18px;
        border: 1px solid rgba(190, 205, 220, 0.95);
        border-radius: 10px;
        background: rgba(255, 255, 255, 0.55);
        overflow: hidden;
      }

      .advanced-search-box > summary {
        list-style: none;
      }

      .advanced-search-box > summary::-webkit-details-marker {
        display: none;
      }

      .advanced-search-summary {
        cursor: pointer;
        display: flex;
        align-items: center;
        gap: 14px;
        padding: 12px 16px;
        background: linear-gradient(
          180deg,
          rgba(255, 255, 255, 0.86) 0%,
          rgba(238, 244, 250, 0.76) 100%
        );
        border-radius: 10px;
        user-select: none;
      }

      .advanced-search-summary:hover {
        background: linear-gradient(
          180deg,
          rgba(255, 255, 255, 0.96) 0%,
          rgba(229, 240, 250, 0.90) 100%
        );
      }

      .advanced-search-summary::before {
        content: '▶';
        font-size: 13px;
        color: #345c7c;
        transition: transform 0.15s ease;
      }

     .advanced-search-box[open] .advanced-search-summary::before {
        content: '▼';
      }

      .advanced-search-title {
        font-size: 17px;
        font-weight: 700;
        color: #1f3b5b;
      }

      .advanced-search-subtitle {
        font-size: 13px;
        color: #667788;
      }

      .advanced-search-content {
        padding: 14px 18px 8px 18px;
        border-top: 1px solid rgba(200, 215, 230, 0.75);
        background: rgba(255, 255, 255, 0.35);
      }
      
      .advanced-search-notes {
        margin: 18px 0 8px 0;
        padding-left: 24px;
        color: #5f6b76;
        font-size: 17px;
        line-height: 1.55;
      }

      .advanced-search-notes li {
        margin-bottom: 8px;
      }

      .advanced-search-notes strong {
        color: #33485c;
      }

      .ortholog-classification-title {
        display: inline-block;
        border: none;
        border-radius: 0;
        background: rgba(255,255,255,0.55);
        padding: 14px 34px;
        font-size: 24px;
        margin-bottom: 18px;
      }
      
      .summary-records-title {
        display: inline-block;
        font-size: 24px;
        font-weight: 400;
        color: #2f2f2f;
        margin-top: 0;
        margin-bottom: 20px;
      }

      .ortholog-classification-img {
        display: block;
        max-width: 1500px;
        width: 100%;
        height: auto;
        background: rgba(255,255,255,0.80);
        border-radius: 8px;
      }

      .classification-heading {
        display: flex;
        justify-content: flex-start;
        align-items: center;
        gap: 14px;

        margin: 0 0 24px 0;

        font-size: 22px;
        background: rgba(255,255,255,0.55);
        padding: 10px 20px;
        width: fit-content;
      }
      
      .classification-heading .help-btn {
        width: 36px;
        height: 36px;
        line-height: 31px;
        font-size: 22px;
        border-radius: 50%;
        padding: 0;
        margin-left: 8px;
        flex: 0 0 auto;
      }

      .classification-heading img {
        width: 56px;
        height: 56px;
        object-fit: contain;
      }

      .category-box {
        border: 1px solid #d9d9d9;
        border-radius: 10px;
        padding: 14px 16px;
        margin-bottom: 16px;
        background: rgba(250,250,250,0.90);
      }

      .disease-title-row {
        display: flex;
        align-items: center;
        gap: 8px;
        margin-top: 20px;
        margin-bottom: 10px;
      }

      .disease-title-row h2 { margin: 0; }

      .coverage-chip {
        display: inline-block;
        min-width: 90px;
        text-align: center;
        padding: 6px 10px;
        margin: 3px 6px 3px 0;
        border-radius: 16px;
        color: #111;
        font-weight: 600;
        font-size: 12px;
      }

      a.coverage-chip,
      a.coverage-chip:hover,
      a.coverage-chip:focus,
      a.coverage-chip:active {
        color: #111 !important;
        text-decoration: none !important;
        outline: none !important;
        box-shadow: none !important;
      }

      .disease-link {
        display: inline-block;
        margin: 4px 6px 4px 0;
        padding: 6px 10px;
        border-radius: 8px;
        background: #f0f0f0;
        text-decoration: none;
        color: #222;
      }

      .model-card {
        border: 1px solid #d9d9d9;
        border-radius: 10px;
        padding: 12px;
        margin-bottom: 10px;
        background: rgba(255,255,255,0.92);
      }

      .top-title { margin-bottom: 4px; }
      .small-muted { color: #666; font-size: 12px; }

      .golden-banner-card-link {
        display: block;
        text-decoration: none !important;
        color: inherit !important;
        cursor: pointer;
        transition: transform 0.15s ease, box-shadow 0.15s ease, border-color 0.15s ease;
      }

      .golden-banner-card-link:hover {
        text-decoration: none !important;
        color: inherit !important;
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(60, 90, 130, 0.12);
        border-color: #bcd3ea;
      }

      .golden-banner-card-link:focus,
      .golden-banner-card-link:active {
        text-decoration: none !important;
        color: inherit !important;
        outline: none;
      }

      .golden-banner-card-wrap {
        flex: 0 0 310px;
        width: 310px;
        text-align: center;
        padding-left: 0;
        padding-right: 0;
      }

      .golden-banner-secondary-btn {
        display: inline-block;
        margin-top: 10px;
        padding: 10px 16px;
        min-width: 285px;

        border-radius: 10px;
        background: #eef2f7;
        border: 1px solid #d6dee8;

        color: #43576b !important;
        text-decoration: none !important;

        font-size: 15px;
        font-weight: 700;
        line-height: 1.2;

        cursor: pointer;
        box-shadow: 0 2px 6px rgba(60, 90, 130, 0.08);
      }

      .golden-banner-secondary-value {
        font-size: 14px;
        font-weight: 500;
        color: #5f7084;
        margin-top: 4px;
        line-height: 1.25;
      }

      .golden-highlight-text {
        display: inline-block;
        font-weight: 700;
        color: #1f3b5b;
        background: #fff3bf;
        padding: 2px 8px;
        border-radius: 6px;
      }

      .golden-collapse-btn,
      .nav-collapse-btn {
        background: #eef3f8;
        border-color: #cfd8e3;
        color: #4a5a6a;
        font-weight: 600;
      }

      .category-header-row { margin-bottom: 6px; }
      .coverage-row { margin-bottom: 8px; }
      .category-details, .subcategory-details { margin-top: 8px; }

      .category-details > summary,
      .subcategory-details > summary {
        cursor: pointer;
        list-style: none;
        outline: none;
        user-select: none;
      }

      .category-details > summary::-webkit-details-marker,
      .subcategory-details > summary::-webkit-details-marker { display: none; }

      .gene-pathway-link {
        color: #337ab7;
        text-decoration: underline;
        cursor: pointer;
        font-weight: 500;
      }

      .gene-pathway-link:hover { color: #23527c; }
      /* Make search result rows look clickable */
      #search_table table.dataTable tbody tr {
        cursor: pointer;
      }

      #search_table table.dataTable tbody td {
        cursor: pointer;
      }

      #search_table table.dataTable tbody tr:hover {
        background-color: rgba(210, 230, 250, 0.45) !important;
      }

      .category-details > summary {
        font-weight: 600;
        color: #2c3e50;
      }

      .subcategory-details > summary {
        font-weight: 600;
        color: #2c3e50;
        display: flex;
        align-items: flex-start;
        justify-content: flex-start;
        gap: 8px;
        width: 100%;
      }

      .subcategory-details > summary .subcat-title {
        flex: 1;
        text-align: left;
      }

      .subcategory-details > summary .subcat-toggle-text {
        flex: 0 0 auto;
        margin-left: auto;
        text-align: right;
        white-space: nowrap;
      }

      .category-details > summary::before { content: '▶ '; font-size: 12px; }
      .category-details[open] > summary::before { content: '▼ '; }
      .subcategory-details > summary::before { content: '▶ '; font-size: 11px; }
      .subcategory-details[open] > summary::before { content: '▼ '; }
      .category-details > summary .toggle-label::after { content: ' Show more'; }
      .category-details[open] > summary .toggle-label::after { content: ' Show less'; }
      .subcategory-details > summary .subcat-toggle-text::after {
        content: 'Show more';
        font-size: 12px;
        color: #666;
        margin-left: 12px;
      }
      .subcategory-details[open] > summary .subcat-toggle-text::after { content: 'Show less'; }
      .category-content { padding-top: 10px; }
      .subcategory-content { padding-top: 8px; padding-left: 8px; }

      .section-title {
        font-size: 28px;
        font-weight: 500;
        margin: 0;
      }
      
      #search_btn {
        font-size: 18px;
        padding: 10px 24px;
        border-radius: 7px;
      }

      #try_example_btn {
        font-size: 17px !important;
        padding: 9px 22px !important;
        border-radius: 7px !important;
      }

      .pathway-header-right {
        display: flex;
        flex-direction: column;
        align-items: flex-end;
        gap: 10px;
        min-width: 420px;
      }

      .pathway-legend {
        display: block;
        width: 420px;
        max-width: 100%;
        height: auto;
      }

      .pathway-gene-highlight {
        position: absolute;
        display: block;
        border: 3px solid #ff3b30;
        background: rgba(255, 59, 48, 0.12);
        box-shadow: 0 0 0 2px rgba(255, 59, 48, 0.18);
        pointer-events: none;
        z-index: 20;
        border-radius: 2px;
      }

      .research-section-title {
        font-size: 28px;
        font-weight: 500;
        margin-top: 28px;
        margin-bottom: 28px;
      }

      .research-unavailable {
        font-size: 20px;
        color: #666;
        margin-left: 34px;
        margin-bottom: 24px;
      }

      .section-divider {
        border-top: 1.5px dashed #444;
        margin-top: 10px;
        margin-bottom: 10px;
        width: 100%;
      }

      .childhood-dementia-mark {
        display: inline-block;
        margin-left: 6px;
        color: #c0392b;
        font-weight: 700;
        cursor: help;
      }

      .help-btn {
        display: inline-block;
        width: 26px;
        height: 26px;
        line-height: 22px;
        text-align: center;
        border-radius: 50%;
        border: 1px solid #999;
        background: #f8f9fa;
        color: #333;
        font-weight: 700;
        font-size: 16px;
        cursor: pointer;
        margin-left: 8px;
        padding: 0;
      }
      
      .golden-banner .help-btn {
        width: 36px;
        height: 36px;
        line-height: 31px;
        font-size: 22px;
        border-radius: 50%;
        padding: 0;
        margin-left: 10px;
      }

      .help-btn:hover { background: #e9ecef; }
      .help-section-title { font-weight: 700; margin-top: 12px; margin-bottom: 8px; }

      .section-helper-text {
        color: #5f6b76;
        font-size: 13px;
        line-height: 1.5;
        margin-top: 2px;
        margin-bottom: 12px;
      }

      .golden-banner {
        background: linear-gradient(
          135deg,
          rgba(247,251,255,0.78) 0%,
          rgba(238,246,255,0.70) 100%
        );
        border: 1px solid #d8e6f5;
        border-radius: 18px;
        padding: 26px 34px 36px 34px;
        margin-bottom: 20px;
        box-shadow: 0 4px 16px rgba(60, 90, 130, 0.10);
      }

      .golden-banner-title {
        font-size: 30px;
        font-weight: 500;
        color: #2f2f2f;
        margin-top: 0;
        margin-bottom: 0;
        line-height: 1.25;
      }

      .golden-banner-subtitle {
        font-size: 18px;
        line-height: 1.35;
        font-weight: 700;
        color: #2f2f2f;
        margin-top: 14px;
        margin-bottom: 28px;
      }

      .golden-banner-grid {
        display: flex;
        flex-wrap: wrap;
        justify-content: center;
        align-items: flex-start;

        max-width: 1160px;
        margin: 30px auto 0 auto;

        gap: 30px 70px;
      }

      .golden-banner-card {
        width: 310px;
        min-height: 275px;

        background: rgba(255,255,255,0.86);
        border: 1px solid #d6e5f5;
        border-radius: 16px;

        padding: 22px 24px 20px 24px;
        text-align: center;

        box-shadow: 0 4px 14px rgba(60, 90, 130, 0.10);
      }

      .golden-banner-img {
        width: 145px;
        height: 95px;
        object-fit: contain;
        margin-bottom: 12px;
      }

      .golden-banner-model {
        font-size: 17px;
        font-style: italic;
        font-weight: 700;
        color: #2f3c4a;
        line-height: 1.25;
        margin-bottom: 14px;
        min-height: 34px;
      }

      .golden-banner-number {
        font-size: 34px;
        font-weight: 700;
        color: #1f3b5b;
        line-height: 1.05;
        margin-bottom: 6px;
      }


      .golden-banner-prop {
        font-size: 15px;
        color: #5f7084;
      }

      .golden-banner-support-label {
        font-size: 15px;
        font-weight: 500;
        color: #4f5f73;
        margin-bottom: 6px;
      }

      .golden-banner-secondary-title {
        font-size: 13px;
        font-weight: 700;
        line-height: 1.2;
      }

      .golden-banner-secondary-value {
        font-size: 12px;
        font-weight: 500;
        color: #5f7084;
        margin-top: 3px;
        line-height: 1.2;
      }

      .search-note {
        margin-top: 12px;
        margin-bottom: 8px;
        padding: 10px 12px;
        border-left: 4px solid #f0ad4e;
        background: #fff8e8;
        color: #5c4b1f;
        font-size: 14px;
        line-height: 1.5;
        border-radius: 4px;
      }

      .search-note-icon { font-weight: 700; margin-right: 8px; }
      .logic-note { margin-top: 8px; margin-bottom: 8px; }

      .relationship-filter-wrap { margin-bottom: 12px; }
      .relationship-filter-title { font-weight: 700; font-size: 16px; margin-bottom: 10px; }
      .relationship-checkbox-group .shiny-options-group {
        display: flex;
        flex-wrap: wrap;
        gap: 18px;
        align-items: center;
      }
      .relationship-checkbox-group .checkbox {
        margin-top: 0;
        margin-bottom: 0;
        font-weight: 600;
        font-size: 15px;
      }
      .relationship-checkbox-group input[type='checkbox'] {
        width: 18px;
        height: 18px;
        margin-right: 8px;
        vertical-align: middle;
      }

      .research-species-block { margin-bottom: 36px; padding-left: 34px; }
      .research-species-title { font-size: 22px; font-weight: 700; margin-bottom: 24px; }
      .research-line { font-size: 20px; margin-left: 34px; margin-bottom: 16px; }
      .model-name-highlight {
        font-weight: 700;
        color: #1f3b5b;
        background: #fff3bf;
        padding: 1px 6px;
        border-radius: 6px;
        display: inline-block;
      }

      .gene-chip {
        display: inline-block;
        margin: 4px 6px 4px 0;
        padding: 6px 10px;
        border-radius: 8px;
        background: #eef3f8;
        text-decoration: none;
        color: #222;
        font-weight: 500;
      }

      .pathway-map-title { font-size: 28px; font-weight: 700; margin-bottom: 20px; }
      .pathway-map-helper-title {
        font-size: 16px;
        font-weight: 600;
        color: #33485c;
        margin-top: -6px;
        margin-bottom: 16px;
        line-height: 1.5;
      }
      .pathway-map-helper-text {
        font-size: 15px;
        color: #5f6b76;
        margin-bottom: 18px;
        line-height: 1.6;
        max-width: 950px;
      }
      .pathway-model-helper-text {
        font-size: 15px;
        color: #5f6b76;
        line-height: 1.6;
        max-width: 980px;
        margin-bottom: 20px;
      }
      .pathway-model-helper-text ul { margin-top: 0; margin-bottom: 0; padding-left: 22px; }
      .pathway-model-helper-text li { margin-bottom: 8px; }
      .pathway-model-grid {
        display: flex;
        flex-wrap: wrap;
        gap: 26px 36px;
        justify-content: center;
        margin-top: 10px;
      }
      .pathway-model-card { width: 280px; text-align: center; }
      .pathway-model-thumb {
        width: 100%;
        max-width: 280px;
        height: 220px;
        object-fit: contain;
        background: #f5f5f8;
        padding: 10px;
        border-radius: 8px;
      }
      .pathway-model-name {
        display: inline-block;
        margin-top: 8px;
        font-size: 20px;
        font-style: italic;
        font-weight: 600;
        color: #222;
        text-decoration: none !important;
      }
      .pathway-model-name:hover { text-decoration: none !important; color: #222; }
      .pathway-map-image {
        width: 100%;
        max-width: 1100px;
        border: 1px solid #888;
        background: white;
        display: block;
      }
      .pathway-map-subtitle { font-size: 22px; font-style: italic; font-weight: 600; margin-bottom: 12px; }
      .pathway-map-wrap { margin-top: 10px; }
      .pathway-map-container {
        position: relative;
        display: inline-block;
        width: 100%;
        max-width: 1100px;
      }
      .pathway-ec-hotspot {
        position: absolute;
        display: block;
        background: rgba(255, 0, 0, 0.00);
        border: 0;
        cursor: pointer;
        z-index: 10;
      }
      .pathway-ec-hotspot:hover {
        background: rgba(255, 0, 0, 0.10);
        outline: 2px solid rgba(255, 0, 0, 0.35);
      }

      .diopt-slider-wrap { position: relative; }
      .diopt-max-marker-layer {
        position: relative;
        height: 24px;
        margin-top: -4px;
        margin-bottom: 8px;
      }
      .diopt-max-marker {
        position: absolute;
        top: 0;
        width: 18px;
        height: 20px;
        transform: translateX(-50%);
        cursor: help;
      }
      .diopt-max-marker-line {
        width: 2px;
        height: 16px;
        background: #c0392b;
        margin: 0 auto;
      }

      @media (max-width: 1100px) {

        .landing-search-layout {
          width: min(94vw, 900px);
          grid-template-columns: 1fr;
          row-gap: 28px;
          margin: 24px auto 30px auto;
        }

        .landing-summary-box {
          max-width: 430px;
          margin-left: auto;
          margin-right: auto;
          border-radius: 34px;
        }

        #landing_one_to_one_summary {
          max-width: 430px;
          margin-left: auto;
          margin-right: auto;
          margin-top: 20px;
        }

        .landing-nav-bar {
          gap: 10px;
        }
        
        .golden-banner-grid {
          max-width: 760px;
          gap: 28px 42px;
        }

        .golden-banner-card-wrap {
          flex: 0 0 290px;
          width: 290px;
        }

        .golden-banner-card {
          width: 290px;
          min-height: 270px;
          padding: 22px 22px 20px 22px;
        }
     
        .golden-banner-img {
          height: 92px;
          margin-bottom: 12px;
        }

        .golden-banner-secondary-btn {
          min-width: 265px;
        }

        .golden-banner-title {
          font-size: 26px;
        }

        .golden-banner-subtitle {
          font-size: 16px;
        }

        .btn.landing-nav-btn {
          min-width: 130px;
          height: 36px;
          font-size: 15px !important;
          padding: 5px 14px !important;
        }

        .btn.landing-nav-btn-wide {
          min-width: 190px;
        }

        .btn.landing-nav-btn-small {
          min-width: 110px;
        }

        .app-header-title {
          font-size: 28px;
        }
      }
      
      @media (max-width: 700px) {

        .search-panel-box {
          padding: 24px 22px 24px 22px;
        }

        .search-panel-box h4 {
          font-size: 26px;
        }
        
        .golden-banner {
          padding: 22px 18px 28px 18px;
        }

        .golden-banner-grid {
          max-width: 100%;
          gap: 24px;
        }

        .golden-banner-card-wrap {
          flex: 0 0 280px;
          width: 280px;
        }

        .golden-banner-card {
          width: 280px;
          min-height: 265px;
          padding: 20px 20px 18px 20px;
        }

        .golden-banner-img {
          height: 88px;
          margin-bottom: 10px;
        }

        .golden-banner-secondary-btn {
          min-width: 255px;
        }

        .golden-banner-title {
          font-size: 24px;
        }

        .golden-banner-subtitle {
          font-size: 15px;
        }

        .golden-banner .help-btn {
          width: 32px;
          height: 32px;
          line-height: 27px;
          font-size: 19px;
        }

        .landing-search-card label {
          font-size: 16px;
        }

        .landing-search-card textarea.form-control {
          min-height: 240px;
          font-size: 17px;
          padding: 12px 14px;
        }

        .landing-search-card textarea.form-control::placeholder {
          font-size: 17px;
        }

        .landing-search-card .section-helper-text {
          font-size: 14px;
        }

        .landing-summary-box {
          max-width: 390px;
          padding: 18px 22px;
        }

        .landing-summary-row {
          font-size: 13px;
          gap: 10px;
        }
      }
    "))),
    
    tags$script(HTML("
      function positionDioptMaxMarkers() {
        var wrap = $('#diopt_slider_wrap');
        var layer = $('#diopt_max_marker_layer');
        if (wrap.length === 0 || layer.length === 0) return;
        var layerLeft = layer.offset().left;
        layer.find('.diopt-max-marker').each(function() {
          var marker = $(this);
          var value = String(marker.data('value'));
          var exactLabel = wrap.find('.irs-grid-text').filter(function() {
            return $.trim($(this).text()) === value;
          }).first();
          if (exactLabel.length > 0) {
            var x = exactLabel.offset().left + exactLabel.outerWidth() / 2 - layerLeft;
            marker.css('left', x + 'px');
          } else {
            var line = wrap.find('.irs-line').first();
            var minText = wrap.find('.irs-min').first().text();
            var maxText = wrap.find('.irs-max').first().text();
            var minVal = parseFloat(minText);
            var maxVal = parseFloat(maxText);
            var thisVal = parseFloat(value);
            if (line.length > 0 && !isNaN(minVal) && !isNaN(maxVal) && !isNaN(thisVal)) {
              var lineLeft = line.offset().left;
              var lineWidth = line.outerWidth();
              var pct = (thisVal - minVal) / (maxVal - minVal);
              var x = lineLeft + pct * lineWidth - layerLeft;
              marker.css('left', x + 'px');
            }
          }
        });
      }
      $(document).on('shiny:bound shiny:value', function() { setTimeout(positionDioptMaxMarkers, 100); });
      $(document).ready(function() { setTimeout(positionDioptMaxMarkers, 500); });
      $(window).on('resize', function() { positionDioptMaxMarkers(); });
      $(document).on('toggle', '#advanced_search_details', function() {
  if (this.open) {
    setTimeout(function() {
      $(window).trigger('resize');
      positionDioptMaxMarkers();
    }, 250);
  }
});
    "))
  ),
  
  tags$div(
    class = "app-header-title",
    "Integrative Inherited Metabolic Disease Gene and Model Organism Explorer"
  ),
  tags$div(
    class = "top-nav-wrap",
    uiOutput("top_nav")
  ),
  uiOutput("page_body")
)

# -------------------------
# 10) SERVER
# -------------------------
server <- function(input, output, session) {
  
  session_start_time <- Sys.time()
  message("[TIMER] browser session connected")
  
  session$onFlushed(function() {
    message(
      "[TIMER] first browser flush finished in ",
      round(as.numeric(difftime(Sys.time(), session_start_time, units = "secs")), 3),
      " sec"
    )
  }, once = TRUE)
  
  start_timer <- function(label) {
    t0 <- Sys.time()
    force(label)
    
    function() {
      message(
        "[TIMER] ",
        label,
        " took ",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3),
        " sec"
      )
    }
  }
  
  rv <- reactiveValues(
    search_results = imd_all[0, ],
    selected_record_id = NULL,
    selected_category = NULL,
    selected_model = NULL,
    selected_golden_model = NULL,
    selected_relationships = c("1:1 Orthology", "1:2 Orthology", "Complex Orthology"),
    selected_pathway_name = NULL,
    selected_pathway_id = NULL,
    selected_pathway_model = NULL,
    selected_click_box = NULL,
    selected_click_genes = NULL,
    highlight_gene = NULL,
    summary_mode = NULL,
    golden_banner_collapsed = FALSE,
    disease_info_collapsed = FALSE,
    research_section_collapsed = FALSE,
    model_section_collapsed = FALSE,
    search_results_visible = FALSE,
    content_mode = "landing",
    landing_mode = "search",
    main_tab_selected = "Search Results",
    
    # store the user's last selected DIOPT threshold
    saved_diopt_threshold = diopt_default
  )
  
  observeEvent(input$diopt_threshold, {
    if (!is.null(input$diopt_threshold) && !is.na(input$diopt_threshold)) {
      rv$saved_diopt_threshold <- as.numeric(input$diopt_threshold)
    }
  }, ignoreInit = TRUE)
  
  current_diopt_threshold <- reactive({
    threshold <- input$diopt_threshold
    
    if (is.null(threshold) || is.na(threshold)) {
      threshold <- rv$saved_diopt_threshold
    }
    
    if (is.null(threshold) || is.na(threshold)) {
      threshold <- diopt_default
    }
    
    as.numeric(threshold)
  })
  
  show_landing_page <- function(mode) {
    rv$content_mode <- "landing"
    rv$landing_mode <- mode
    rv$selected_record_id <- NULL
    rv$selected_category <- NULL
    rv$selected_model <- NULL
    rv$selected_golden_model <- NULL
    rv$selected_click_box <- NULL
    rv$selected_click_genes <- NULL
    rv$summary_mode <- NULL
    reset_selected_pathway()
  }
  
  output$top_nav <- renderUI({
    
    nav_button <- function(input_id, label, mode, extra_class = "") {
      
      active_class <- if (
        identical(rv$content_mode, "landing") &&
        identical(rv$landing_mode, mode)
      ) {
        " active-mode"
      } else {
        ""
      }
      
      actionButton(
        input_id,
        label,
        class = paste(
          "landing-nav-btn",
          extra_class,
          active_class
        )
      )
    }
    
    tags$div(
      class = "landing-nav-bar",
      
      nav_button(
        input_id = "nav_disease_search",
        label = "Disease Search",
        mode = "search"
      ),
      
      nav_button(
        input_id = "nav_browse_models",
        label = "Browse by Models",
        mode = "models"
      ),
      
      nav_button(
        input_id = "nav_icimd",
        label = "ICIMD Classification",
        mode = "icimd"
      ),
      
      nav_button(
        input_id = "nav_kegg",
        label = "IMD-KEGG Classification",
        mode = "kegg",
        extra_class = "landing-nav-btn-wide"
      ),
      
      nav_button(
        input_id = "nav_more_info",
        label = "More Info",
        mode = "more",
        extra_class = "landing-nav-btn-small"
      )
    )
  })
  
  observeEvent(input$nav_disease_search, {
    show_landing_page("search")
  }, ignoreInit = TRUE)
  
  observeEvent(input$nav_browse_models, {
    show_landing_page("models")
  }, ignoreInit = TRUE)
  
  observeEvent(input$nav_icimd, {
    show_landing_page("icimd")
  }, ignoreInit = TRUE)
  
  observeEvent(input$nav_kegg, {
    show_landing_page("kegg")
  }, ignoreInit = TRUE)
  
  observeEvent(input$nav_more_info, {
    show_landing_page("more")
  }, ignoreInit = TRUE)
  
  output$page_body <- renderUI({
    tagList(
      tags$div(
        id = "landing_page_wrap",
        style = if (identical(rv$content_mode, "landing")) "display:block;" else "display:none;",
        uiOutput("landing_content")
      ),
      tags$div(
        id = "results_page_wrap",
        class = "results-page-wrap",
        style = if (identical(rv$content_mode, "results")) "display:block;" else "display:none;",
        tabsetPanel(
          id = "main_tabs",
          selected = rv$main_tab_selected,
          
          tabPanel(
            title = "Search Results",
            value = "Search Results",
            uiOutput("search_results_section")
          ),
          
          tabPanel(
            title = "IMD Disease Page",
            value = "Disease Page",
            uiOutput("disease_page")
          ),
          
          tabPanel(
            title = "IMD Pathway Map",
            value = "Pathway Map",
            uiOutput("pathway_map_page")
          )
        )
      )
    )
  })
  
  observeEvent(input$icimd_classification_help, {
    showModal(
      modalDialog(
        title = "How to use ICIMD disease classification",
        easyClose = TRUE,
        size = "m",
        footer = modalButton("Close"),
        
        tags$ul(
          tags$li("Browse IMD records using the ICIMD disease classification hierarchy."),
          tags$li(paste0("Each category shows any ortholog support across ", length(all_models), " model organisms.")),
          tags$li("The clickable coloured badges show the percentage and number of IMD genes in that category with ortholog support in each model organism."),
          tags$li("Click 'Show more' to expand a category and view its subcategories and disease entries."),
          tags$li("Click a disease entry to open the IMD Disease Page for detailed disease information and model-organism orthology results.")
        )
      )
    )
  })
  
  observeEvent(input$kegg_classification_help, {
    showModal(
      modalDialog(
        title = "How to use IMD-KEGG pathway classification",
        easyClose = TRUE,
        size = "m",
        footer = modalButton("Close"),
        
        tags$ul(
          tags$li("Explore IMD-associated genes in the context of KEGG pathway biology."),
          tags$li("KEGG pathways are grouped into broad functional categories and can be expanded step by step."),
          tags$li("Each pathway lists the IMD-associated human genes mapped to that pathway."),
          tags$li("Use the 'Map' button to visualise IMD genes on the corresponding KEGG pathway background."),
          tags$li("Select a gene to open the IMD Disease Page, where disease details and model-organism orthology support are shown.")
        )
      )
    )
  })
  
  output$landing_content <- renderUI({
    mode <- rv$landing_mode
    
    if (identical(mode, "models")) {
      return(tags$div(class = "browse-model-page", uiOutput("golden_gene_banner")))
    }
    
    if (identical(mode, "icimd")) {
      return(
        tags$div(
          class = "classification-page",
          
          tags$div(
            class = "classification-heading",
            tags$img(src = "assets/IEMbase.jpg"),
            tags$span("Disease classification"),
            actionButton(
              "icimd_classification_help",
              "?",
              class = "help-btn classification-help-btn",
              title = "How to use disease classification"
            )
          ),
          
          uiOutput("category_ui")
        )
      )
    }
    
    if (identical(mode, "kegg")) {
      return(
        tags$div(
          class = "classification-page",
          
          tags$div(
            class = "classification-heading",
            tags$img(src = "assets/KEGG.gif"),
            tags$span("Pathway classification"),
            actionButton(
              "kegg_classification_help",
              "?",
              class = "help-btn classification-help-btn",
              title = "How to use pathway classification"
            )
          ),
          
          uiOutput("category_ui")
        )
      )
    }
    
    if (identical(mode, "more")) {
      return(uiOutput("more_info_page"))
    }
    
    tags$div(
      class = "landing-search-layout",
      tags$div(
        class = "search-panel-box landing-search-card",
        h4("Search"),
        textAreaInput(
          "search_text",
          label = "Search by Human Symbol, ICIEM Name, Alternative Name, or Disease Abbreviation",
          placeholder = paste(
            "Try one of these examples:",
            "• PAH",
            "• OTC deficiency",
            "• CPS1",
            "",
            "You may enter multiple queries, one per line.",
            sep = "\n"
          ),
          rows = 8,
          width = "100%"
        ),
        tags$details(
          id = "advanced_search_details",
          class = "advanced-search-box",
          
          tags$summary(
            class = "advanced-search-summary",
            tags$span(class = "advanced-search-title", "Advanced Search"),
            tags$span(
              class = "advanced-search-subtitle",
              "Set the global DIOPT threshold for orthology-based results"
            )
          ),
          
          tags$div(
            class = "advanced-search-content",
            
            tags$div(
              style = "margin-top:8px; margin-bottom:10px;",
              
              tags$div(
                style = "display:flex; align-items:center; gap:8px; margin-bottom:4px;",
                tags$strong("DIOPT score threshold", style = "font-size:18px;"),
                actionButton(
                  "diopt_score_help",
                  "?",
                  class = "help-btn",
                  title = "What is DIOPT score?"
                )
              ),
              
              tags$div(
                id = "diopt_slider_wrap",
                class = "diopt-slider-wrap",
                sliderInput(
                  "diopt_threshold",
                  label = NULL,
                  min = min(available_diopt_thresholds, na.rm = TRUE),
                  max = max(available_diopt_thresholds, na.rm = TRUE),
                  value = isolate(rv$saved_diopt_threshold),
                  step = 1,
                  ticks = TRUE,
                  width = "100%"
                ),
                make_diopt_max_markers(diopt_species_max)
              ),
              
              tags$ul(
                class = "advanced-search-notes",
                tags$li("Default threshold is 7."),
                tags$li("Increasing the threshold makes orthology support more stringent and retains only stronger DIOPT-supported relationships."),
                tags$li("Decreasing the threshold is less stringent and includes weaker DIOPT-supported relationships."),
                tags$li(
                  tags$strong("This setting is applied globally: "),
                  "the selected DIOPT threshold is used for all orthology-based search and browse results in this session, including Browse by Models, classification browsing, and IMD Disease Page model outputs."
                )
              )
            )
          )
        ),
        tags$div(
          style = "display:flex; justify-content:space-between; align-items:flex-start; margin-top:6px;",
          actionButton("search_btn", "Search", class = "btn-primary"),
          actionButton(
            "try_example_btn",
            "Try Example",
            class = "btn btn-default btn-sm",
            style = "background:#eef3f8; border-color:#cfd8e3; color:#4a5a6a; padding:4px 12px; margin-top:2px;"
          )
        )
      ),
      uiOutput("landing_one_to_one_summary")
    )
  })
  
  output$more_info_page <- renderUI({
    tags$div(
      class = "more-info-page",
      tags$div(
        class = "summary-records-box",
        tags$div(
          class = "summary-records-title",
          "Summary of IMD records"
        ),
        tags$div(
          class = "summary-records-btn-row",
          actionButton(
            "show_all_imd_records",
            "Show all IMD Records",
            class = "btn btn-default summary-records-btn"
          ),
          actionButton(
            "show_childhood_dementia",
            "Show Childhood Dementia",
            class = "btn btn-default summary-records-btn"
          )
        )
      ),
      tags$div(class = "ortholog-classification-title", "Classification of Orthologs"),
      if (nzchar(ortholog_classification_src)) {
        tags$img(src = ortholog_classification_src, class = "ortholog-classification-img")
      } else {
        tags$div(
          class = "search-note",
          "Ortholog_classification image was not found. Put Ortholog_classification.png in the app folder or in www/."
        )
      }
    )
  })
  
  get_golden_model_diseases <- function(display_model_name) {
    model_key <- case_when(
      display_model_name == "Mus musculus" ~ "Mouse",
      display_model_name == "Danio rerio" ~ "Zebrafish",
      display_model_name == "Drosophila melanogaster" ~ "Fly",
      display_model_name == "Caenorhabditis elegans" ~ "Worm",
      display_model_name == "Schizosaccharomyces pombe" ~ "Fission Yeast",
      display_model_name == "Saccharomyces cerevisiae" ~ "Budding Yeast",
      TRUE ~ NA_character_
    )
    
    req(!is.na(model_key))
    
    golden_genes <- orth_all_current() %>%
      filter(
        Model == model_key,
        clean_text(Orthology_Relationship) == "1:1"
      ) %>%
      distinct(Human_Symbol) %>%
      pull(Human_Symbol)
    
    imd_all %>%
      filter(Gene %in% golden_genes) %>%
      distinct(record_id, .keep_all = TRUE) %>%
      arrange(Disease_Category, Disease_Subcategory, ICIEM_Name) %>%
      mutate(
        Selected_Orthology_Relationship = "1:1",
        Relationship_Group = "1:1 Orthology"
      )
  }
  
  get_ortholog_model_diseases <- function(display_model_name) {
    model_key <- case_when(
      display_model_name == "Mus musculus" ~ "Mouse",
      display_model_name == "Danio rerio" ~ "Zebrafish",
      display_model_name == "Drosophila melanogaster" ~ "Fly",
      display_model_name == "Caenorhabditis elegans" ~ "Worm",
      display_model_name == "Schizosaccharomyces pombe" ~ "Fission Yeast",
      display_model_name == "Saccharomyces cerevisiae" ~ "Budding Yeast",
      TRUE ~ NA_character_
    )
    
    req(!is.na(model_key))
    
    ortholog_genes <- orth_all_current() %>%
      filter(
        Model == model_key,
        Has_Ortholog
      ) %>%
      distinct(Human_Symbol) %>%
      pull(Human_Symbol)
    
    imd_all %>%
      filter(Gene %in% ortholog_genes) %>%
      distinct(record_id, .keep_all = TRUE) %>%
      arrange(Disease_Category, Disease_Subcategory, ICIEM_Name) %>%
      mutate(
        Selected_Orthology_Relationship = "With ortholog support",
        Relationship_Group = "Ortholog support"
      )
  }
  
  output$search_results_section <- renderUI({
    
    done <- start_timer("output$search_results_section")
    on.exit(done(), add = TRUE)
    
    tagList(
      tags$div(id = "search_results_top"),
      h3("Search Results"),
      uiOutput("search_filter_info"),
      actionButton("clear_search_filter", "Reset filter"),
      br(), br(),
      
      if (isTRUE(rv$search_results_visible)) {
        tagList(
          uiOutput("relationship_filter_ui"),
          br(),
          DTOutput("search_table"),
          uiOutput("search_result_note"),
          br()
        )
      }
    )
  })
  
  observeEvent(input$toggle_disease_info_section, {
    rv$disease_info_collapsed <- !isTRUE(rv$disease_info_collapsed)
  })
  
  observeEvent(input$toggle_research_section, {
    rv$research_section_collapsed <- !isTRUE(rv$research_section_collapsed)
  })
  
  observeEvent(input$toggle_model_section, {
    rv$model_section_collapsed <- !isTRUE(rv$model_section_collapsed)
  })
  
  observeEvent(rv$selected_record_id, {
    rv$disease_info_collapsed <- FALSE
    rv$research_section_collapsed <- FALSE
    rv$model_section_collapsed <- FALSE
  }, ignoreInit = TRUE)
  
  observeEvent(input$clear_nav_mode, {
    updateRadioButtons(
      session,
      "nav_mode",
      selected = character(0)
    )
  })
  
  orth_all_current <- reactive({
    threshold <- current_diopt_threshold()
    
    threshold_key <- as.character(as.numeric(threshold))
    x <- orth_by_threshold[[threshold_key]]
    
    if (is.null(x)) {
      return(orth_all_raw[0, ])
    }
    
    x
  })
  
  category_coverage_full_current <- reactive({
    
    done <- start_timer("category_coverage_full_current")
    on.exit(done(), add = TRUE)
    
    ortholog_presence_by_model <- orth_all_current() %>%
      mutate(
        Human_Symbol = clean_text(Human_Symbol),
        Model = clean_text(Model),
        Orthology_Relationship = clean_text(Orthology_Relationship)
      ) %>%
      filter(Human_Symbol != "", Model != "") %>%
      group_by(Human_Symbol, Model) %>%
      summarise(
        Has_Ortholog = any(Has_Ortholog, na.rm = TRUE),
        .groups = "drop"
      )
    
    category_gene_model_table <- category_gene_table %>%
      tidyr::crossing(Model = all_models) %>%
      left_join(
        ortholog_presence_by_model,
        by = c("Gene" = "Human_Symbol", "Model" = "Model")
      ) %>%
      mutate(
        Has_Ortholog = replace_na(Has_Ortholog, FALSE)
      )
    
    category_coverage <- category_gene_model_table %>%
      group_by(Disease_Category, Model) %>%
      summarise(
        A = n_distinct(Gene),
        B = n_distinct(Gene[Has_Ortholog]),
        coverage_pct = round(ifelse(A > 0, 100 * B / A, 0), 1),
        .groups = "drop"
      )
    
    tidyr::expand_grid(
      Disease_Category = all_categories,
      Model = all_models
    ) %>%
      left_join(category_coverage, by = c("Disease_Category", "Model")) %>%
      mutate(
        A = replace_na(A, 0),
        B = replace_na(B, 0),
        coverage_pct = replace_na(coverage_pct, 0)
      )
  })
  
  observe({
    lapply(model_display_order, function(display_model_name) {
      local({
        this_model <- display_model_name
        this_id <- paste0("ortholog_model_", safe_id(this_model))
        
        observeEvent(input[[this_id]], {
          rv$selected_category <- NULL
          rv$selected_model <- NULL
          rv$selected_golden_model <- paste0(this_model, " | all ortholog support")
          rv$selected_click_box <- NULL
          rv$selected_click_genes <- NULL
          rv$selected_record_id <- NULL
          rv$summary_mode <- NULL
          rv$selected_relationships <- c("1:1 Orthology", "1:2 Orthology", "Complex Orthology")
          
          rv$search_results <- get_ortholog_model_diseases(this_model)
          
          rv$search_results_visible <- TRUE
          rv$content_mode <- "results"
          rv$main_tab_selected <- "Search Results"
          updateTabsetPanel(session, "main_tabs", selected = "Search Results")
          scroll_to_search_results(300)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  
  observeEvent(input$golden_genes_help, {
    showModal(
      modalDialog(
        title = "IMD genes: how to interpret this section",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
        tags$ul(
          tags$li("Orthology support is a guide for model prioritisation, not the only basis for model selection."),
          tags$li("1:1 orthology suggests a clearer gene match, but does not guarantee full biological consistency."),
          tags$li("Gene function, phenotype, and disease relevance may still differ across species."),
          tags$li("Complex orthology relationships (for example, one-to-many or many-to-many) require careful interpretation."),
          tags$li("Model selection should also consider biological context, known phenotypes, and experimental feasibility.")
        )
      )
    )
  })
  
  observeEvent(input$diopt_score_help, {
    showModal(
      modalDialog(
        title = "What is DIOPT score?",
        easyClose = TRUE,
        size = "l",
        footer = modalButton("Close"),
        
        tags$div(
          tags$p(
            "DIOPT is an integrative orthology prediction resource. It combines ortholog predictions from multiple published tools to help users identify high-confidence orthologs, or to explore broader sets of possible orthologs for a gene of interest."
          ),
          
          tags$p(
            "DIOPT integrates ortholog predictions for human, mouse, fly, worm, zebrafish, fission yeast, and budding yeast from resources including Ensembl Compara, HomoloGene, Inparanoid, Isobase, OMA, OrthoMCL, Phylome, RoundUp, and TreeFam."
          ),
          
          tags$p(
            "The DIOPT score is a simple support score. It indicates how many independent orthology prediction tools support a given orthologous gene-pair relationship. A higher DIOPT score generally means stronger support from multiple prediction methods."
          ),
          
          tags$p(
            "In this tool, the DIOPT score threshold controls how stringent the orthology filtering is. Increasing the threshold keeps only more strongly supported ortholog relationships, while decreasing the threshold includes weaker or less consistently supported predictions."
          ),
          
          tags$p(
            "DIOPT also provides additional information such as protein and domain alignments, including amino acid identity, which can help users choose among multiple possible orthologs."
          ),
          
          tags$p(
            tags$b("More details: "),
            tags$a(
              href = "https://www.flyrnai.org/DIOPT_documentation_v10.html",
              target = "_blank",
              "DIOPT documentation"
            )
          )
        )
      )
    )
  })
  
  observeEvent(input$toggle_golden_banner, {
    rv$golden_banner_collapsed <- !isTRUE(rv$golden_banner_collapsed)
  })
  
  golden_gene_stats_current <- reactive({
    x <- orth_all_current() %>%
      mutate(
        Human_Symbol = clean_text(Human_Symbol),
        Orthology_Relationship = clean_text(Orthology_Relationship)
      ) %>%
      inner_join(
        imd_all %>%
          filter(Gene != "") %>%
          distinct(Gene),
        by = c("Human_Symbol" = "Gene")
      ) %>%
      filter(Orthology_Relationship == "1:1") %>%
      group_by(Model) %>%
      summarise(
        golden_gene_n = n_distinct(Human_Symbol),
        .groups = "drop"
      ) %>%
      mutate(
        Display_Model = unname(model_display_map[Model])
      ) %>%
      select(Display_Model, golden_gene_n)
    
    tibble(Display_Model = model_display_order) %>%
      left_join(x, by = "Display_Model") %>%
      mutate(
        golden_gene_n = replace_na(golden_gene_n, 0L),
        total_imd_gene_n = total_imd_gene_n,
        golden_gene_pct = ifelse(
          total_imd_gene_n > 0,
          100 * golden_gene_n / total_imd_gene_n,
          0
        ),
        golden_gene_label = paste0(
          golden_gene_n, " / ", total_imd_gene_n,
          " (", sprintf("%.1f", golden_gene_pct), "%)"
        ),
        figure_src = vapply(Display_Model, function(x) {
          out <- get_model_figure_rel(x)
          if (is.null(out)) "" else out
        }, character(1))
      )
  })
  
  
  ortholog_gene_stats_current <- reactive({
    x <- orth_all_current() %>%
      mutate(
        Human_Symbol = clean_text(Human_Symbol),
        Orthology_Relationship = clean_text(Orthology_Relationship)
      ) %>%
      inner_join(
        imd_all %>%
          filter(Gene != "") %>%
          distinct(Gene),
        by = c("Human_Symbol" = "Gene")
      ) %>%
      filter(Has_Ortholog) %>%
      group_by(Model) %>%
      summarise(
        ortholog_gene_n = n_distinct(Human_Symbol),
        .groups = "drop"
      ) %>%
      mutate(
        Display_Model = unname(model_display_map[Model])
      ) %>%
      select(Display_Model, ortholog_gene_n)
    
    tibble(Display_Model = model_display_order) %>%
      left_join(x, by = "Display_Model") %>%
      mutate(
        ortholog_gene_n = replace_na(ortholog_gene_n, 0L),
        total_imd_gene_n = total_imd_gene_n,
        ortholog_gene_pct = ifelse(
          total_imd_gene_n > 0,
          100 * ortholog_gene_n / total_imd_gene_n,
          0
        ),
        ortholog_gene_label = paste0(
          ortholog_gene_n, " / ", total_imd_gene_n,
          " (", sprintf("%.1f", ortholog_gene_pct), "%)"
        )
      )
  })
  
  
  output$landing_one_to_one_summary <- renderUI({
    df <- golden_gene_stats_current()
    threshold_value <- if (!is.null(input$diopt_threshold)) input$diopt_threshold else diopt_default
    
    tags$div(
      class = "landing-summary-box",
      tags$div(
        class = "landing-summary-title",
        paste0("Total IMD genes: ", total_imd_gene_n)
      ),
      tags$div(
        class = "landing-summary-threshold",
        HTML(paste0("DIOPT Score &ge; ", threshold_value))
      ),
      tags$div(class = "landing-summary-subtitle", "One-to-one ortholog"),
      lapply(seq_len(nrow(df)), function(i) {
        tags$div(
          class = "landing-summary-row",
          tags$div(class = "landing-summary-species", df$Display_Model[i]),
          tags$div(
            class = "landing-summary-value",
            paste0(df$golden_gene_n[i], " (", sprintf("%.1f", df$golden_gene_pct[i]), "%)")
          )
        )
      })
    )
  })
  
  observe({
    lapply(model_display_order, function(display_model_name) {
      local({
        this_model <- display_model_name
        this_id <- paste0("golden_model_", safe_id(this_model))
        
        observeEvent(input[[this_id]], {
          rv$selected_category <- NULL
          rv$selected_model <- NULL
          rv$selected_golden_model <- this_model
          rv$selected_click_box <- NULL
          rv$selected_click_genes <- NULL
          rv$selected_record_id <- NULL
          rv$summary_mode <- NULL
          rv$selected_relationships <- "1:1 Orthology"
          
          rv$search_results <- get_golden_model_diseases(this_model)
          
          rv$search_results_visible <- TRUE
          rv$content_mode <- "results"
          rv$main_tab_selected <- "Search Results"
          updateTabsetPanel(session, "main_tabs", selected = "Search Results")
          scroll_to_search_results(300)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  get_category_model_diseases <- function(category_name, model_name) {
    req(category_name, model_name)
    
    model_map <- orth_all_current() %>%
      filter(Model == model_name, Has_Ortholog) %>%
      transmute(
        Gene = Human_Symbol,
        Selected_Model = Model,
        Selected_Orthology_Relationship = Orthology_Relationship
      ) %>%
      distinct()
    
    imd_all %>%
      filter(Disease_Category == category_name) %>%
      left_join(model_map, by = "Gene") %>%
      mutate(
        Selected_Orthology_Relationship = replace_na(Selected_Orthology_Relationship, ""),
        Relationship_Group = case_when(
          Selected_Orthology_Relationship == "1:1" ~ "1:1 Orthology",
          Selected_Orthology_Relationship == "1:2" ~ "1:2 Orthology",
          Selected_Orthology_Relationship == "" ~ "No ortholog",
          TRUE ~ "Complex Orthology"
        )
      ) %>%
      filter(Relationship_Group != "No ortholog") %>%
      arrange(Disease_Subcategory, ICIEM_Name)
  }
  
  observeEvent(input$category_chip_click, {
    req(input$category_chip_click$category)
    req(input$category_chip_click$model)
    
    this_cat <- clean_text(input$category_chip_click$category)
    this_model <- clean_text(input$category_chip_click$model)
    
    rv$selected_category <- this_cat
    rv$selected_model <- this_model
    rv$selected_click_box <- NULL
    rv$selected_click_genes <- NULL
    rv$selected_relationships <- c("1:1 Orthology", "1:2 Orthology", "Complex Orthology")
    
    rv$search_results <- get_category_model_diseases(this_cat, this_model)
    rv$search_results_visible <- TRUE
    rv$content_mode <- "results"
    
    rv$main_tab_selected <- "Search Results"
    updateTabsetPanel(session, "main_tabs", selected = "Search Results")
    scroll_to_search_results(300)
  }, ignoreInit = TRUE)
  
  
  observeEvent(input$disease_click, {
    req(input$disease_click$record_id)
    
    rid <- suppressWarnings(as.integer(input$disease_click$record_id))
    req(!is.na(rid))
    
    rv$selected_record_id <- rid
    rv$content_mode <- "results"
    reset_selected_pathway()
    
    rv$main_tab_selected <- "Disease Page"
    updateTabsetPanel(session, "main_tabs", selected = "Disease Page")
  }, ignoreInit = TRUE)
  
  
  observeEvent(input$kegg_gene_click, {
    req(input$kegg_gene_click$record_id)
    
    rid <- suppressWarnings(as.integer(input$kegg_gene_click$record_id))
    req(!is.na(rid))
    
    rv$selected_record_id <- rid
    rv$content_mode <- "results"
    
    rv$main_tab_selected <- "Disease Page"
    updateTabsetPanel(session, "main_tabs", selected = "Disease Page")
  }, ignoreInit = TRUE)
  
  
  observeEvent(input$pathway_map_pick, {
    req(input$pathway_map_pick$pathway_id)
    req(input$pathway_map_pick$pathway_name)
    
    rv$selected_pathway_name <- clean_text(input$pathway_map_pick$pathway_name)
    rv$selected_pathway_id <- clean_text(input$pathway_map_pick$pathway_id)
    rv$selected_pathway_model <- NULL
    rv$highlight_gene <- NULL
    rv$content_mode <- "results"
    
    rv$main_tab_selected <- "Pathway Map"
    updateTabsetPanel(session, "main_tabs", selected = "Pathway Map")
  }, ignoreInit = TRUE)
  
  
  observeEvent(input$pathway_model_pick, {
    req(input$pathway_model_pick$model)
    
    rv$selected_pathway_model <- clean_text(input$pathway_model_pick$model)
    rv$content_mode <- "results"
    
    rv$main_tab_selected <- "Pathway Map"
    updateTabsetPanel(session, "main_tabs", selected = "Pathway Map")
  }, ignoreInit = TRUE)
  
  output$golden_gene_banner <- renderUI({
    
    done <- start_timer("output$golden_gene_banner")
    on.exit(done(), add = TRUE)
    
    if (isTRUE(rv$golden_banner_collapsed)) {
      return(
        tags$div(
          class = "golden-banner",
          
          tags$div(
            style = "display:flex; align-items:center; justify-content:space-between; gap:12px; margin-bottom:4px;",
            
            tags$div(
              style = "display:flex; align-items:center; gap:8px;",
              tags$div(
                class = "golden-banner-title",
                "IMD genes across model organisms"
              ),
              actionButton(
                "golden_genes_help",
                "?",
                class = "help-btn",
                title = "How to interpret IMD genes across model organisms"
              )
            ),
            
            actionButton(
              "toggle_golden_banner",
              "Expand",
              class = "btn btn-default btn-sm golden-collapse-btn"
            )
          ),
          
          tags$div(
            class = "golden-banner-subtitle",
            "Expand to view model-level IMD orthology support."
          )
        )
      )
    }
    
    banner_df <- golden_gene_stats_current() %>%
      left_join(
        ortholog_gene_stats_current() %>%
          select(Display_Model, ortholog_gene_n, ortholog_gene_label),
        by = "Display_Model"
      )
    
    tags$div(
      class = "golden-banner",
      
      tags$div(
        style = "display:flex; align-items:center; justify-content:space-between; gap:12px; margin-bottom:4px;",
        
        tags$div(
          style = "display:flex; align-items:center; gap:8px;",
          tags$div(
            class = "golden-banner-title",
            "IMD genes across model organisms"
          ),
          actionButton(
            "golden_genes_help",
            "?",
            class = "help-btn",
            title = "How to interpret IMD genes across model organisms"
          )
        ),
        
        actionButton(
          "toggle_golden_banner",
          if (isTRUE(rv$golden_banner_collapsed)) "Expand" else "Collapse",
          class = "btn btn-default btn-sm golden-collapse-btn"
        )
      ),
      
      if (!isTRUE(rv$golden_banner_collapsed)) {
        tagList(
          tags$div(
            class = "golden-banner-subtitle",
            tags$span(
              class = "golden-highlight-text",
              "Unique IMD-causing genes with orthology support"
            ),
            " in each model. Click the main card to view 1:1-supported entries, or use the secondary button to view all IMD genes with ortholog support in that model."
          ),
          
          tags$div(
            class = "golden-banner-grid",
            lapply(seq_len(nrow(banner_df)), function(i) {
              row <- banner_df[i, ]
              golden_click_id <- paste0("golden_model_", safe_id(row$Display_Model))
              ortholog_click_id <- paste0("ortholog_model_", safe_id(row$Display_Model))
              
              tags$div(
                class = "golden-banner-card-wrap",
                
                tags$a(
                  href = "javascript:void(0);",
                  class = "golden-banner-card golden-banner-card-link",
                  onclick = sprintf(
                    "Shiny.setInputValue('%s', Math.random(), {priority: 'event'})",
                    golden_click_id
                  ),
                  if (nzchar(row$figure_src)) {
                    tags$img(src = row$figure_src, class = "golden-banner-img")
                  },
                  tags$div(class = "golden-banner-model", row$Display_Model),
                  tags$div(class = "golden-banner-number", row$golden_gene_n),
                  tags$div(class = "golden-banner-support-label", "1:1 ortholog support"),
                  tags$div(class = "golden-banner-prop", row$golden_gene_label)
                ),
                
                tags$a(
                  href = "javascript:void(0);",
                  class = "golden-banner-secondary-btn",
                  onclick = sprintf(
                    "Shiny.setInputValue('%s', Math.random(), {priority: 'event'})",
                    ortholog_click_id
                  ),
                  title = "View IMD genes with any ortholog support in this model",
                  tags$div(class = "golden-banner-secondary-title", "Total ortholog support"),
                  tags$div(class = "golden-banner-secondary-value", row$ortholog_gene_label)
                )
              )
            })
          )
        )
      }
    )
  })
  
  output$relationship_filter_ui <- renderUI({
    
    done <- start_timer("output$relationship_filter_ui")
    on.exit(done(), add = TRUE)
    
    if (!is.null(rv$selected_category) && !is.null(rv$selected_model) && nrow(rv$search_results) > 0) {
      tags$div(
        class = "relationship-filter-wrap",
        tags$div(class = "relationship-filter-title", "Human : Model Organism"),
        div(
          class = "relationship-checkbox-group",
          checkboxGroupInput(
            "relationship_filter",
            label = NULL,
            choices = c("1:1 Orthology", "1:2 Orthology", "Complex Orthology"),
            selected = rv$selected_relationships,
            inline = TRUE
          )
        )
      )
    }
  })
  
  
  scroll_to_search_results <- function(delay = 250) {
    shinyjs::runjs(sprintf("
    setTimeout(function() {
      var el = document.getElementById('search_results_top');
      if (el) {
        el.scrollIntoView({behavior: 'smooth', block: 'start'});
      }
    }, %d);
  ", delay))
  }
  
  observeEvent(input$back_to_pathway_models, {
    rv$selected_pathway_model <- NULL
    rv$main_tab_selected <- "Pathway Map"
    updateTabsetPanel(session, "main_tabs", selected = "Pathway Map")
  })
  
  observeEvent(input$relationship_filter, {
    rv$selected_relationships <- input$relationship_filter
  }, ignoreNULL = FALSE)
  
  observeEvent(input$try_example_btn, {
    example_queries <- c(
      "PAH",
      "PKU",
      "CPS1",
      "OTC deficiency",
      "ASS1"
    )
    
    chosen_example <- sample(example_queries, 1)
    
    updateTextAreaInput(
      session,
      "search_text",
      value = chosen_example
    )
  })
  
  observeEvent(input$clear_search_filter, {
    rv$selected_category <- NULL
    rv$selected_model <- NULL
    rv$selected_golden_model <- NULL
    rv$selected_click_box <- NULL
    rv$selected_click_genes <- NULL
    rv$selected_relationships <- c("1:1 Orthology", "1:2 Orthology", "Complex Orthology")
    rv$search_results <- imd_all[0, ]
    rv$summary_mode <- NULL
    rv$search_results_visible <- FALSE
    rv$content_mode <- "landing"
    rv$landing_mode <- "search"
  })
  
  observeEvent(input$gene_to_pathway, {
    dx <- selected_disease()
    req(nrow(dx) == 1)
    
    gene_symbol <- dx$Gene[[1]]
    pw <- get_gene_pathways(gene_symbol)
    
    if (nrow(pw) == 0) {
      showModal(
        modalDialog(
          title = "No pathway map found",
          paste("No KEGG pathway map was found for gene", gene_symbol),
          easyClose = TRUE,
          footer = modalButton("Close")
        )
      )
      return(NULL)
    }
    
    # if only one pathway, jump directly
    if (nrow(pw) == 1) {
      rv$selected_pathway_name <- pw$KEGG_Pathway[[1]]
      rv$selected_pathway_id   <- pw$Pathway_KEGG_ID[[1]]
      rv$selected_pathway_model <- NULL
      rv$highlight_gene <- gene_symbol
      rv$content_mode <- "results"
      rv$main_tab_selected <- "Pathway Map"
      updateTabsetPanel(session, "main_tabs", selected = "Pathway Map")
      return(NULL)
    }
    
    # if multiple pathways, let user choose
    pathway_buttons <- lapply(seq_len(nrow(pw)), function(i) {
      this_id <- paste0("gene_pathway_choice_", i)
      
      actionButton(
        inputId = this_id,
        label = paste0(pw$KEGG_Pathway[i], " (", pw$Pathway_KEGG_ID[i], ")"),
        class = "btn btn-default",
        style = "display:block; margin-bottom:8px; width:100%; text-align:left;"
      )
    })
    
    showModal(
      modalDialog(
        title = paste("Select pathway map for", gene_symbol),
        easyClose = TRUE,
        footer = modalButton("Close"),
        do.call(tagList, pathway_buttons)
      )
    )
    
    lapply(seq_len(nrow(pw)), function(i) {
      local({
        ii <- i
        this_id <- paste0("gene_pathway_choice_", ii)
        
        observeEvent(input[[this_id]], {
          removeModal()
          rv$selected_pathway_name <- pw$KEGG_Pathway[[ii]]
          rv$selected_pathway_id   <- pw$Pathway_KEGG_ID[[ii]]
          rv$selected_pathway_model <- NULL
          rv$highlight_gene <- gene_symbol
          rv$content_mode <- "results"
          rv$main_tab_selected <- "Pathway Map"
          updateTabsetPanel(session, "main_tabs", selected = "Pathway Map")
        }, ignoreInit = TRUE, once = TRUE)
      })
    })
  })
  
  
  observeEvent(input$show_childhood_dementia, {
    rv$selected_category <- NULL
    rv$selected_model <- NULL
    rv$selected_golden_model <- NULL
    rv$search_results_visible <- TRUE
    rv$content_mode <- "results"
    rv$selected_click_box <- NULL
    rv$selected_click_genes <- NULL
    rv$selected_record_id <- NULL
    rv$selected_relationships <- c("1:1 Orthology", "1:2 Orthology", "Complex Orthology")
    rv$summary_mode <- "childhood_dementia"
    rv$search_results <- get_childhood_dementia_records()
    
    rv$main_tab_selected <- "Search Results"
    updateTabsetPanel(session, "main_tabs", selected = "Search Results")
  })
  
  observeEvent(input$show_all_imd_records, {
    rv$selected_category <- NULL
    rv$selected_model <- NULL
    rv$selected_golden_model <- NULL
    rv$search_results_visible <- TRUE
    rv$content_mode <- "results"
    rv$selected_click_box <- NULL
    rv$selected_click_genes <- NULL
    rv$selected_record_id <- NULL
    rv$selected_relationships <- c("1:1 Orthology", "1:2 Orthology", "Complex Orthology")
    rv$summary_mode <- "all_imd_records"
    rv$search_results <- imd_all
    
    rv$main_tab_selected <- "Search Results"
    updateTabsetPanel(session, "main_tabs", selected = "Search Results")
  })
  
  # initial summary
  output$summary_table <- renderDT({
    if (is.null(rv$summary_mode) || identical(rv$summary_mode, "")) {
      df <- imd_all[0, ]
    } else if (identical(rv$summary_mode, "all")) {
      df <- imd_all
    } else if (identical(rv$summary_mode, "childhood_dementia")) {
      df <- get_childhood_dementia_records()
    } else {
      df <- imd_all[0, ]
    }
    
    df <- df %>%
      select(
        record_id,
        Disease_Category,
        Disease_Subcategory,
        ICIEM_Name,
        Gene,
        Disease_Abbreviation,
        Alternative_Names
      )
    
    datatable(
      df,
      selection = "single",
      rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE)
    )
  })
  
  output$search_filter_info <- renderUI({
    if (!is.null(rv$selected_click_box) && !is.null(rv$selected_pathway_id)) {
      n_hits <- nrow(rv$search_results)
      gene_text <- if (!is.null(rv$selected_click_genes) && nzchar(rv$selected_click_genes)) {
        paste0(" | Gene(s): ", rv$selected_click_genes)
      } else {
        ""
      }
      
      return(
        tags$p(
          tags$b("Current filter: "),
          paste0(
            rv$selected_pathway_id,
            " | ",
            n_hits,
            " disease record(s)",
            gene_text
          )
        )
      )
    }
    
    if (!is.null(rv$selected_golden_model)) {
      n_hits <- nrow(rv$search_results)
      
      if (str_detect(rv$selected_golden_model, fixed(" | all ortholog support"))) {
        this_model <- str_remove(rv$selected_golden_model, fixed(" | all ortholog support"))
        return(
          tags$p(
            tags$b("Current filter: "),
            paste0(
              "Ortholog support | ",
              this_model,
              " | ",
              n_hits,
              " disease record(s) with ortholog support."
            )
          )
        )
      }
      
      return(
        tags$p(
          tags$b("Current filter: "),
          paste0(
            "IMD genes with 1:1 orthology support | ",
            rv$selected_golden_model,
            " | ",
            n_hits,
            " disease record(s) with 1:1 orthology support."
          )
        )
      )
    }
    
    if (!is.null(rv$selected_category) && !is.null(rv$selected_model)) {
      n_hits <- nrow(rv$search_results)
      return(
        tags$p(
          tags$b("Current filter: "),
          paste0(
            rv$selected_category,
            " | ",
            rv$selected_model,
            " | ",
            n_hits,
            " disease record(s) with ortholog support. Use the orthology filter below to refine the list."
          )
        )
      )
    }
    
    tags$p(class = "small-muted", "No category filter or search applied.")
  })
  
  output$search_result_note <- renderUI({
    
    done <- start_timer("output$search_result_note")
    on.exit(done(), add = TRUE)
    
    df <- filtered_search_results()
    
    if (nrow(df) == 0) {
      return(NULL)
    }
    
    # 1) Results coming from disease-classification model browse
    if (!is.null(rv$selected_category) && !is.null(rv$selected_model)) {
      return(
        tags$div(
          class = "search-note",
          tags$span(class = "search-note-icon", "!"),
          "Please review the results carefully before selecting an entry."
        )
      )
    }
    
    # 2) Results coming from pathway-map click
    if (!is.null(rv$selected_click_box) && !is.null(rv$selected_pathway_id)) {
      return(
        tags$div(
          class = "search-note",
          tags$span(class = "search-note-icon", "!"),
          "Please review the results carefully before selecting an entry."
        )
      )
    }
    
    # 3) Results coming from IMD genes banner
    if (!is.null(rv$selected_golden_model)) {
      return(
        tags$div(
          class = "search-note",
          tags$span(class = "search-note-icon", "!"),
          "Please review the results carefully before selecting an entry."
        )
      )
    }
    
    # 4) Results coming from Summary of IMD records buttons
    if (!is.null(rv$summary_mode) && rv$summary_mode %in% c("all_imd_records", "childhood_dementia")) {
      return(
        tags$div(
          class = "search-note",
          tags$span(class = "search-note-icon", "!"),
          "Please review the results carefully before selecting an entry."
        )
      )
    }
    
    # 5) Results coming from direct search
    first_hit <- clean_text(df$ICIEM_Name[[1]])
    if (first_hit == "") {
      first_hit <- clean_text(df$Gene[[1]])
    }
    
    tags$div(
      class = "search-note",
      tags$span(class = "search-note-icon", "!"),
      HTML(paste0(
        "<b>The first hit for your search is:</b> ",
        htmlEscape(first_hit),
        ". Please review the results carefully before selecting an entry. ",
        "The search function uses disease-relevant matching and may return approximate results. ",
        "If this is not the entry you intended, please refine the gene or disease name and search again."
      ))
    )
  })
  
  observeEvent(input$summary_table_rows_selected, {
    idx <- input$summary_table_rows_selected
    if (length(idx) == 1) {
      rv$selected_record_id <- imd_all$record_id[idx]
      rv$content_mode <- "results"
      rv$main_tab_selected <- "Disease Page"
      updateTabsetPanel(session, "main_tabs", selected = "Disease Page")
    }
  })
  
  observeEvent(input$disease_pathway_pick, {
    req(input$disease_pathway_pick$pathway_id)
    req(input$disease_pathway_pick$pathway_name)
    req(rv$selected_record_id)
    
    dx <- selected_disease()
    gene_symbol <- clean_text(dx$Gene[[1]])
    
    rv$selected_pathway_id <- clean_text(input$disease_pathway_pick$pathway_id)
    rv$selected_pathway_name <- clean_text(input$disease_pathway_pick$pathway_name)
    rv$selected_pathway_model <- NULL
    rv$selected_click_box <- NULL
    rv$selected_click_genes <- NULL
    rv$highlight_gene <- gene_symbol
    rv$content_mode <- "results"
    
    rv$main_tab_selected <- "Pathway Map"
    updateTabsetPanel(session, "main_tabs", selected = "Pathway Map")
  })
  
  observeEvent(input$search_btn, {
    rv$selected_category <- NULL
    rv$selected_model <- NULL
    rv$selected_golden_model <- NULL
    rv$selected_click_box <- NULL
    rv$selected_click_genes <- NULL
    rv$selected_record_id <- NULL
    rv$summary_mode <- NULL
    rv$search_results <- do_search(input$search_text, imd_all)
    rv$search_results_visible <- TRUE
    rv$content_mode <- "results"
    
    rv$main_tab_selected <- "Search Results"
    updateTabsetPanel(session, "main_tabs", selected = "Search Results")
  })
  
  filtered_search_results <- reactive({
    df <- rv$search_results
    
    if ("Relationship_Group" %in% names(df) && !is.null(rv$selected_model)) {
      df <- df %>%
        filter(Relationship_Group %in% rv$selected_relationships)
    }
    
    df
  })
  
  output$search_table <- renderDT({
    
    done <- start_timer("output$search_table")
    on.exit(done(), add = TRUE)
    
    df <- filtered_search_results()
    
    show_cols <- c(
      "record_id",
      "Disease_Category",
      "Disease_Subcategory",
      "ICIEM_Name",
      "Gene",
      "Disease_Abbreviation",
      "Alternative_Names"
    )
    
    if ("Selected_Orthology_Relationship" %in% names(df)) {
      show_cols <- c(show_cols, "Selected_Orthology_Relationship")
    }
    
    df <- df %>%
      select(any_of(show_cols))
    
    datatable(
      df,
      selection = "single",
      rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE)
    )
  })
  
  observeEvent(input$back_to_pathway_list, {
    rv$selected_pathway_name <- NULL
    rv$selected_pathway_id <- NULL
    rv$selected_pathway_model <- NULL
    rv$selected_click_box <- NULL
    rv$selected_click_genes <- NULL
    rv$highlight_gene <- NULL
    rv$content_mode <- "results"
    
    rv$main_tab_selected <- "Pathway Map"
    updateTabsetPanel(session, "main_tabs", selected = "Pathway Map")
  })
  
  observeEvent(input$search_table_rows_selected, {
    idx <- input$search_table_rows_selected
    df <- filtered_search_results()
    
    if (length(idx) == 1 && nrow(df) >= idx) {
      rv$selected_record_id <- df$record_id[idx]
      rv$content_mode <- "results"
      reset_selected_pathway()
      rv$main_tab_selected <- "Disease Page"
      updateTabsetPanel(session, "main_tabs", selected = "Disease Page")
    }
  })
  
  
  # category / subcategory / disease navigation
  output$category_ui <- renderUI({
    
    done <- start_timer("output$category_ui")
    on.exit(done(), add = TRUE)
    
    
    nav_mode_current <- dplyr::case_when(
      identical(rv$landing_mode, "icimd") ~ "disease",
      identical(rv$landing_mode, "kegg") ~ "pathway",
      TRUE ~ ""
    )
    
    if (!nzchar(nav_mode_current)) {
      return(NULL)
    }
    
    # =====================================
    # A) DISEASE CLASSIFICATION NAVIGATION
    # =====================================
    if (nav_mode_current == "disease") {
      cat_panels <- lapply(all_categories, function(cat_name) {
        sub_df <- imd_all %>%
          filter(Disease_Category == cat_name) %>%
          distinct(Disease_Subcategory, ICIEM_Name, record_id)
        
        subcat_order <- sub_df %>%
          distinct(Disease_Subcategory) %>%
          mutate(
            sub_major = suppressWarnings(as.numeric(stringr::str_extract(Disease_Subcategory, "^\\d+"))),
            sub_minor = suppressWarnings(as.numeric(stringr::str_extract(Disease_Subcategory, "(?<=\\.)\\d+"))),
            sub_major = ifelse(is.na(sub_major), 9999, sub_major),
            sub_minor = ifelse(is.na(sub_minor), 9999, sub_minor)
          ) %>%
          arrange(sub_major, sub_minor, Disease_Subcategory) %>%
          pull(Disease_Subcategory)
        
        sub_df <- sub_df %>%
          mutate(
            Disease_Subcategory = factor(Disease_Subcategory, levels = subcat_order)
          ) %>%
          arrange(Disease_Subcategory, ICIEM_Name)
        
        sub_list <- split(sub_df, as.character(sub_df$Disease_Subcategory))
        sub_list <- sub_list[subcat_order]
        
        cov_df <- category_coverage_full_current() %>%
          filter(Disease_Category == cat_name) %>%
          arrange(match(Model, all_models))
        
        chips <- lapply(seq_len(nrow(cov_df)), function(i) {
          pct <- cov_df$coverage_pct[i]
          model <- cov_df$Model[i]
          col <- coverage_colour(pct)
          
          tags$a(
            href = "javascript:void(0);",
            class = "coverage-chip",
            style = paste0(
              "background:", col, ";",
              "display:inline-block;",
              "text-decoration:none;",
              "border:none;"
            ),
            onclick = sprintf(
              "Shiny.setInputValue('category_chip_click', {category: '%s', model: '%s', nonce: Math.random()}, {priority: 'event'})",
              js_escape(cat_name),
              js_escape(model)
            ),
            paste0(model, ": ", pct, "% (", cov_df$B[i], "/", cov_df$A[i], ")")
          )
        })
        
        sub_tags <- lapply(names(sub_list), function(subcat) {
          disease_df <- sub_list[[subcat]]
          
          disease_links <- lapply(seq_len(nrow(disease_df)), function(i) {
            rid <- disease_df$record_id[i]
            dnm <- disease_df$ICIEM_Name[i]
            
            tags$a(
              href = "javascript:void(0);",
              class = "disease-link",
              onclick = sprintf(
                "Shiny.setInputValue('disease_click', {record_id: %s, nonce: Math.random()}, {priority: 'event'})",
                rid
              ),
              dnm
            )
          })
          
          tags$details(
            class = "subcategory-details",
            tags$summary(
              tags$span(class = "subcat-title", subcat),
              tags$span(class = "subcat-toggle-text")
            ),
            tags$div(
              class = "subcategory-content",
              disease_links
            )
          )
        })
        
        tags$div(
          class = "category-box",
          tags$div(
            class = "category-header-row",
            tags$h4(class = "top-title", cat_name)
          ),
          tags$div(class = "coverage-row", chips),
          tags$details(
            class = "category-details",
            tags$summary(
              tags$span(class = "toggle-label")
            ),
            tags$div(
              class = "category-content",
              sub_tags
            )
          )
        )
      })
      
      return(do.call(tagList, cat_panels))
    }
    
    # =====================================
    # B) PATHWAY CLASSIFICATION NAVIGATION
    # =====================================
    pathway_category_order <- c(
      "Metabolism",
      "Genetic Information Processing",
      "Environmental Information Processing",
      "Cellular Processes",
      "Organismal Systems",
      "Human Diseases"
    )
    
    pathway_category_order <- pathway_category_order[
      pathway_category_order %in% unique(kegg_imd_tbl$Category)
    ]
    
    pathway_panels <- lapply(pathway_category_order, function(cat_name) {
      sub_df <- kegg_imd_tbl %>%
        filter(Category == cat_name)
      
      sub_order <- sub_df %>%
        distinct(Subcategory) %>%
        arrange(Subcategory) %>%
        pull(Subcategory)
      
      sub_tags <- lapply(sub_order, function(subcat_name) {
        pathway_df <- sub_df %>%
          filter(Subcategory == subcat_name) %>%
          distinct(KEGG_Pathway, Pathway_KEGG_ID, Human_Gene_Symbol, record_id) %>%
          arrange(KEGG_Pathway, Pathway_KEGG_ID, Human_Gene_Symbol)
        
        pathway_keys <- pathway_df %>%
          distinct(KEGG_Pathway, Pathway_KEGG_ID) %>%
          arrange(KEGG_Pathway, Pathway_KEGG_ID)
        
        pathway_tags <- lapply(seq_len(nrow(pathway_keys)), function(i) {
          pw_name <- pathway_keys$KEGG_Pathway[i]
          pw_id   <- pathway_keys$Pathway_KEGG_ID[i]
          
          gene_df <- pathway_df %>%
            filter(KEGG_Pathway == pw_name, Pathway_KEGG_ID == pw_id) %>%
            distinct(Human_Gene_Symbol, .keep_all = TRUE) %>%
            arrange(Human_Gene_Symbol)
          
          gene_tags <- lapply(seq_len(nrow(gene_df)), function(j) {
            g <- gene_df$Human_Gene_Symbol[j]
            rid <- gene_df$record_id[j]
            
            tags$a(
              href = "javascript:void(0);",
              class = "gene-chip",
              onclick = sprintf(
                "Shiny.setInputValue('kegg_gene_click', {record_id: %s, nonce: Math.random()}, {priority: 'event'})",
                rid
              ),
              g
            )
          })
          
          pw_btn_id <- paste0("pathway_map_btn_", safe_id(cat_name, subcat_name, pw_id))
          
          tags$details(
            class = "subcategory-details",
            tags$summary(
              tags$span(
                class = "subcat-title",
                paste0(pw_name, " (", pw_id, ")")
              ),
              tags$a(
                href = "javascript:void(0);",
                class = "btn btn-default btn-xs",
                style = "margin-left:10px; margin-right:10px; padding:2px 10px;",
                onclick = sprintf(
                  "Shiny.setInputValue('pathway_map_pick', {pathway_id: '%s', pathway_name: '%s', nonce: Math.random()}, {priority: 'event'})",
                  js_escape(pw_id),
                  js_escape(pw_name)
                ),
                "Map"
              ),
              tags$span(class = "subcat-toggle-text")
            ),
            tags$div(
              class = "subcategory-content",
              gene_tags
            )
          )
        })
        
        tags$details(
          class = "subcategory-details",
          tags$summary(
            tags$span(class = "subcat-title", subcat_name),
            tags$span(class = "subcat-toggle-text")
          ),
          tags$div(
            class = "subcategory-content",
            pathway_tags
          )
        )
      })
      
      tags$div(
        class = "category-box",
        tags$div(
          class = "category-header-row",
          tags$h4(class = "top-title", cat_name)
        ),
        tags$details(
          class = "category-details",
          tags$summary(
            tags$span(class = "toggle-label")
          ),
          tags$div(
            class = "category-content",
            sub_tags
          )
        )
      )
    })
    
    do.call(tagList, pathway_panels)
  })
  
  
  
  
  
  selected_disease <- reactive({
    req(rv$selected_record_id)
    imd_all %>% filter(record_id == rv$selected_record_id)
  })
  
  selected_disease_pathways <- reactive({
    req(rv$selected_record_id)
    
    dx <- selected_disease()
    gene_symbol <- clean_text(dx$Gene[[1]])
    
    if (gene_symbol == "") {
      return(tibble())
    }
    
    get_gene_pathways(gene_symbol)
  })
  
  reset_selected_pathway <- function() {
    rv$selected_pathway_name <- NULL
    rv$selected_pathway_id <- NULL
    rv$selected_pathway_model <- NULL
    rv$selected_click_box <- NULL
    rv$selected_click_genes <- NULL
    rv$highlight_gene <- NULL
  }
  
  observeEvent(input$model_logic_help, {
    showModal(
      modalDialog(
        title = "Model Selection: Working Logic and Matrix",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
        tags$div(
          class = "help-section-title",
          "Working Logic"
        ),
        tags$img(
          src = "assets/Working_Logic.png",
          style = "width: 100%; max-width: 100%; margin-bottom: 12px;"
        ),
        
        tags$div(
          class = "logic-note",
          tags$p(tags$b("How the logic works")),
          tags$ul(
            tags$li("1:1 orthology is prioritized because it provides the clearest one-to-one correspondence for functional modelling."),
            tags$li("1:2 orthology is retained as usable, but interpretation should consider potential redundancy among paralogues."),
            tags$li("2:1 and more complex orthology relationships are deprioritized because they reduce modelling specificity and complicate functional interpretation."),
            tags$li("When multiple models remain suitable after orthology-based prioritization, the model selection matrix is used to rank them according to fitness for purpose.")
          )
        ),
        
        tags$div(
          class = "help-section-title",
          "Model Choice Matrix"
        ),
        tags$img(
          src = "assets/Model_choice_matrix.png",
          style = "width: 100%; max-width: 100%;"
        )
      )
    )
  })
  
  observeEvent(input$diopt_threshold, {
    if (!is.null(rv$selected_golden_model)) {
      
      if (grepl("\\| all ortholog support$", rv$selected_golden_model)) {
        model_name <- sub(" \\| all ortholog support$", "", rv$selected_golden_model)
        rv$search_results <- get_ortholog_model_diseases(model_name)
      } else {
        rv$search_results <- get_golden_model_diseases(rv$selected_golden_model)
      }
    }
    
    if (!is.null(rv$selected_category) && !is.null(rv$selected_model)) {
      rv$search_results <- get_category_model_diseases(
        rv$selected_category,
        rv$selected_model
      )
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$pathway_gene_click, {
    req(rv$selected_pathway_id)
    
    clicked_box <- clean_text(input$pathway_gene_click$box_id)
    if (clicked_box == "") return(NULL)
    
    rv$selected_category <- NULL
    rv$selected_model <- NULL
    rv$selected_click_box <- clicked_box
    
    clickmap_model_name <- get_clickmap_model_name(rv$selected_pathway_model)
    
    genes <- get_box_genes(
      model_name = clickmap_model_name,
      pathway_id = rv$selected_pathway_id,
      box_id = clicked_box
    )
    
    rv$selected_click_genes <- paste(genes, collapse = ", ")
    
    rv$search_results <- get_box_disease_records(
      model_name = clickmap_model_name,
      pathway_id = rv$selected_pathway_id,
      box_id = clicked_box
    )
    
    rv$search_results_visible <- TRUE
    rv$content_mode <- "results"
    rv$main_tab_selected <- "Search Results"
    updateTabsetPanel(session, "main_tabs", selected = "Search Results")
  }, ignoreInit = TRUE)
  
  
  output$selected_info <- renderUI({
    
    done <- start_timer("output$selected_info")
    on.exit(done(), add = TRUE)
    
    if (!identical(input$main_tabs, "Disease Page")) {
      return(NULL)
    }
    
    if (is.null(rv$selected_record_id)) {
      return(tags$div(class = "small-muted", "No disease selected yet."))
    }
    
    dx <- selected_disease()
    
    tags$div(
      tags$b(dx$ICIEM_Name[[1]]),
      tags$br(),
      paste0("Gene: ", dx$Gene[[1]]),
      tags$br(),
      paste0("Category: ", dx$Disease_Category[[1]])
    )
  })
  
  output$pathway_map_page <- renderUI({
    
    done <- start_timer("output$pathway_map_page")
    on.exit(done(), add = TRUE)
    
    if (!identical(input$main_tabs, "Pathway Map")) {
      return(NULL)
    }
    
    done <- start_timer("output$pathway_map_page")
    on.exit(done(), add = TRUE)
    
    done <- start_timer("output$pathway_map_page")
    on.exit(done(), add = TRUE)
    # -------------------------------------------------
    # 0) no disease selected and no pathway selected
    # -------------------------------------------------
    if (is.null(rv$selected_record_id) && (is.null(rv$selected_pathway_id) || is.null(rv$selected_pathway_name))) {
      return(tags$p("Please select a disease first, then choose one of its related pathways here."))
    }
    
    # -------------------------------------------------
    # 1) disease selected, but pathway not chosen yet
    # -------------------------------------------------
    if (!is.null(rv$selected_record_id) && (is.null(rv$selected_pathway_id) || is.null(rv$selected_pathway_name))) {
      dx <- selected_disease()
      gene_symbol <- clean_text(dx$Gene[[1]])
      pw <- selected_disease_pathways()
      
      if (nrow(pw) == 0) {
        return(tags$p(paste("No KEGG pathway map was found for gene", gene_symbol)))
      }
      
      pathway_buttons <- lapply(seq_len(nrow(pw)), function(i) {
        pw_name <- pw$KEGG_Pathway[i]
        pw_id   <- pw$Pathway_KEGG_ID[i]
        
        tags$a(
          href = "javascript:void(0);",
          class = "btn btn-default",
          style = "display:block; margin-bottom:12px; width:100%; max-width:900px; text-align:left; font-size:18px; padding:12px 18px;",
          onclick = sprintf(
            "Shiny.setInputValue('disease_pathway_pick', {pathway_id: '%s', pathway_name: '%s', nonce: Math.random()}, {priority: 'event'})",
            gsub("'", "\\\\'", pw_id),
            gsub("'", "\\\\'", pw_name)
          ),
          paste0(pw_name, " (", pw_id, ")")
        )
      })
      
      return(
        tagList(
          tags$div(
            class = "pathway-map-title",
            paste0("Pathways associated with ", gene_symbol)
          ),
          tags$div(
            class = "pathway-map-helper-title",
            paste0("Selected gene: ", gene_symbol)
          ),
          tags$div(
            class = "pathway-map-helper-text",
            paste0(
              "The pathways below include ", gene_symbol,
              ". Choose one pathway to explore how this gene is represented across model organisms and view the corresponding pathway map."
            )
          ),
          tags$div(
            style = "margin-top:8px; max-width:950px;",
            pathway_buttons
          )
        )
      )
    }
    
    # -------------------------------------------------
    # 2) pathway selected, but model not chosen yet
    # -------------------------------------------------
    if (is.null(rv$selected_pathway_model)) {
      model_cards <- lapply(model_display_order, function(model_name) {
        fig_src <- get_model_figure_rel(model_name)
        
        tags$div(
          class = "pathway-model-card",
          
          tags$a(
            href = "javascript:void(0);",
            onclick = sprintf(
              "Shiny.setInputValue('pathway_model_pick', {model: '%s', nonce: Math.random()}, {priority: 'event'})",
              js_escape(model_name)
            ),
            if (!is.null(fig_src)) {
              tags$img(src = fig_src, class = "pathway-model-thumb")
            } else {
              tags$div(
                class = "pathway-model-thumb",
                style = "display:flex;align-items:center;justify-content:center;color:#777;",
                "Image not found"
              )
            }
          ),
          
          tags$a(
            href = "javascript:void(0);",
            onclick = sprintf(
              "Shiny.setInputValue('pathway_model_pick', {model: '%s', nonce: Math.random()}, {priority: 'event'})",
              js_escape(model_name)
            ),
            class = "pathway-model-name",
            model_name
          )
        )
      })
      
      return(
        tagList(
          tags$div(
            style = "display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:10px;",
            
            tags$div(
              class = "pathway-map-title",
              "Explore pathway coverage across model organisms"
            ),
            
            actionButton(
              "back_to_pathway_list",
              "Back to pathway",
              class = "btn btn-default btn-sm",
              style = "background:#eef3f8; border-color:#cfd8e3; color:#4a5a6a; font-weight:600;"
            )
          ),
          
          tags$div(
            class = "pathway-map-helper-title",
            paste0(rv$selected_pathway_name, " (", rv$selected_pathway_id, ")")
          ),
          
          tags$div(
            class = "pathway-model-helper-text",
            tags$ul(
              tags$li(
                "Review all IMD-related genes mapped to this pathway and the biochemical processes in which they are involved."
              ),
              tags$li(
                "Select a model organism below to assess pathway-level orthology support and the coverage of IMD genes in that organism."
              )
            )
          ),
          
          tags$div(
            class = "pathway-model-grid",
            model_cards
          )
        )
      )
    }
    
    # -------------------------------------------------
    # 3) pathway + model selected: show detailed map
    # -------------------------------------------------
    map_src <- get_model_pathway_map_rel(rv$selected_pathway_model, rv$selected_pathway_id)
    
    clickmap_model_name <- get_clickmap_model_name(rv$selected_pathway_model)
    
    clickmap_tbl <- get_pathway_gene_clickmap()
    
    local_boxes_raw <- clickmap_tbl %>%
      filter(
        Model == clickmap_model_name,
        Pathway_KEGG_ID == rv$selected_pathway_id
      )
    
    if (nrow(local_boxes_raw) == 0) {
      local_boxes_raw <- clickmap_tbl %>%
        filter(
          Model == "Mus musculus",
          Pathway_KEGG_ID == rv$selected_pathway_id
        )
    }
    
    local_boxes <- local_boxes_raw %>%
      mutate(
        Box_ID = paste(Pathway_KEGG_ID, x1, y1, x2, y2, sep = "__")
      ) %>%
      group_by(
        Pathway_KEGG_ID, Box_ID,
        x1, y1, x2, y2, image_width, image_height
      ) %>%
      summarise(
        Genes = paste(sort(unique(Gene)), collapse = ", "),
        Gene_List = list(sort(unique(Gene))),
        .groups = "drop"
      ) %>%
      mutate(
        left_pct = 100 * x1 / image_width,
        top_pct = 100 * y1 / image_height,
        width_pct = 100 * (x2 - x1) / image_width,
        height_pct = 100 * (y2 - y1) / image_height
      )
    
    highlight_boxes <- tibble()
    
    if (!is.null(rv$highlight_gene) && nzchar(clean_text(rv$highlight_gene))) {
      highlight_boxes <- local_boxes %>%
        rowwise() %>%
        filter(rv$highlight_gene %in% Gene_List) %>%
        ungroup()
    }
    
    tagList(
      tags$div(
        style = "display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:12px;",
        
        tags$div(
          tags$div(
            class = "pathway-map-title",
            paste0(rv$selected_pathway_name, " (", rv$selected_pathway_id, ")")
          ),
          tags$div(
            class = "pathway-map-subtitle",
            rv$selected_pathway_model
          ),
          if (!is.null(rv$highlight_gene) && nzchar(clean_text(rv$highlight_gene))) {
            tags$div(
              style = "font-size:16px; color:#c0392b; margin-bottom:8px;",
              paste("Highlighted gene:", rv$highlight_gene)
            )
          }
        ),
        
        tags$div(
          class = "pathway-header-right",
          
          tags$div(
            style = "display:flex; gap:8px;",
            actionButton(
              "back_to_pathway_models",
              "Back to models",
              class = "btn btn-default btn-sm"
            ),
            actionButton(
              "back_to_pathway_list",
              "Back to pathway",
              class = "btn btn-default btn-sm"
            )
          ),
          
          tags$img(
            src = "assets/Legend.png",
            class = "pathway-legend"
          )
        )
      ),
      
      if (!is.null(map_src)) {
        tags$div(
          class = "pathway-map-container",
          tags$img(
            src = map_src,
            class = "pathway-map-image",
            style = "display:block; width:100%; height:auto;"
          ),
          lapply(seq_len(nrow(local_boxes)), function(i) {
            b <- local_boxes[i, ]
            
            tags$a(
              href = "javascript:void(0);",
              class = "pathway-ec-hotspot",
              title = paste("Gene(s):", b$Genes),
              onclick = sprintf(
                "Shiny.setInputValue('pathway_gene_click', {box_id: '%s', nonce: Math.random()}, {priority: 'event'})",
                gsub("'", "\\\\'", b$Box_ID)
              ),
              style = paste0(
                "left:", b$left_pct, "%;",
                "top:", b$top_pct, "%;",
                "width:", b$width_pct, "%;",
                "height:", b$height_pct, "%;"
              )
            )
          }),
          lapply(seq_len(nrow(highlight_boxes)), function(i) {
            b <- highlight_boxes[i, ]
            
            tags$div(
              class = "pathway-gene-highlight",
              title = paste("Highlighted gene in box:", rv$highlight_gene),
              style = paste0(
                "left:", b$left_pct, "%;",
                "top:", b$top_pct, "%;",
                "width:", b$width_pct, "%;",
                "height:", b$height_pct, "%;"
              )
            )
          })
        )
      } else {
        tags$p("No pathway map was found for this model organism and pathway.")
      }
    )
  })
  
  make_disease_section_header <- function(title, help_ui = NULL, toggle_id, collapsed,
                                          margin_top = "16px", margin_bottom = "8px") {
    tags$div(
      style = paste0(
        "display:flex;",
        "align-items:center;",
        "justify-content:space-between;",
        "gap:12px;",
        "margin-top:", margin_top, ";",
        "margin-bottom:", margin_bottom, ";"
      ),
      
      tags$div(
        style = "display:flex; align-items:center; gap:8px;",
        tags$div(class = "section-title", title),
        help_ui
      ),
      
      actionButton(
        toggle_id,
        if (isTRUE(collapsed)) "Expand" else "Collapse",
        class = "btn btn-default btn-sm golden-collapse-btn"
      )
    )
  }
  
  output$disease_page <- renderUI({
    
    done <- start_timer("output$disease_page")
    on.exit(done(), add = TRUE)
    
    if (!identical(input$main_tabs, "Disease Page")) {
      return(NULL)
    }
    
    done <- start_timer("output$disease_page")
    on.exit(done(), add = TRUE)
    
    
    done <- start_timer("output$disease_page")
    on.exit(done(), add = TRUE)
    
    if (is.null(rv$selected_record_id)) {
      return(tags$p("Please select a disease from the Search Results tab or browse by category to continue."))
    }
    
    dx <- selected_disease()
    md <- get_disease_models(dx, orth_all_current())
    fmt_model_list_html <- function(x) {
      x <- clean_text(x)
      x <- x[x != ""]
      
      if (length(x) == 0) {
        return("none")
      }
      
      tagList(
        lapply(seq_along(x), function(i) {
          tagList(
            tags$span(class = "model-name-highlight", x[i]),
            if (i < length(x)) ", "
          )
        })
      )
    }
    
    ortholog_models <- md %>%
      filter(Has_Ortholog) %>%
      pull(Model)
    
    one_to_one_models <- md %>%
      filter(clean_text(Orthology_Relationship) == "1:1") %>%
      pull(Model)
    
    one_to_two_models <- md %>%
      filter(clean_text(Orthology_Relationship) == "1:2") %>%
      pull(Model)
    
    complex_models <- md %>%
      filter(
        Has_Ortholog,
        !clean_text(Orthology_Relationship) %in% c("1:1", "1:2")
      ) %>%
      pull(Model)
    
    model_summary_ui <- if (length(ortholog_models) == 0) {
      tags$p(paste0("No ortholog support was identified across the ", nrow(md), " model organisms."))
    } else {
      tags$p(
        "This gene has ortholog support in ",
        tags$b(length(ortholog_models)),
        " of ",
        tags$b(nrow(md)),
        " model organisms. ",
        
        if (length(one_to_one_models) > 0) {
          tagList(
            "Models with 1:1 orthology: ",
            fmt_model_list_html(one_to_one_models),
            ". "
          )
        } else {
          "No model shows 1:1 orthology. "
        },
        
        if (length(one_to_two_models) > 0) {
          tagList(
            "Models with 1:2 orthology: ",
            fmt_model_list_html(one_to_two_models),
            ". "
          )
        } else {
          NULL
        },
        
        if (length(complex_models) > 0) {
          tagList(
            "Models with complex orthology: ",
            fmt_model_list_html(complex_models),
            "."
          )
        } else {
          "No model shows complex orthology."
        }
      )
    }
    research_df <- get_research_in_models(dx$Gene[[1]])
    is_cd_gene <- is_childhood_dementia_gene(dx$Gene[[1]])
    
    info_tbl <- tags$table(
      class = "table table-bordered",
      tags$tbody(
        tags$tr(tags$th("ICIEM Name"), tags$td(dx$ICIEM_Name[[1]])),
        tags$tr(
          tags$th("Gene"),
          tags$td(
            actionLink(
              inputId = "gene_to_pathway",
              label = dx$Gene[[1]],
              class = "gene-pathway-link",
              title = "Click to open related pathway map"
            ),
            if (is_cd_gene) {
              tags$span(
                class = "childhood-dementia-mark",
                title = "Related to Childhood Dementia",
                "*"
              )
            }
          )
        ),
        tags$tr(tags$th("Disease Category"), tags$td(dx$Disease_Category[[1]])),
        tags$tr(tags$th("Disease Subcategory"), tags$td(dx$Disease_Subcategory[[1]])),
        tags$tr(tags$th("Alternative Names"), tags$td(dx$Alternative_Names[[1]])),
        tags$tr(tags$th("Disease Abbreviation"), tags$td(dx$Disease_Abbreviation[[1]])),
        tags$tr(tags$th("Mode of Inheritance"), tags$td(dx$Mode_of_Inheritance[[1]])),
        tags$tr(tags$th("Treatable?"), tags$td(dx$Treatable[[1]])),
        tags$tr(tags$th("Prevalence"), tags$td(dx$Prevalence[[1]])),
        tags$tr(tags$th("OMIM"), tags$td(dx$OMIM[[1]])),
        tags$tr(tags$th("OMIM Link"), tags$td(HTML(make_clickable_link(dx$OMIM_Link[[1]], "Open")))),
        tags$tr(tags$th("OrphaCode"), tags$td(dx$OrphaCode[[1]])),
        tags$tr(tags$th("OrphaCode Link"), tags$td(HTML(make_clickable_link(dx$OrphaCode_Link[[1]], "Open")))),
        tags$tr(tags$th("IEMBASE Link"), tags$td(HTML(make_clickable_link(dx$IEMBASE_Link[[1]], "Open")))),
        tags$tr(tags$th("GeneReviews"), tags$td(dx$GeneReviews[[1]])),
        tags$tr(tags$th("GeneReviews Link"), tags$td(HTML(make_clickable_link(dx$GeneReviews_Link[[1]], "Open"))))
      )
    )
    
    research_section <- tagList(
      make_disease_section_header(
        title = "Published IMD studies in model organisms",
        help_ui = tags$span(
          class = "help-btn",
          title = "adapted from Alliance of Genome Resources",
          "?"
        ),
        toggle_id = "toggle_research_section",
        collapsed = rv$research_section_collapsed,
        margin_top = "28px",
        margin_bottom = "28px"
      ),
      
      if (!isTRUE(rv$research_section_collapsed)) {
        if (is.null(research_df) || nrow(research_df) == 0) {
          tags$div(class = "research-unavailable", "unavailable")
        } else if (all(research_df$Status == "unavailable")) {
          tags$div(class = "research-unavailable", "unavailable")
        } else {
          species_order <- unique(research_df$Species[research_df$Status == "available"])
          
          species_blocks <- lapply(seq_along(species_order), function(i) {
            sp <- species_order[i]
            sp_df <- research_df %>%
              dplyr::filter(Species == sp, Status == "available")
            
            tags$div(
              class = "research-species-block",
              tags$div(
                class = "research-species-title",
                paste0("Species ", i, ": ", sp)
              ),
              lapply(seq_len(nrow(sp_df)), function(j) {
                tags$div(
                  class = "research-paper-block",
                  tags$div(
                    class = "research-paper-title",
                    HTML(paste0("<b>Title:</b> ", htmlEscape(sp_df$Title[j])))
                  ),
                  tags$div(
                    class = "research-paper-ref",
                    HTML(paste0("<b>References:</b> ", make_pubmed_link(sp_df$Reference[j])))
                  )
                )
              })
            )
          })
          
          tagList(species_blocks)
        }
      }
    )
    
    
    model_cards <- lapply(seq_len(nrow(md)), function(i) {
      row <- md[i, ]
      
      bg <- case_when(
        row$Priority == "Prioritized" ~ "#e8f5e9",   # green for 1:1
        row$Priority == "Retained" ~ "#fff8e1",      # yellow for 1:2
        row$Priority == "Deprioritized" ~ "#f8d7da", # red/pink for everything else
        TRUE ~ "#f5f5f5"
      )
      
      model_label <- row$Model
      if (model_label == "Worm") model_label <- "Worms"
      if (model_label == "Fly") model_label <- "Flies"
      if (model_label == "Mouse") model_label <- "Mice"
      if (model_label == "Zebrafish") model_label <- "Zebrafish"
      if (model_label == "Fission Yeast") model_label <- "Fission Yeast"
      if (model_label == "Budding Yeast") model_label <- "Budding Yeast"
      
      has_support <- isTRUE(row$Has_Ortholog)
      
      tags$div(
        class = "model-card",
        style = paste0("background:", bg, ";"),
        tags$h4(row$Model),
        
        if (has_support && clean_text(row$Model_DB_Link) != "") {
          tags$p(
            tags$b("Model database: "),
            HTML(make_clickable_link(row$Model_DB_Link, "Open"))
          )
        },
        
        if (has_support) {
          tagList(
            tags$p(tags$b("Orthology relationship: "), row$Orthology_Relationship),
            tags$p(
              tags$b("Humans in component: "),
              row$Human_Component_With_DIOPT
            ),
            tags$p(
              tags$b(paste0(model_label, " in component: ")),
              row$Model_Component_With_DIOPT
            )
          )
        } else {
          tags$p(
            class = "small-muted",
            "No ortholog support identified for this gene in this model."
          )
        }
      )
    })
    
    tagList(
      make_disease_section_header(
        title = dx$ICIEM_Name[[1]],
        help_ui = tags$span(
          class = "help-btn",
          title = "Disease information is adapted from IEMBase",
          "?"
        ),
        toggle_id = "toggle_disease_info_section",
        collapsed = rv$disease_info_collapsed,
        margin_top = "0px",
        margin_bottom = "12px"
      ),
      
      if (!isTRUE(rv$disease_info_collapsed)) {
        info_tbl
      },
      
      tags$div(class = "section-divider"),
      
      research_section,
      
      tags$div(class = "section-divider"),
      
      make_disease_section_header(
        title = "Model Selection Output",
        help_ui = actionButton("model_logic_help", "?", class = "help-btn"),
        toggle_id = "toggle_model_section",
        collapsed = rv$model_section_collapsed,
        margin_top = "16px",
        margin_bottom = "8px"
      ),
      
      if (!isTRUE(rv$model_section_collapsed)) {
        tagList(
          tags$div(
            style = paste0(
              "font-size:14px;",
              "color:#4f5f73;",
              "line-height:1.5;",
              "margin-top:2px;",
              "margin-bottom:10px;"
            ),
            tags$p(
              style = "margin-bottom:4px;",
              tags$b("Current DIOPT threshold used: "),
              paste0("≥ ", input$diopt_threshold)
            ),
            tags$p(
              style = "margin-bottom:0;",
              "For each gene, the DIOPT score indicates how many orthology prediction tools support its assignment to the shown orthology relationship."
            )
          ),
          
          model_summary_ui,
          model_cards
        )
      }
    )
  })
}

shinyApp(ui, server)









