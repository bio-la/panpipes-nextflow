#!/usr/bin/env Rscript

# Code adapted for Nextflow - M. Gonzales 
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(magrittr)
  library(optparse)
  library(jsonlite)
  library(data.table)
})

options(stringsAsFactors = FALSE)
options(bitmaptype = "cairo")

# -------------------------- helpers --------------------------

scale_y_origin <- function(...) {
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), ...)
}

panpipes_theme <- function(font_size = 10) {
  theme(
    text = element_text(size = font_size),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "grey", linetype = 3, size = 0.5),
    axis.text = element_text(size = font_size),
    axis.title = element_text(size = font_size),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )
}

# Turn a list/vector of variables into a data.frame with modality + variable
# Accepts:
#   - character vector like c("rna:leiden_res1","sample_id")
#   - named list like list(rna=c("rna:leiden_res1"), prot=c("prot:leiden_res1"))
parse_vars <- function(x) {
  if (is.null(x)) return(NULL)

  if (is.list(x)) {
    uniq <- unique(unlist(x))
  } else {
    uniq <- unique(as.character(x))
  }
  if (length(uniq) == 0) return(NULL)

  # Split "mod:var" and normalize column names to use "_" instead of ":"
  parts <- str_split(uniq, ":", simplify = TRUE)
  mod <- ifelse(parts[,1] %in% c("rna","prot","atac","rep"), parts[,1], "multimodal")
  tibble(
    mod = mod,
    variable = gsub(":", "_", uniq)
  )
}

# Ensure factors for listed variables; silently skip missing cols
parse_cell_metadata <- function(cmtd, cat_vars, grp_vars) {
  colnames(cmtd) <- gsub(":", "_", colnames(cmtd))
  if ("rep_has_ir" %in% colnames(cmtd)) {
    cmtd$rep_has_ir <- tidyr::replace_na(cmtd$rep_has_ir, FALSE)
  }
  if (!is.null(cat_vars) && nrow(cat_vars) > 0) {
    vars <- intersect(cat_vars$variable, colnames(cmtd))
    if (length(vars)) cmtd <- cmtd %>% mutate(across(all_of(vars), as.factor))
  }
  if (!is.null(grp_vars) && nrow(grp_vars) > 0) {
    vars <- intersect(grp_vars$variable, colnames(cmtd))
    if (length(vars)) cmtd <- cmtd %>% mutate(across(all_of(vars), as.factor))
  }
  cmtd
}

# -------------------------- CLI --------------------------

option_list <- list(
  make_option(c("--mtd_object"), help = "TSV with per-cell metadata (mdata.obs)"),
  make_option(c("--grouping_vars_json"), default = NULL,
              help = "JSON: list or vector of grouping vars (e.g. '{\"rna\":[\"rna:leiden_res1\"],\"all\":[\"sample_id\"]}' or '[\"sample_id\"]')"),
  make_option(c("--categorical_vars_json"), default = NULL,
              help = "JSON: list or vector of categorical vars"),
  make_option(c("--continuous_vars_json"), default = NULL,
              help = "JSON: list or vector of continuous vars"),
  make_option(c("--do_categorical_barplots"), default = "true", help = "true/false"),
  make_option(c("--do_categorical_stacked_barplots"), default = "true", help = "true/false"),
  make_option(c("--do_continuous_violin"), default = "true", help = "true/false"),
  make_option(c("--outdir"), default = ".", help = "Base output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))
message("Running with options:")
print(opt)

# -------------------------- IO --------------------------

if (is.null(opt$mtd_object)) {
  stop("Missing --mtd_object")
}
cmtd <- data.table::fread(opt$mtd_object, sep = "\t", na.strings = "")

# Parse JSON helpers (empty/NULL-safe)
parse_json <- function(s) {
  if (is.null(s)) return(NULL)
  s <- trimws(s)
  if (s == "" || s == "{}" || s == "[]" || tolower(s) %in% c("null","none")) return(NULL)
  jsonlite::fromJSON(s, simplifyVector = TRUE)
}

grouping_raw   <- parse_json(opt$grouping_vars_json)
categorical_raw<- parse_json(opt$categorical_vars_json)
continuous_raw <- parse_json(opt$continuous_vars_json)

# If 'all' key exists, merge into each modality
merge_all <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || is.null(x$all)) return(x)
  allv <- x$all
  x[["all"]] <- NULL
  lapply(x, function(v) unique(c(v, allv)))
}
grouping_raw    <- merge_all(grouping_raw)
categorical_raw <- merge_all(categorical_raw)
continuous_raw  <- merge_all(continuous_raw)

grp_vars <- parse_vars(grouping_raw)
cat_vars <- parse_vars(categorical_raw)
cont_vars<- parse_vars(continuous_raw)

cmtd <- parse_cell_metadata(cmtd, cat_vars, grp_vars)

# choose active modalities present in the tables
modalities_cat <- if (!is.null(cat_vars)) unique(cat_vars$mod) else character(0)
modalities_grp <- if (!is.null(grp_vars)) unique(grp_vars$mod) else character(0)
modalities_cont<- if (!is.null(cont_vars)) unique(cont_vars$mod) else character(0)

# Coerce booleans
to_bool <- function(x) tolower(as.character(x)) %in% c("true","1","yes","y")
DO_BAR       <- to_bool(opt$do_categorical_barplots)
DO_STACKED   <- to_bool(opt$do_categorical_stacked_barplots)
DO_VIOLIN    <- to_bool(opt$do_continuous_violin)

# Ensure outdir exists
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)
setwd(opt$outdir)

# -------------------------- Plots --------------------------

# Barplots per categorical variable
if (DO_BAR && !is.null(cat_vars) && nrow(cat_vars) > 0) {
  for (mod in unique(modalities_cat)) {
    if (!dir.exists(file.path(mod))) dir.create(file.path(mod))
    message(sprintf("Plotting barplots for modality %s", mod))
    uniq_cat_vars <- cat_vars %>% filter(mod == !!mod) %>% pull(variable)
    for (cat_v in uniq_cat_vars) {
      if (!cat_v %in% colnames(cmtd)) next
      plot_dat <- cmtd %>%
        mutate(!!cat_v := factor(!!rlang::sym(cat_v))) %>%
        group_by(!!sym(cat_v)) %>%
        summarise(n_cells = n(), .groups = "drop") %>%
        drop_na()
      if (nrow(plot_dat) == 0) next
      flip_axes <- max(nchar(as.character(plot_dat[[cat_v]]))) > 10
      p1 <- ggplot(plot_dat, aes_string(x = cat_v, y = "n_cells", fill = cat_v)) +
        geom_col() +
        ggtitle(cat_v) +
        scale_y_origin() +
        panpipes_theme() +
        theme(legend.position = "none")
      if (flip_axes) p1 <- p1 + coord_flip()
      outfile <- file.path(mod, paste0("bar_", cat_v, ".png"))
      message("Saving: ", outfile)
      cowplot::save_plot(filename = outfile, plot = p1, ncol = 1, nrow = 1,
                         base_asp = 1.3, base_height = 5)
    }
  }
}

# Stacked barplots by grouping var
if (DO_STACKED && !is.null(cat_vars) && nrow(cat_vars) > 0 && !is.null(grp_vars) && nrow(grp_vars) > 0) {
  for (mod in unique(modalities_cat)) {
    message(sprintf("Plotting stacked barplots for modality %s", mod))
    uniq_cat_vars <- cat_vars %>% filter(mod == !!mod) %>% pull(variable)
    for (gv in grp_vars$variable) {
      if (!gv %in% colnames(cmtd)) next
      dir.create(file.path(mod, gv), recursive = TRUE, showWarnings = FALSE)
      for (cat_v in uniq_cat_vars) {
        if (gv == cat_v || !cat_v %in% colnames(cmtd)) next
        plot_dat <- cmtd %>%
          select(all_of(c(gv, cat_v))) %>%
          drop_na() %>%
          group_by(!!sym(gv), !!sym(cat_v)) %>%
          summarise(n_cells = n(), .groups = "drop_last") %>%
          group_by(!!sym(gv)) %>%
          mutate(proportion = n_cells / sum(n_cells)) %>%
          ungroup()
        if (nrow(plot_dat) == 0) next
        flip_axes <- max(nchar(as.character(plot_dat[[gv]]))) > 10
        p1 <- ggplot(plot_dat, aes_string(x = gv, y = "n_cells", fill = cat_v)) +
          geom_col() + scale_y_origin() + panpipes_theme(font_size = 10)
        p2 <- ggplot(plot_dat, aes_string(x = gv, y = "proportion", fill = cat_v)) +
          geom_col() + scale_y_origin() + panpipes_theme(font_size = 10)
        if (flip_axes) { p1 <- p1 + coord_flip(); p2 <- p2 + coord_flip() }
        pg <- cowplot::plot_grid(p1, p2)
        outfile <- file.path(mod, gv, paste0("stackedbar_", cat_v, ".png"))
        message("Saving: ", outfile)
        cowplot::save_plot(filename = outfile, plot = pg, ncol = 2, nrow = 1,
                           base_height = 5, base_asp = ifelse(flip_axes, 1.5, 1.1))
      }
    }
  }
}

# Violin plots for continuous variables grouped by grouping vars
if (DO_VIOLIN && !is.null(cont_vars) && nrow(cont_vars) > 0 && !is.null(grp_vars) && nrow(grp_vars) > 0) {
  for (mod in unique(modalities_cont)) {
    message(sprintf("Plotting violin plots for modality %s", mod))
    uniq_cont_vars <- cont_vars %>% filter(mod == !!mod) %>% pull(variable)
    for (gv in grp_vars$variable) {
      if (!gv %in% colnames(cmtd)) next
      dir.create(file.path(mod, gv), recursive = TRUE, showWarnings = FALSE)
      for (cont_v in uniq_cont_vars) {
        if (gv == cont_v || !cont_v %in% colnames(cmtd)) next
        plot_dat <- cmtd %>% filter(!is.na(!!sym(cont_v)) & !is.na(!!sym(gv)))
        if (nrow(plot_dat) == 0) next
        flip_axes <- max(nchar(as.character(plot_dat[[gv]]))) > 10
        p1 <- ggplot(plot_dat, aes_string(x = gv, y = cont_v, fill = gv)) +
          geom_violin() + scale_y_origin() + panpipes_theme(font_size = 10)
        if (flip_axes) p1 <- p1 + coord_flip()
        outfile <- file.path(mod, gv, paste0("violin_", cont_v, ".png"))
        message("Saving: ", outfile)
        cowplot::save_plot(filename = outfile, plot = p1, nrow = 1, ncol = 1,
                           base_height = 5, base_asp = ifelse(flip_axes, 1.5, 1.1))
      }
    }
  }
}
