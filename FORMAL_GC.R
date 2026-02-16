# ===============================
# CALCULATING GC and ∆GC CONTENT
# ===============================

# This script is designed with bacterial hosts and bacteriophages in mind
# it can be applied to other virus/host pairs :)

# Set base and output directories
# Load metadata - should have columns "sample_id" (organism name), "file_path" (file path for the specific FASTA sequence),
#                 "sample_type" (i.e. phage or host), "host_genus", "host_species", "lifestyle" (Virulent or Temperate),
#                 and "trna_group" (High, Low, or None)

base_dir      <- "/your/file/path/FASTAs"   # assumes FASTA files are all saved to a single parent folder
metadata_xlsx <- "/your/file/path/metadata.xlsx" 
out_dir       <- "/your/file/path/Results/GC"

# Load required packages

req_pkgs <- c("readxl", "dplyr", "tidyr", "purrr", "stringr", "seqinr", "ggplot2")

not_installed <- req_pkgs[!req_pkgs %in% installed.packages()[,"Package"]]
if(length(not_installed)) install.packages(not_installed, dependencies = TRUE)
invisible(lapply(req_pkgs, library, character.only = TRUE))

# ==================================================
# Helper functions for calculating and comparing GC
# ==================================================

stop_if_missing_cols <- function(df, cols) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) stop("Metadata is missing required column(s): ", paste(miss, collapse = ", "))
}

resolve_path <- function(base_dir, file_path) {
  if (is.na(file_path) || !nzchar(file_path)) return(NA_character_)
  # treat as absolute if starts with "/" (mac/linux) or drive letter (windows)
  if (grepl("^/", file_path) || grepl("^[A-Za-z]:\\\\", file_path)) return(file_path)
  file.path(base_dir, file_path)
}

clean_sample_id_from_path <- function(path) {
  x <- tools::file_path_sans_ext(basename(path))
  x <- sub("_per_gene_GC$", "", x)
  x <- sub(" CDS$", "", x)
  x
}

# ==========================
# Calculation of GC Content
# ==========================

calculate_gc_content <- function(dna_seq) {
  c(
    Total_GC = seqinr::GC(dna_seq),
    GC1      = seqinr::GC1(dna_seq),
    GC2      = seqinr::GC2(dna_seq),
    GC3      = seqinr::GC3(dna_seq)
  )
}

calc_gc_from_fasta <- function(fasta_file, sample_id) {
  x <- seqinr::read.fasta(file = fasta_file, seqtype = "DNA")
  gc <- lapply(x, calculate_gc_content)
  df <- as.data.frame(do.call(rbind, gc))
  df$Gene      <- names(x)
  df$sample_id <- sample_id
  df <- df[, c("Gene", "sample_id", "Total_GC", "GC1", "GC2", "GC3")]
  df
}

# ===================================
# Functions for statistical analyses
# ===================================

# Function to compare GC of phage to GC of host
compare_phages_to_host <- function(host_gc_df, phage_gc_df) {
  phage_list <- unique(phage_gc_df$sample_id)
  out <- vector("list", length(phage_list) * 4)
  k <- 1
  for (phage in phage_list) {
    phage_data <- phage_gc_df %>% filter(sample_id == phage)
    for (gc_type in c("Total_GC", "GC1", "GC2", "GC3")) {
      if (nrow(phage_data) > 0 && !all(is.na(phage_data[[gc_type]]))) {
        p_value <- wilcox.test(host_gc_df[[gc_type]], phage_data[[gc_type]])$p.value
      } else {
        p_value <- NA_real_
      }
      out[[k]] <- data.frame(
        Phage   = phage,
        GC_Type = gc_type,
        p_value = p_value
      )
      k <- k + 1
    }
  }
  bind_rows(out)
}

# Find the mean values of GC content for each host and the collective phages
gc_means <- function(df) {
  df %>%
    summarise(
      Total_GC = mean(Total_GC, na.rm = TRUE),
      GC1      = mean(GC1,      na.rm = TRUE),
      GC2      = mean(GC2,      na.rm = TRUE),
      GC3      = mean(GC3,      na.rm = TRUE)
    )
}

# ============================================
# Load the metadata into R and standardize it
# ============================================

meta <- readxl::read_excel(metadata_xlsx)

# REQUIRED COLUMNS
stop_if_missing_cols(meta, c("file_path", "sample_type", "host_genus"))

# OPTIONAL but RECOMMENDED COLUMNS
if (!("sample_id" %in% colnames(meta))) meta$sample_id <- NA_character_
if (!("trna_group" %in% colnames(meta))) meta$trna_group <- NA_character_

# normalize key columns
meta <- meta %>%
  mutate(
    sample_type = tolower(sample_type),
    host_genus  = as.character(host_genus),
    trna_group  = as.character(trna_group),
    abs_path    = vapply(file_path, function(p) resolve_path(base_dir, p), character(1)),
    sample_id   = ifelse(is.na(sample_id) | !nzchar(sample_id),
                         vapply(abs_path, clean_sample_id_from_path, character(1)),
                         as.character(sample_id))
  ) %>%
  filter(!is.na(abs_path), nzchar(abs_path))

# ==========================================================
# Compute per-gene GC content for all organisms in metadata
# ==========================================================

gc_out_dir <- file.path(out_dir, "GC")
dir.create(gc_out_dir, recursive = TRUE, showWarnings = FALSE)

# Compute + write per-sample per-gene GC
all_per_gene <- list()
all_per_sample_means <- list()

for (i in seq_len(nrow(meta))) {
  row <- meta[i, ]
  path <- row$abs_path
  sid  <- row$sample_id
  if (!file.exists(path)) {
    warning("Missing file: ", path)
    next
  }
  message("GC: ", row$host_genus, " | ", row$sample_type, " | ", sid)
  df <- calc_gc_from_fasta(path, sample_id = sid) %>%
    mutate(
      host_genus  = row$host_genus,
      sample_type = row$sample_type,
      trna_group  = ifelse(row$sample_type == "phage", row$trna_group, NA_character_)
    )
  genus_dir <- file.path(gc_out_dir, row$host_genus)
  dir.create(genus_dir, recursive = TRUE, showWarnings = FALSE)
  out_per_gene <- file.path(genus_dir, paste0(sid, "_per_gene_GC.csv"))
  write.csv(df, out_per_gene, row.names = FALSE)
  means <- df %>%
    summarise(
      Total_GC = mean(Total_GC, na.rm = TRUE),
      GC1      = mean(GC1,      na.rm = TRUE),
      GC2      = mean(GC2,      na.rm = TRUE),
      GC3      = mean(GC3,      na.rm = TRUE)
    ) %>%
    mutate(
      sample_id   = sid,
      host_genus  = row$host_genus,
      sample_type = row$sample_type,
      trna_group  = ifelse(row$sample_type == "phage", row$trna_group, NA_character_)
    ) %>%
    select(sample_id, host_genus, sample_type, trna_group, Total_GC, GC1, GC2, GC3)
  all_per_gene[[length(all_per_gene) + 1]] <- df
  all_per_sample_means[[length(all_per_sample_means) + 1]] <- means
}

per_gene_gc <- bind_rows(all_per_gene)
per_sample_mean_gc <- bind_rows(all_per_sample_means)

write.csv(per_gene_gc, file.path(gc_out_dir, "ALL_samples_per_gene_GC.csv"), row.names = FALSE)
write.csv(per_sample_mean_gc, file.path(gc_out_dir, "ALL_samples_mean_GC.csv"), row.names = FALSE)

# ======================================
# Outputs and statistics per host genus
# ======================================

# Set your preferred genus order
genus_order <- c("Host1","Host2","...","HostN")

# Only analyze genera present
genera_present <- intersect(genus_order, unique(per_gene_gc$host_genus))
if (length(genera_present) == 0) genera_present <- sort(unique(per_gene_gc$host_genus))

for (g in genera_present) {
  genus_dir <- file.path(gc_out_dir, g)
  dir.create(genus_dir, recursive = TRUE, showWarnings = FALSE)
  host_df  <- per_gene_gc %>% filter(host_genus == g, sample_type == "host")
  phage_df <- per_gene_gc %>% filter(host_genus == g, sample_type == "phage")
  # Save combined per-gene files
  write.csv(host_df,  file.path(genus_dir, paste0(g, "_hosts_per_gene_GC.csv")), row.names = FALSE)
  write.csv(phage_df, file.path(genus_dir, paste0(g, "_phages_per_gene_GC.csv")), row.names = FALSE)
  # If either is missing, skip stats
  if (nrow(host_df) == 0 || nrow(phage_df) == 0) next
  # Host vs each phage
  ph_vs_host <- compare_phages_to_host(host_df, phage_df)
  write.csv(ph_vs_host, file.path(genus_dir, paste0(g, "_phage_vs_host_wilcox.csv")), row.names = FALSE)
}

# ============================
# Mean Summary per host genus
# ============================

genus_summary <- bind_rows(lapply(genera_present, function(g) {
  host_df  <- per_gene_gc %>% filter(host_genus == g, sample_type == "host")
  phage_df <- per_gene_gc %>% filter(host_genus == g, sample_type == "phage")
  bind_rows(
    if (nrow(host_df)  > 0) gc_means(host_df)  %>% mutate(Genus = g, Type = "Host")  else NULL,
    if (nrow(phage_df) > 0) gc_means(phage_df) %>% mutate(Genus = g, Type = "Phages") else NULL
  )
}))

write.csv(genus_summary, file.path(gc_out_dir, "genus_gc_means.csv"), row.names = FALSE)

# =======================================================
# Violin plot to compare host to phages, faceted by host
# =======================================================

# plots are formatted to have y-axis in percentages, no gridlines, and clean labels

library(scales)  

plot_df <- per_gene_gc %>%
  mutate(
    Genus = factor(host_genus, levels = genus_order),
    Set   = ifelse(sample_type == "host", "Host", "Phages"),
    Set   = factor(Set, levels = c("Host","Phages"))
  ) %>%
  filter(!is.na(Genus)) %>%
  pivot_longer(
    cols      = c("Total_GC","GC1","GC2","GC3"),
    names_to  = "GC_Type",
    values_to = "GC_fraction"
  ) %>%
  mutate(
    GC_Type = factor(GC_Type, levels = c("GC1","GC2","GC3","Total_GC")),
    GC_Type = gsub("_", " ", GC_Type),
    Genus   = factor(gsub("_", " ", as.character(Genus)), levels = gsub("_"," ", genus_order))
  )
p_gc <- ggplot(plot_df, aes(x = Set, y = GC_fraction, fill = Set)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, color = "black") +
  facet_grid(GC_Type ~ Genus) +
  labs(x = "", y = "GC content (%)") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "Host"   = "#F8766D",
    "Phages" = "#00B0F6"
  )) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    panel.spacing = unit(0.6, "lines"),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16)
  )
ggsave(
  file.path(gc_out_dir, "GC_host_vs_phage_by_genus_violin.png"),
  plot = p_gc, width = 18, height = 8, dpi = 300
)

print(p_gc)

# =========================================
# Calculation of ∆GC content (phage - host)
# =========================================

# =================
# Helper functions
# =================

standardize_text <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("_", " ") %>%  # "no underscores"
    stringr::str_squish()
}

require_columns <- function(df, cols, df_name = "df") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(df_name, " is missing required columns: ", paste(miss, collapse = ", "))
  }
  df
}

# ============================================
# Re-load and clean metadata for ∆GC analysis
# ============================================

meta <- readxl::read_excel(metadata_xlsx) %>%
  mutate(across(where(is.character), standardize_text)) %>%
  transmute(
    sample_id,
    sample_type,
    host_genus,
    host_species,
    lifestyle,
    trna_group
  ) %>%
  distinct(sample_id, .keep_all = TRUE)

# bind a column for host_species to the per_gene_gc object
require_columns(per_gene_gc, c("sample_id", "sample_type", "Total_GC", "GC1", "GC2", "GC3"), "per_gene_gc")

per_gene_gc2 <- per_gene_gc %>%
  mutate(across(where(is.character), standardize_text)) %>%
  # Drop metadata-like columns so they only come from `meta`
  select(-any_of(c("host_genus","host_species","lifestyle","trna_group"))) %>%
  left_join(
    meta %>% select(sample_id, host_genus, host_species, lifestyle, trna_group),
    by = "sample_id",
    relationship = "many-to-one"
  )

# ===============================================================
# Calculate ∆GC for each phage compared to its host (phage-host)
# ===============================================================

hosts_long <- per_gene_gc2 %>%
  filter(sample_type == "host", !is.na(host_species)) %>%
  select(host_species, host_genus, any_of(gc_types))

phages_long <- per_gene_gc2 %>%
  filter(sample_type == "phage") %>%
  select(sample_id, host_species, host_genus, lifestyle, trna_group, any_of(gc_types)) %>%
  filter(!is.na(sample_id), !is.na(host_species), !is.na(host_genus))

gc_types <- c("GC1","GC2","GC3","Total_GC")

delta_gc <- phages_long %>%
  group_by(sample_id, host_species, host_genus, lifestyle, trna_group) %>%
  group_modify(~{
    phage_df <- .x
    keys     <- .y
    hs <- keys$host_species[[1]]
    host_df <- hosts_long %>%
      dplyr::filter(host_species == hs)
    if (nrow(host_df) == 0) {
      return(tibble::tibble(
        GC_Type       = gc_types,
        Host_Mean_GC  = NA_real_,
        Phage_Mean_GC = sapply(gc_types, function(gt) mean(phage_df[[gt]], na.rm = TRUE)),
        Delta_GC      = NA_real_,
        p_value       = NA_real_,
        n_host_genes  = 0L,
        n_phage_genes = nrow(phage_df)
      ))
    }
    out <- lapply(gc_types, function(gt) {
      host_vals  <- host_df[[gt]]
      phage_vals <- phage_df[[gt]]
      host_mean  <- mean(host_vals,  na.rm = TRUE)
      phage_mean <- mean(phage_vals, na.rm = TRUE)
      delta      <- phage_mean - host_mean
      pval <- tryCatch(
        stats::wilcox.test(phage_vals, host_vals)$p.value,
        error = function(e) NA_real_
      )
      tibble::tibble(
        GC_Type       = gt,
        Host_Mean_GC  = host_mean,
        Phage_Mean_GC = phage_mean,
        Delta_GC      = delta,
        p_value       = pval,
        n_host_genes  = sum(!is.na(host_vals)),
        n_phage_genes = sum(!is.na(phage_vals))
      )
    })
    dplyr::bind_rows(out)
  }) %>%
  ungroup() %>%
  dplyr::rename(Phage = sample_id) %>%
  dplyr::mutate(
    host_genus = factor(host_genus, levels = genus_order),
    GC_Type    = factor(GC_Type, levels = c("GC1","GC2","GC3","Total_GC")),
    lifestyle  = factor(lifestyle, levels = c("Virulent","Temperate")),
    trna_group = factor(trna_group, levels = c("High","Low","None"))
  )

readr::write_csv(delta_gc, file.path(out_dir, "deltaGC_per_phage_vs_hostspecies.csv"))

# =================================================================
# Summarize to the genus level (takes into account the possibility
#                         of multiple host strains within a genus)
# =================================================================

# For making heatmaps

# based on tRNA group
delta_gc_trna_summary <- delta_gc %>%
  filter(!is.na(trna_group), !is.na(host_genus), !is.na(Delta_GC)) %>%
  group_by(host_genus, trna_group, GC_Type) %>%
  summarise(
    n_phages     = n_distinct(Phage),
    mean_deltaGC = mean(Delta_GC, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(delta_gc_trna_summary, file.path(out_dir, "deltaGC_mean_by_genus_trna_group.csv"))

# based on lifestyle
delta_gc_life_summary <- delta_gc %>%
  filter(!is.na(lifestyle), !is.na(host_genus), !is.na(Delta_GC)) %>%
  group_by(host_genus, lifestyle, GC_Type) %>%
  summarise(
    n_phages     = n_distinct(Phage),
    mean_deltaGC = mean(Delta_GC, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(delta_gc_life_summary, file.path(out_dir, "deltaGC_mean_by_genus_lifestyle.csv"))

# ==============
# Plot heatmaps
# ==============

heat_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(size = 14, face = "bold")
  )

# tRNA group heatmap
p_trna <- ggplot(delta_gc_trna_summary,
                 aes(x = trna_group, y = host_genus, fill = mean_deltaGC)) +
  geom_tile(color = "white") +
  facet_wrap(~ GC_Type, nrow = 1) +   # one row across
  scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#EF4444", midpoint = 0,
                       name = "Mean ΔGC") +
  labs(x = "tRNA group", y = "Host genus",
       title = "Mean ΔGC (Phage − Host strain) by tRNA group and genus") +
  heat_theme

ggsave(file.path(out_dir, "deltaGC_heatmap_trna_by_genus.png"),
       p_trna, width = 12, height = 5, dpi = 300)

# lifestyle heatmap
p_life <- ggplot(delta_gc_life_summary,
                 aes(x = lifestyle, y = host_genus, fill = mean_deltaGC)) +
  geom_tile(color = "white") +
  facet_wrap(~ GC_Type, nrow = 1) +   # one row across
  scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#EF4444", midpoint = 0,
                       name = "Mean ΔGC") +
  labs(x = "Lifestyle", y = "Host genus",
       title = "Mean ΔGC (Phage − Host strain) by lifestyle and genus") +
  heat_theme

ggsave(file.path(out_dir, "deltaGC_heatmap_lifestyle_by_genus.png"),
       p_life, width = 12, height = 5, dpi = 300)


