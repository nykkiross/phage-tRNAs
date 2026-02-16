# =========================
# CALCULATING tAI and ∆tAI
# =========================

# This script is designed with bacterial hosts and bacteriophages in mind
# it can be applied to other virus/host pairs :)

# Set base and output directories
# Load metadata - should have columns "sample_id" (organism name), "file_path" (file path for the specific FASTA sequence),
#                 "sample_type" (i.e. phage or host), "host_genus", "host_species", "lifestyle" (Virulent or Temperate),
#                 and "trna_group" (High, Low, or None)

base_dir      <- "/your/file/path/FASTAs"   # assumes FASTA files are all saved to a single parent folder
metadata_xlsx <- "/your/file/path/metadata.xlsx" 
metadata <- metadata %>%
  mutate(file_path = file.path(base_dir, file_path))

out_dir       <- "/your/file/path/Results/tAI"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load required packages

req_pkgs <- c(
  "dplyr","readr","tidyr","purrr","stringr","forcats",
  "ggplot2","glmmTMB","broom.mixed","rstatix","emmeans","tibble",
  "rlang")

not_installed <- req_pkgs[!req_pkgs %in% installed.packages()[,"Package"]]
if(length(not_installed)) install.packages(not_installed, dependencies = TRUE)
invisible(lapply(req_pkgs, library, character.only = TRUE))

theme_set(theme_minimal(base_size = 14))

# =======================================
# Calculate RSCU per gene (not averaged)
# =======================================

# necessary as the RSCU comparisons are performed at the genome level, but for tAI we need gene level

calculate_gene_rscu <- function(fasta_file, sample_id, sample_type, host_species, lifestyle = NA, trna_group = NA) {
  seqs <- read.fasta(fasta_file, seqtype = "DNA")
  rscu_list <- lapply(seqs, function(seq) {
    gene_name <- attr(seq, "name")
    codon_counts <- uco(seq, index = "rscu")
    tibble(
      sample_id = sample_id,
      gene_id = gene_name,
      codon = tolower(names(codon_counts)),
      RSCU = as.numeric(codon_counts),
      sample_type = sample_type,
      host_species = host_species,
      lifestyle = lifestyle,
      trna_group = trna_group
    )
  })
  bind_rows(rscu_list)
}

all_gene_rscu <- metadata %>%
  mutate(
    rscu = pmap(list(file_path, sample_id, sample_type, host_species, lifestyle, trna_group),
                ~ calculate_gene_rscu(..1, ..2, ..3, ..4, ..5, ..6))
  ) %>%
  pull(rscu) %>%
  bind_rows()

write_csv(all_gene_rscu, "rscu_tidy_per_gene.csv")

# ============================================
# Load and process tRNAscan-SE 2.0 .txt files
# ============================================

# Assumes file has columns: Name, tRNA_Type, Anti_Codon, Score
process_trnascan_txt <- function(file_path, species_name, sample_type = "host") {
  df <- read_tsv(
    file_path,
    skip = 2,
    col_names = c("Name", "tRNA_number", "Begin", "End", 
                  "tRNA_Type", "Anti_Codon", "Intron_Begin", "Intron_End", 
                  "Score", "Isotype_CM", "Isotype_Score", "Note"),
    show_col_types = FALSE
  )
  
  df <- df %>%
    select(Name, tRNA_Type, Anti_Codon, Score) %>%
    mutate(
      host_species = species_name,
      sample_type = sample_type,
      Anti_Codon = tolower(Anti_Codon),
      Score = as.numeric(Score)
    ) %>%
    filter(!is.na(Score))  # 🔑 removes the dashed divider row
  
  return(df)
}

host_files <- list("your/file/paths/tRNAscan.txt",
                   "your/file/paths/tRNAscan2.txt",
                   "...","your/file/paths/tRNAscanN.txt")

host_trna_list <- lapply(names(host_files), function(species) {
  process_trnascan_txt(host_files[[species]], species)
})
names(host_trna_list) <- names(host_files)

# =============================
# Generate anticodon libraries
# =============================

generate_anticodon_dictionary <- function(trna_data) {
  tapply(as.numeric(trna_data$Score), trna_data$Anti_Codon, sum, na.rm = TRUE)
}

host_dicts <- lapply(host_trna_list, generate_anticodon_dictionary)

# ================================================
# Reverse complement for codon/anticodon matching
# ================================================

reverse_complement <- function(codon) {
  comp <- c("a" = "t", "t" = "a", "g" = "c", "c" = "g")
  paste0(rev(comp[strsplit(codon, NULL)[[1]]]), collapse = "")
}

# =========================
# tAI Calculation per gene
# =========================

# Calculates codon weights and per-gene tAI
# Normalizes per species and produces genome-level summary (avg across genomes)
calculate_tai_with_rscu <- function(rscu_data, anticodon_dict, min_weight = 0.01) {
  rscu_data$tAI_weight <- mapply(function(codon, rscu_value) {
    anticodon <- reverse_complement(codon)
    matching_scores <- anticodon_dict[anticodon]
    if (!is.na(matching_scores)) {
      (rscu_value + 1e-6) * pmax(matching_scores, min_weight)
    } else {
      min_weight
    }
  }, rscu_data$codon, rscu_data$RSCU, SIMPLIFY = TRUE)
  gene_tai <- rscu_data %>%
    group_by(sample_id, gene_id, sample_type, host_species, lifestyle, trna_group) %>%
    summarize(
      raw_tAI = exp(mean(log(pmax(tAI_weight, min_weight)), na.rm = TRUE)),
      .groups = "drop"
    )
  max_tai <- max(gene_tai$raw_tAI, na.rm = TRUE)
  gene_tai <- gene_tai %>% mutate(tAI = raw_tAI / max_tai)
  genome_tai <- gene_tai %>%
    group_by(sample_id, sample_type, host_species, lifestyle, trna_group) %>%
    summarize(
      mean_raw_tAI = mean(raw_tAI, na.rm = TRUE),
      mean_tAI     = mean(tAI, na.rm = TRUE),
      .groups = "drop"
    )
  list(
    gene_level   = gene_tai,
    genome_level = genome_tai
  )
}

# Loop over host species to calculate tAI and ∆tAI
all_results <- list()

for (species in names(host_dicts)) {
  message("Processing host: ", species)
  subset_rscu <- all_gene_rscu %>% filter(host_species == species)
  results <- calculate_tai_with_rscu(subset_rscu, host_dicts[[species]])
  genome_res <- results$genome_level
  host_mean <- mean(genome_res$mean_tAI[genome_res$sample_type == "host"], na.rm = TRUE)
  genome_res <- genome_res %>%
    mutate(delta_tAI = ifelse(sample_type == "phage", mean_tAI - host_mean, NA))
  all_results[[species]] <- genome_res
}

combined_results <- bind_rows(all_results)

# Optional: per-gene resolution of tAI
for (species in names(host_dicts)) {
  message("Processing host: ", species)
  subset_rscu <- all_gene_rscu %>% filter(host_species == species)
  results <- calculate_tai_with_rscu(subset_rscu, host_dicts[[species]])
  gene_res <- results$gene_level
  host_mean <- mean(gene_res$tAI[gene_res$sample_type == "host"], na.rm = TRUE)
  gene_res <- gene_res %>%
    mutate(delta_tAI = ifelse(sample_type == "phage", tAI - host_mean, NA))
  all_results[[species]] <- gene_res
}

combined_results <- bind_rows(all_results)

write_csv(combined_results, "tai_results_per_gene_all_hosts.csv")

# =================================================
# Load libraries and data for statistical analyses
# and pre-process gene-level data
# =================================================

tai_results <- read_csv("/your/file/path/Results/tAI/tai_results_all_hosts.csv")
gene_tai_results <- read_csv("/your/file/path/Results/tAI/tai_results_per_gene_all_hosts.csv")

eps <- 1e-6
gene_tai_results <- gene_tai_results %>%
  mutate(
    lifestyle  = factor(lifestyle, levels = c("Temperate","Virulent")),
    trna_group = factor(trna_group, levels = c("None","Low","High")),
    # keep inside (0,1) for beta regression
    tAI_beta   = pmin(pmax(tAI, eps), 1 - eps),
    # rescale ΔtAI into (0,1) for beta regression
    delta_tAI_rescaled = ifelse(!is.na(delta_tAI),
                                (delta_tAI - min(delta_tAI, na.rm = TRUE) + eps) /
                                  (max(delta_tAI, na.rm = TRUE) - min(delta_tAI, na.rm = TRUE) + 2*eps),
                                NA)
  )

# ===============================================
# Within genus models - tRNA group and lifestyle
# ===============================================

per_host_tai <- list()
per_host_delta <- list()

for (sp in unique(gene_tai_results$host_species)) {
  df <- gene_tai_results %>% filter(sample_type == "phage", host_species == sp)
  if (length(unique(na.omit(df$lifestyle))) > 1 || length(unique(na.omit(df$trna_group))) > 1) {
    mod_tai <- tryCatch(
      glmmTMB(tAI_beta ~ lifestyle + trna_group + (1 | sample_id),
              family = beta_family(), data = df),
      error = function(e) NULL
    )
    if (!is.null(mod_tai)) {
      per_host_tai[[sp]] <- broom.mixed::tidy(mod_tai, effects = "fixed") %>%
        mutate(host_species = sp, response = "tAI")
    }
  }
  df_delta <- df %>% filter(!is.na(delta_tAI_rescaled))
  if (nrow(df_delta) > 0 &&
      (length(unique(na.omit(df_delta$lifestyle))) > 1 ||
       length(unique(na.omit(df_delta$trna_group))) > 1)) {
    mod_delta <- tryCatch(
      glmmTMB(delta_tAI_rescaled ~ lifestyle + trna_group + (1 | sample_id),
              family = beta_family(), data = df_delta),
      error = function(e) NULL
    )
    if (!is.null(mod_delta)) {
      per_host_delta[[sp]] <- broom.mixed::tidy(mod_delta, effects = "fixed") %>%
        mutate(host_species = sp, response = "delta_tAI")
    }
  }
}

per_host_tai <- bind_rows(per_host_tai)
per_host_delta <- bind_rows(per_host_delta)

write_csv(per_host_tai, file.path(out_dir, "mixedmodel_perhost_tAI.csv"))
write_csv(per_host_delta, file.path(out_dir, "mixedmodel_perhost_delta_tAI.csv"))

# ========================================
# Global model - tRNA group and lifestyle
# ========================================

gene_tai_results <- gene_tai_results %>%
  mutate(host_species = relevel(factor(host_species), ref = "E. coli"))

m_global_tai <- glmmTMB(
  tAI_beta ~ lifestyle + trna_group + host_species + (1 | sample_id),
  family = beta_family(),
  data = gene_tai_results %>% filter(sample_type == "phage")
)

tidy_global_tai <- broom.mixed::tidy(m_global_tai, effects = "fixed")
write_csv(tidy_global_tai, file.path(out_dir, "mixedmodel_global_tAI.csv"))

m_global_delta <- glmmTMB(
  delta_tAI_rescaled ~ lifestyle + trna_group + (1|sample_id),
  family = beta_family(),
  data = gene_tai_results %>% filter(sample_type == "phage", !is.na(delta_tAI))
)

tidy_global_delta <- broom.mixed::tidy(m_global_delta, effects = "fixed")
write_csv(tidy_global_delta, file.path(out_dir, "mixedmodel_global_delta_tAI.csv"))

# =====================
# Plotting preparation
# =====================

# failsafe - bring in genus from metadata
meta_genus <- readr::read_csv(
  "/Users/nross1994/Documents/Doore Lab/codon usage phages/metadata/all phage metadata.csv",
  show_col_types = FALSE
) |>
  dplyr::select(sample_id, host_genus) |>
  dplyr::distinct()

# Helper - expands abbreviated genus to full genus and prevents issues further down the pipeline
expand_genus <- function(host_species) {
  abbr_map <- c(
    "A." = "HostA",
    "B." = NA_character_, # if more than 1 host has the same abbreviation, use this
    "..." = "...",
    "N." = "HostN"
  )
  # species-based disambiguation
  species_map <- list(
    "1"      = "B1",
    "2"      = "B2",
    "..."    = "...",
    "N"      = "BN"
  )
  out <- rep(NA_character_, length(host_species))
  sp <- trimws(host_species)
  first_token <- sub("\\s+.*$", "", sp)
  spelled_out <- grepl("^[A-Z][a-z]+$", first_token)
  out[spelled_out] <- first_token[spelled_out]
  abbr_token <- sub("\\s+.*$", "", sp)
  is_abbr <- grepl("^[A-Z]\\.$", abbr_token)
  if (any(is_abbr)) {
    initial <- abbr_token[is_abbr]
    genus_guess <- abbr_map[initial]
    second_token <- sub("^[^\\s]+\\s+", "", sp[is_abbr])
    second_token <- sub("^([A-Za-z]+).*", "\\1", second_token)  # strip subspecies
    needs_disambig <- initial == "S."
    if (any(needs_disambig)) {
      sec <- tolower(second_token[needs_disambig])
      genus_from_species <- vapply(sec, function(w) {
        # try direct
        g <- species_map[[w]]
        if (is.null(g)) NA_character_ else g
      }, character(1))
      genus_guess[needs_disambig] <- genus_from_species
    }
    out[is_abbr & is.na(out)] <- genus_guess
  }
  
  out
}

# Add host_genus to tai_results
tai_results <- tai_results |>
  dplyr::left_join(meta_genus, by = "sample_id") |>
  dplyr::mutate(
    host_genus = dplyr::coalesce(host_genus, expand_genus(host_species))
  )

# Add host_genus to gene_tai_results
gene_tai_results <- gene_tai_results |>
  dplyr::left_join(meta_genus, by = "sample_id") |>
  dplyr::mutate(
    host_genus = dplyr::coalesce(host_genus, expand_genus(host_species))
  )

# ===============================
# PLOTTING tAI & ΔtAI COMPARISONS
# ===============================

# Make sure factor order is consistent
tai_results <- tai_results %>%
  mutate(
    lifestyle   = factor(lifestyle, levels = c("Temperate","Virulent")),
    trna_group  = factor(trna_group, levels = c("None","Low","High")),
    sample_type = factor(sample_type, levels = c("host","phage")),
    host_genus  = factor(
      host_genus,
      levels = c("Host1", "Host2", "...", "HostN")))

# Settings - can be changed
palette_lifestyle <- c(Virulent = "#1f78b4", Temperate = "#33a02c")
palette_trna      <- c(None = "#5e3c99", Low = "#b2abd2", High = "#fdb863")

theme_large <- theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 16, face = "bold"),
    strip.text  = element_text(size = 16, face = "bold"),
    plot.title  = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )

clean_labels <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(x)
}

clean_theme <- theme(
  panel.grid = element_blank(),   # remove gridlines
  legend.title = element_text(size = 13),
  legend.text  = element_text(size = 12)
)

# Host vs. phages tAI Violin plot
p_host_vs_phage <- ggplot(tai_results, aes(x = host_genus, y = mean_tAI)) +
  geom_violin(data = filter(tai_results, sample_type == "phage"),
              aes(fill = "Phages"), trim = FALSE) +
  geom_boxplot(data = filter(tai_results, sample_type == "phage"),
               width = 0.15, alpha = 0.5, outlier.shape = NA) +
  geom_point(data = filter(tai_results, sample_type == "host"),
             aes(color = "Host"), size = 3, shape = 18) +
  scale_fill_manual(values = c("Phages" = "steelblue")) +
  scale_color_manual(values = c("Host" = "red")) +
  labs(title = "Host vs Phage tAI by Host Genus", y = "tAI", x = "Host Genus",
       fill = NULL, color = NULL) +
  theme_large +
  clean_theme

ggsave(file.path(out_dir, "violin_tAI_host_vs_phage.png"),
       p_host_vs_phage, width = 11, height = 7, dpi = 300)

# ∆tAI by tRNA group - faceted by genus
p_trna_delta <- tai_results %>%
  filter(sample_type == "phage") %>%
  ggplot(aes(x = trna_group, y = delta_tAI, fill = trna_group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55) +
  facet_wrap(~ host_genus, scales = "free_y",
             labeller = labeller(host_genus = clean_labels)) +
  scale_fill_manual(values = palette_trna, labels = clean_labels,
                    name = clean_labels("lifestyle")) +
  scale_x_discrete(labels = clean_labels) +
  coord_cartesian(ylim = c(-0.5, 0.4)) +
  labs(title = "ΔtAI vs tRNA Group", x = "tRNA Group", y = "ΔtAI") +
  theme_large +
  clean_theme

ggsave(file.path(out_dir, "violin_delta_tAI_by_trna_group.png"),
       p_trna_delta, width = 12, height = 8, dpi = 300)

# ∆tAI by lifestyle - faceted by genus
p_lifestyle_delta <- tai_results %>%
  filter(sample_type == "phage") %>%
  ggplot(aes(x = lifestyle, y = delta_tAI, fill = lifestyle)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55) +
  facet_wrap(~ host_genus, scales = "free_y",
             labeller = labeller(host_genus = clean_labels)) +
  scale_fill_manual(values = palette_lifestyle, labels = clean_labels,
                    name = clean_labels("trna_group")) +
  scale_x_discrete(labels = clean_labels) +
  coord_cartesian(ylim = c(-0.5, 0.4)) +
  labs(title = "ΔtAI vs Lifestyle", x = "Lifestyle", y = "ΔtAI") +
  theme_large +
  clean_theme

ggsave(file.path(out_dir, "violin_delta_tAI_by_lifestyle.png"),
       p_lifestyle_delta, width = 12, height = 8, dpi = 300)

# Global ∆tAI by tRNA group
p_trna_delta_global <- tai_results %>%
  filter(sample_type == "phage") %>%
  ggplot(aes(x = trna_group, y = delta_tAI, fill = trna_group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55) +
  scale_fill_manual(values = palette_trna, labels = clean_labels,
                    name = clean_labels("trna_group")) +
  scale_x_discrete(labels = clean_labels) +
  labs(title = "Global ΔtAI by tRNA Group", x = "tRNA Group", y = "ΔtAI") +
  theme_large +
  clean_theme

ggsave(file.path(out_dir, "violin_delta_tAI_global_by_trna_group.png"),
       p_trna_delta_global, width = 8, height = 6, dpi = 300)

# Global ∆tAI by lifestyle
p_life_delta_global <- tai_results %>%
  filter(sample_type == "phage") %>%
  ggplot(aes(x = lifestyle, y = delta_tAI, fill = lifestyle)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55) +
  scale_fill_manual(values = palette_lifestyle, labels = clean_labels,
                    name = clean_labels("lifestyle")) +
  scale_x_discrete(labels = clean_labels) +
  labs(title = "Global ΔtAI by Lifestyle", x = "Lifestyle", y = "ΔtAI") +
  theme_large +
  clean_theme

ggsave(file.path(out_dir, "violin_delta_tAI_global_by_lifestyle.png"),
       p_life_delta_global, width = 8, height = 6, dpi = 300)

# OPTIONAL: tAI by tRNA group
p_trna_tai <- tai_results %>%
  filter(sample_type == "phage") %>%
  ggplot(aes(x = trna_group, y = mean_tAI, fill = trna_group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55) +
  facet_wrap(~ host_genus, scales = "free_y", 
             labeller = labeller(host_genus = clean_labels)) +
  scale_fill_manual(values = palette_trna, labels = clean_labels,
                    name = clean_labels("trna_group")) +
  scale_x_discrete(labels = clean_labels) +
  coord_cartesian(ylim = c(0, 0.7)) +
  labs(title = "tRNA Group vs Phage tAI", x = "tRNA Group", y = "tAI") +
  theme_large +
  clean_theme

ggsave(file.path(out_dir, "violin_tAI_by_trna_group.png"),
       p_trna_tai, width = 12, height = 8, dpi = 300)

# OPTIONAL: tAI by lifestyle
p_lifestyle_tai <- tai_results %>%
  filter(sample_type == "phage") %>%
  ggplot(aes(x = lifestyle, y = mean_tAI, fill = lifestyle)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55) +
  facet_wrap(~ host_genus, scales = "free_y",
             labeller = labeller(host_genus = clean_labels)) +
  scale_fill_manual(values = palette_lifestyle, labels = clean_labels,
                    name = clean_labels("lifestyle")) +
  scale_x_discrete(labels = clean_labels) +
  coord_cartesian(ylim = c(0, 0.7)) +
  labs(title = "Phage Lifestyle vs tAI", x = "Lifestyle", y = "tAI") +
  theme_large +
  clean_theme

ggsave(file.path(out_dir, "violin_tAI_by_lifestyle.png"),
       p_lifestyle_tai, width = 12, height = 8, dpi = 300)

# =========================================
# Statistical analyses - ∆tAI ~ 0 by genus
# =========================================

# If host_genus is missing in tai_results, bring it in from metadata
if (!"host_genus" %in% names(tai_results)) {
  metadata <- read_csv("/Users/nross1994/Documents/Doore Lab/codon usage phages/metadata/all phage metadata.csv") %>%
    select(sample_id, host_genus)
  tai_results <- tai_results %>%
    left_join(metadata, by = "sample_id")
}

phage_only <- tai_results %>%
  filter(sample_type == "phage", !is.na(delta_tAI), !is.na(host_genus))

genus_counts <- phage_only %>% dplyr::count(host_genus, name = "n_phages")
phage_only_kept <- phage_only %>%
  inner_join(genus_counts %>% filter(n_phages >= 3), by = "host_genus")

# Wilcoxon one-sample (median ΔtAI vs 0) per genus
stats_delta_by_genus <- phage_only_kept %>%
  group_by(host_genus) %>%
  wilcox_test(delta_tAI ~ 1, mu = 0, alternative = "two.sided") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  # add n / median / IQR for effect-size context
  left_join(
    phage_only_kept %>%
      group_by(host_genus) %>%
      summarise(
        n = n(),
        median_delta = median(delta_tAI, na.rm = TRUE),
        iqr_delta = IQR(delta_tAI, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "host_genus"
  ) %>%
  arrange(p.adj)

readr::write_csv(stats_delta_by_genus,
                 file.path(out_dir, "stats_delta_tAI_phage_vs_zero_by_genus.csv"))

# Summary stats
delta_summary_by_genus <- phage_only %>%
  group_by(host_genus) %>%
  summarise(
    n = n(),
    mean_delta   = mean(delta_tAI, na.rm = TRUE),
    sd_delta     = sd(delta_tAI, na.rm = TRUE),
    median_delta = median(delta_tAI, na.rm = TRUE),
    iqr_delta    = IQR(delta_tAI, na.rm = TRUE),
    min_delta    = min(delta_tAI, na.rm = TRUE),
    max_delta    = max(delta_tAI, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(delta_summary_by_genus,
                 file.path(out_dir, "summary_delta_tAI_by_genus.csv"))

# =======================
# Per-genus mixed models
# =======================

# Ensure host_genus present; attach if missing
if (!"host_genus" %in% names(gene_tai_results)) {
  metadata_path <- "/Users/nross1994/Documents/Doore Lab/codon usage phages/metadata/all phage metadata.csv"
  metadata <- read_csv(metadata_path, show_col_types = FALSE) %>%
    select(sample_id, host_genus)
  gene_tai_results <- gene_tai_results %>%
    left_join(metadata, by = "sample_id")
}

# Standardize factor levels
gene_tai_results <- gene_tai_results %>%
  mutate(
    lifestyle  = tolower(lifestyle),
    trna_group = tolower(trna_group),
    sample_type = tolower(sample_type),
    host_genus  = as.character(host_genus),
    tAI_beta = pmin(pmax(tAI, 1e-6), 1 - 1e-6)
  )

# Helper to build formula given available predictors in a genus
.build_formula <- function(resp, has_life, has_trna) {
  rhs <- c()
  if (has_life) rhs <- c(rhs, "lifestyle")
  if (has_trna) rhs <- c(rhs, "trna_group")
  rhs <- c(rhs, "(1|sample_id)")
  as.formula(paste(resp, "~", paste(rhs, collapse = " + ")))
}

per_genus_tai   <- list()
per_genus_delta <- list()

# Loop per host_genus
for (g in sort(unique(na.omit(gene_tai_results$host_genus)))) {
  df_g <- gene_tai_results %>%
    filter(sample_type == "phage", host_genus == g) %>%
    mutate(
      lifestyle  = droplevels(factor(lifestyle,  levels = c("temperate","virulent"))),
      trna_group = droplevels(factor(trna_group, levels = c("none","low","high")))
    )
  
  if (nrow(df_g) < 50 || length(unique(df_g$sample_id)) < 3) next
  
  has_life <- nlevels(df_g$lifestyle)  > 1
  has_trna <- nlevels(df_g$trna_group) > 1
  
  if (has_life || has_trna) {
    form_tai <- .build_formula("tAI_beta", has_life, has_trna)
    mod_tai <- tryCatch(
      glmmTMB(form_tai, family = beta_family(), data = df_g),
      error = function(e) NULL
    )
    if (!is.null(mod_tai)) {
      per_genus_tai[[g]] <- tidy(mod_tai, effects = "fixed") %>%
        mutate(host_genus = g, response = "tAI")
    }
  }
  
  df_delta <- df_g %>%
    filter(!is.na(delta_tAI)) %>%
    mutate(
      delta_min = min(delta_tAI, na.rm = TRUE),
      delta_max = max(delta_tAI, na.rm = TRUE),
      delta_tAI_rescaled = (delta_tAI - delta_min + 1e-6) /
        (delta_max - delta_min + 2e-6)
    ) %>%
    select(-delta_min, -delta_max)
  
  if (nrow(df_delta) >= 50) {
    has_life_d <- nlevels(droplevels(df_delta$lifestyle))  > 1
    has_trna_d <- nlevels(droplevels(df_delta$trna_group)) > 1
    
    if (has_life_d || has_trna_d) {
      form_delta <- .build_formula("delta_tAI_rescaled", has_life_d, has_trna_d)
      mod_delta <- tryCatch(
        glmmTMB(form_delta, family = beta_family(), data = df_delta),
        error = function(e) NULL
      )
      if (!is.null(mod_delta)) {
        per_genus_delta[[g]] <- tidy(mod_delta, effects = "fixed") %>%
          mutate(host_genus = g, response = "delta_tAI")
      }
    }
  }
}

per_genus_tai_df   <- bind_rows(per_genus_tai)
per_genus_delta_df <- bind_rows(per_genus_delta)

if (nrow(per_genus_tai_df)) {
  write_csv(per_genus_tai_df, file.path(out_dir, "mixedmodel_per_genus_tAI.csv"))
}
if (nrow(per_genus_delta_df)) {
  write_csv(per_genus_delta_df, file.path(out_dir, "mixedmodel_per_genus_delta_tAI.csv"))
}

write_csv(per_genus_tai_df,   file.path(out_dir, "mixedmodel_per_genus_tAI.csv"))
write_csv(per_genus_delta_df, file.path(out_dir, "mixedmodel_per_genus_delta_tAI.csv"))

# Summary stats - tAI per genus
# phage only
tai_summary_genus <- tai_results %>%
  filter(sample_type == "phage") %>%
  group_by(host_genus) %>%
  summarise(
    n        = n(),
    mean_tAI = mean(mean_tAI, na.rm = TRUE),
    sd_tAI   = sd(mean_tAI, na.rm = TRUE),
    .groups  = "drop"
  )
write_csv(tai_summary_genus,
          file.path(out_dir, "tai_summary_stats_by_genus.csv"))

# hosts and phages
tai_summary_genus <- tai_results %>%
  group_by(host_genus, sample_type) %>%
  summarise(
    n         = n(),
    median_tAI = median(mean_tAI, na.rm = TRUE),
    sd_tAI     = sd(mean_tAI, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = host_genus,
    names_from  = sample_type,
    values_from = c(n, median_tAI, sd_tAI),
    names_sep   = "_"
  )
write_csv(tai_summary_genus,
          file.path(out_dir, "tai_summary_stats_by_genus_host_phage.csv"))

# ==========================
# Global statistical models
# ==========================

# Make sure factors + beta-safe columns exist
eps <- 1e-6
gene_tai_results <- gene_tai_results %>%
  mutate(
    lifestyle  = factor(tolower(lifestyle),  levels = c("temperate","virulent")),
    trna_group = factor(tolower(trna_group), levels = c("none","low","high")),
    sample_type = factor(tolower(sample_type), levels = c("host","phage")),
    host_genus  = as.factor(host_genus),
    # squeeze tAI into (0,1)
    tAI_beta = pmin(pmax(tAI, eps), 1 - eps)
  )

df_phage <- gene_tai_results %>% filter(sample_type == "phage")

has_life <- nlevels(droplevels(df_phage$lifestyle))  > 1
has_trna <- nlevels(droplevels(df_phage$trna_group)) > 1

form_tai <- as.formula(
  paste0("tAI_beta ~ ",
         paste(c(if (has_life) "lifestyle", if (has_trna) "trna_group", "host_genus"), collapse = " + "),
         " + (1|sample_id)")
)

m_global_tai <- glmmTMB(form_tai, family = beta_family(), data = df_phage)

# Predicted means by lifestyle
emm_global_tai_life <- emmeans(m_global_tai, ~ lifestyle, type = "response") %>%
  as_tibble() %>%
  mutate(response = "tAI")

# Predicted means by tRNA group
emm_global_tai_trna <- emmeans(m_global_tai, ~ trna_group, type = "response") %>%
  as_tibble() %>%
  mutate(response = "tAI")

# Optional: pairwise contrasts
contr_global_tai_life <- contrast(emmeans(m_global_tai, ~ lifestyle, type = "response"), method = "pairwise") %>% as_tibble()
contr_global_tai_trna <- contrast(emmeans(m_global_tai, ~ trna_group, type = "response"), method = "pairwise") %>% as_tibble()

write_csv(emm_global_tai_life, file.path(out_dir, "emm_global_tAI_by_lifestyle.csv"))
write_csv(emm_global_tai_trna, file.path(out_dir, "emm_global_tAI_by_trna.csv"))
write_csv(contr_global_tai_life, file.path(out_dir, "emm_global_tAI_contrasts_lifestyle.csv"))
write_csv(contr_global_tai_trna, file.path(out_dir, "emm_global_tAI_contrasts_trna.csv"))

# Rescale ΔtAI to (0,1) globally for beta family
eps <- 1e-6
df_phage <- gene_tai_results %>% dplyr::filter(sample_type == "phage")

delta_min_g <- min(df_phage$delta_tAI, na.rm = TRUE)
delta_max_g <- max(df_phage$delta_tAI, na.rm = TRUE)

df_delta <- df_phage %>%
  dplyr::filter(!is.na(delta_tAI)) %>%
  dplyr::mutate(delta_tAI_rescaled = (delta_tAI - delta_min_g + eps) /
                  (delta_max_g - delta_min_g + 2*eps))

has_life_d <- nlevels(droplevels(df_delta$lifestyle))  > 1
has_trna_d <- nlevels(droplevels(df_delta$trna_group)) > 1

form_delta <- as.formula(
  paste0("delta_tAI_rescaled ~ ",
         paste(c(if (has_life_d) "lifestyle", if (has_trna_d) "trna_group", "host_genus"),
               collapse = " + "),
         " + (1|sample_id)")
)

m_global_delta <- glmmTMB(form_delta, family = beta_family(), data = df_delta)

# helper to back-transform from rescaled (0..1) to original ΔtAI
unscale_delta <- function(mu_rescaled) {
  mu_rescaled * (delta_max_g - delta_min_g + 2*eps) - eps + delta_min_g
}

# lifestyle emmeans
emm_global_delta_life <- emmeans(m_global_delta, ~ lifestyle, type = "response") %>%
  as_tibble() %>%
  mutate(
    resp = response,  # copy numeric column
    mean_delta  = unscale_delta(resp),
    lower_delta = unscale_delta(asymp.LCL),
    upper_delta = unscale_delta(asymp.UCL),
    measure = "delta_tAI"
  ) %>%
  select(lifestyle, resp, asymp.LCL, asymp.UCL, mean_delta, lower_delta, upper_delta, measure)

# tRNA group emmeans
emm_global_delta_trna <- emmeans(m_global_delta, ~ trna_group, type = "response") %>%
  as_tibble() %>%
  mutate(
    resp = response,
    mean_delta  = unscale_delta(resp),
    lower_delta = unscale_delta(asymp.LCL),
    upper_delta = unscale_delta(asymp.UCL),
    measure = "delta_tAI"
  ) %>%
  select(trna_group, resp, asymp.LCL, asymp.UCL, mean_delta, lower_delta, upper_delta, measure)

write_csv(emm_global_delta_life, file.path(out_dir, "emm_global_delta_tAI_by_lifestyle.csv"))
write_csv(emm_global_delta_trna, file.path(out_dir, "emm_global_delta_tAI_by_trna.csv"))

# ================================================================
# OPTIONAL ADDITIONAL TESTS (were not included in our manuscript)
# ================================================================

per_genus_emm_tai   <- list()
per_genus_emm_delta <- list()


# Standardize emmeans column names to: pred, lcl, ucl
get_emm_df <- function(emm_obj) {
  df <- as.data.frame(emm_obj, stringsAsFactors = FALSE)

  pred_col <- intersect(c("response", "emmean", "prob"), names(df))
  lcl_col  <- intersect(c("asymp.LCL", "lower.CL", "LCL"), names(df))
  ucl_col  <- intersect(c("asymp.UCL", "upper.CL", "UCL"), names(df))
  
  if (length(pred_col) == 0) stop("No predicted mean column in emmeans output.")
  if (length(lcl_col)  == 0) stop("No lower CI column in emmeans output.")
  if (length(ucl_col)  == 0) stop("No upper CI column in emmeans output.")
  
  names(df)[names(df) == pred_col[1]] <- "pred"
  names(df)[names(df) == lcl_col[1]]  <- "lcl"
  names(df)[names(df) == ucl_col[1]]  <- "ucl"
  tibble::as_tibble(df)
}

eps <- 1e-6

# Ensure lowercase factors
gene_tai_results <- gene_tai_results %>%
  mutate(
    lifestyle  = factor(tolower(lifestyle),  levels = c("temperate","virulent")),
    trna_group = factor(tolower(trna_group), levels = c("none","low","high")),
    sample_type = factor(tolower(sample_type), levels = c("host","phage")),
    tAI_beta = pmin(pmax(tAI, eps), 1 - eps)
  )

# Helper to build formulas 
.build_formula <- function(resp, has_life, has_trna) {
  rhs <- c()
  if (has_life) rhs <- c(rhs, "lifestyle")
  if (has_trna) rhs <- c(rhs, "trna_group")
  rhs <- c(rhs, "(1|sample_id)")
  as.formula(paste(resp, "~", paste(rhs, collapse = " + ")))
}

per_genus_emm_tai   <- list()
per_genus_emm_delta <- list()
per_genus_log       <- list()

# You can relax/tighten these thresholds:
MIN_ROWS    <- 20   # require at least this many phage-gene rows
MIN_GENOMES <- 2    # require at least this many distinct phage genomes

for (g in sort(unique(na.omit(gene_tai_results$host_genus)))) {
  log_msgs <- c()
  df_g <- gene_tai_results %>%
    filter(sample_type == "phage", host_genus == g) %>%
    mutate(
      lifestyle  = droplevels(lifestyle),
      trna_group = droplevels(trna_group)
    )
  
  n_rows <- nrow(df_g)
  n_genomes <- n_distinct(df_g$sample_id)
  if (n_rows < MIN_ROWS) {
    log_msgs <- c(log_msgs, sprintf("SKIP %s: n_rows=%d < %d", g, n_rows, MIN_ROWS))
  }
  if (n_genomes < MIN_GENOMES) {
    log_msgs <- c(log_msgs, sprintf("SKIP %s: n_genomes=%d < %d", g, n_genomes, MIN_GENOMES))
  }
  
  if (length(log_msgs) == 0) {
    has_life <- nlevels(df_g$lifestyle)  > 1
    has_trna <- nlevels(df_g$trna_group) > 1
    
    if (has_life || has_trna) {
      form_tai_g <- .build_formula("tAI_beta", has_life, has_trna)
      mod_tai_g <- tryCatch(glmmTMB(form_tai_g, family = beta_family(), data = df_g), error = function(e) e)
      
      if (inherits(mod_tai_g, "error")) {
        log_msgs <- c(log_msgs, sprintf("tAI model failed for %s: %s", g, mod_tai_g$message))
      } else {
        if (has_life) {
          tmp <- tryCatch(emmeans(mod_tai_g, ~ lifestyle, type = "response"), error = function(e) e)
          if (!inherits(tmp, "error")) {
            df_tmp <- get_emm_df(tmp)
            per_genus_emm_tai[[paste0(g, "_life")]] <- df_tmp %>%
              dplyr::mutate(host_genus = as.character(g),
                            measure = "tAI",
                            predictor = "lifestyle") %>%
              dplyr::rename(mean = pred) %>%
              dplyr::select(host_genus, measure, predictor, lifestyle, mean, lcl, ucl)
          } else {
            log_msgs <- c(log_msgs, sprintf("tAI emmeans(lifestyle) failed for %s: %s", g, tmp$message))
          }
        }
        if (has_trna) {
          tmp <- tryCatch(emmeans(mod_tai_g, ~ trna_group, type = "response"), error = function(e) e)
          if (!inherits(tmp, "error")) {
            df_tmp <- get_emm_df(tmp)
            per_genus_emm_tai[[paste0(g, "_trna")]] <- df_tmp %>%
              dplyr::mutate(host_genus = as.character(g),
                            measure = "tAI",
                            predictor = "trna_group") %>%
              dplyr::rename(mean = pred) %>%
              dplyr::select(host_genus, measure, predictor, trna_group, mean, lcl, ucl)
          } else {
            log_msgs <- c(log_msgs, sprintf("tAI emmeans(trna_group) failed for %s: %s", g, tmp$message))
          }
        }
      }
    } else {
      log_msgs <- c(log_msgs, sprintf("tAI: %s has only one level for both predictors (lifestyle/trna_group).", g))
    }
    
    df_delta_g <- df_g %>% filter(!is.na(delta_tAI))
    if (nrow(df_delta_g) == 0) {
      log_msgs <- c(log_msgs, sprintf("delta_tAI: %s has no non-NA ΔtAI.", g))
    } else {
      dmin <- suppressWarnings(min(df_delta_g$delta_tAI, na.rm = TRUE))
      dmax <- suppressWarnings(max(df_delta_g$delta_tAI, na.rm = TRUE))
      if (!is.finite(dmin) || !is.finite(dmax) || dmin == dmax) {
        log_msgs <- c(log_msgs, sprintf("delta_tAI: %s has no variation (min==max) or non-finite.", g))
      } else {
        df_delta_g <- df_delta_g %>%
          mutate(delta_tAI_rescaled = (delta_tAI - dmin + eps) / (dmax - dmin + 2*eps))
        
        has_life_d <- nlevels(df_delta_g$lifestyle)  > 1
        has_trna_d <- nlevels(df_delta_g$trna_group) > 1
        
        if (has_life_d || has_trna_d) {
          form_delta_g <- .build_formula("delta_tAI_rescaled", has_life_d, has_trna_d)
          mod_delta_g <- tryCatch(glmmTMB(form_delta_g, family = beta_family(), data = df_delta_g), error = function(e) e)
          if (inherits(mod_delta_g, "error")) {
            log_msgs <- c(log_msgs, sprintf("delta_tAI model failed for %s: %s", g, mod_delta_g$message))
          } else {
            unscale_g <- function(mu) mu * (dmax - dmin + 2*eps) - eps + dmin
            if (has_life_d) {
              tmp <- tryCatch(emmeans(mod_delta_g, ~ lifestyle, type = "response"), error = function(e) e)
              if (!inherits(tmp, "error")) {
                df_tmp <- get_emm_df(tmp)
                per_genus_emm_delta[[paste0(g, "_life")]] <- df_tmp %>%
                  dplyr::mutate(
                    host_genus  = as.character(g),
                    measure     = "delta_tAI",
                    predictor   = "lifestyle",
                    mean_delta  = unscale_g(pred),
                    lower_delta = unscale_g(lcl),
                    upper_delta = unscale_g(ucl)
                  ) %>%
                  dplyr::select(host_genus, measure, predictor, lifestyle,
                                mean_delta, lower_delta, upper_delta)
              } else {
                log_msgs <- c(log_msgs, sprintf("delta emmeans(lifestyle) failed for %s: %s", g, tmp$message))
              }
            }
            if (has_trna_d) {
              tmp <- tryCatch(emmeans(mod_delta_g, ~ trna_group, type = "response"), error = function(e) e)
              if (!inherits(tmp, "error")) {
                df_tmp <- get_emm_df(tmp)
                per_genus_emm_delta[[paste0(g, "_trna")]] <- df_tmp %>%
                  dplyr::mutate(
                    host_genus  = as.character(g),
                    measure     = "delta_tAI",
                    predictor   = "trna_group",
                    mean_delta  = unscale_g(pred),
                    lower_delta = unscale_g(lcl),
                    upper_delta = unscale_g(ucl)
                  ) %>%
                  dplyr::select(host_genus, measure, predictor, trna_group,
                                mean_delta, lower_delta, upper_delta)
              } else {
                log_msgs <- c(log_msgs, sprintf("delta emmeans(trna_group) failed for %s: %s", g, tmp$message))
              }
            }
          }
        } else {
          log_msgs <- c(log_msgs, sprintf("delta_tAI: %s has single level for both predictors.", g))
        }
      }
    }
  }
  per_genus_log[[g]] <- tibble(host_genus = g, note = paste(log_msgs, collapse = " | "))
}

emm_per_genus_tai   <- bind_rows(per_genus_emm_tai)
emm_per_genus_delta <- bind_rows(per_genus_emm_delta)
log_df              <- bind_rows(per_genus_log)

if (nrow(emm_per_genus_tai))   write_csv(emm_per_genus_tai,   file.path(out_dir, "emm_per_genus_tAI.csv"))
if (nrow(emm_per_genus_delta)) write_csv(emm_per_genus_delta, file.path(out_dir, "emm_per_genus_delta_tAI.csv"))
write_csv(log_df, file.path(out_dir, "emm_per_genus_log.csv"))


