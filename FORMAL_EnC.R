# ===========================================
# CALCULATING EnC, ∆EnC, and EnC/GC3 CONTENT
# ===========================================

# This script is designed with bacterial hosts and bacteriophages in mind
# it can be applied to other virus/host pairs :)

# Set base and output directories
# Load metadata - should have columns "sample_id" (organism name), "file_path" (file path for the specific FASTA sequence),
#                 "sample_type" (i.e. phage or host), "host_genus", "host_species", "lifestyle" (Virulent or Temperate),
#                 and "trna_group" (High, Low, or None)

base_dir      <- "/your/file/path/FASTAs"   # assumes FASTA files are all saved to a single parent folder
metadata_xlsx <- "/your/file/path/metadata.xlsx" 
out_dir       <- "/your/file/path/Results/EnC"

# Load required packages

req_pkgs <- c("readxl", "dplyr", "tidyr", "purrr", "stringr", "seqinr", "ggplot2", "FSA")

not_installed <- req_pkgs[!req_pkgs %in% installed.packages()[,"Package"]]
if(length(not_installed)) install.packages(not_installed, dependencies = TRUE)
invisible(lapply(req_pkgs, library, character.only = TRUE))

# there may be issues with dplyr being masked, so load that again to avoid
library(dplyr)

# ===================================================
# Helper functions for calculating and comparing EnC
# ===================================================

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

standardize_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_trim(x)
  x
}

safe_read_fasta <- function(path) {
  read.fasta(path, seqtype = "DNA")
}

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_sum <- file.path(out_dir, "summaries")
if (!dir.exists(out_sum)) dir.create(out_sum, recursive = TRUE)

kw_eta2 <- function(kw_obj, n) {
  stat <- tryCatch(unname(as.numeric(kw_obj$statistic)), error = function(e) NA_real_)
  if (is.null(kw_obj) || is.na(stat) || is.null(n) || is.na(n) || n <= 1) return(NA_real_)
  stat / (n - 1)
}

cliffs_delta <- function(x, g) {
  g <- droplevels(factor(g))
  if (nlevels(g) != 2) return(NA_real_)
  x1 <- x[g == levels(g)[1]]
  x2 <- x[g == levels(g)[2]]
  if (length(x1) == 0 || length(x2) == 0) return(NA_real_)
  m <- outer(x1, x2, "-")
  (sum(m > 0) - sum(m < 0)) / (length(x1) * length(x2))
}

safe_kw <- function(formula, data) {
  out <- try(suppressWarnings(kruskal.test(formula, data = data)), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  stat <- tryCatch(unname(as.numeric(out$statistic)), error = function(e) NA_real_)
  if (length(stat) == 0 || is.na(stat)) return(NULL)
  out
}

# ================
# EnC calculation
# ================

# takes into account potential small-n groups

compute_sample_Ks <- function(genome, pseudocount = 1.0) {
  first_counts <- uco(genome[[1]], index = "eff", frame = 0)
  total_counts <- rep(0, length(first_counts))
  names(total_counts) <- toupper(names(first_counts))
  for (seq in genome) {
    cc <- uco(seq, index = "eff", frame = 0)
    cc <- cc[match(names(total_counts), toupper(names(cc)))]
    cc[is.na(cc)] <- 0
    total_counts <- total_counts + cc
  }
  codon_table <- list(
    K = c("AAA", "AAG"), N = c("AAC", "AAT"), Q = c("CAA", "CAG"),
    H = c("CAC", "CAT"), D = c("GAC", "GAT"), E = c("GAA", "GAG"),
    Y = c("TAC", "TAT"), C = c("TGC", "TGT"), F = c("TTC", "TTT"), L2 = c("TTA", "TTG"),
    I = c("ATA", "ATC", "ATT"),
    T = c("ACA", "ACC", "ACG", "ACT"), P = c("CCA", "CCC", "CCG", "CCT"),
    A = c("GCA", "GCC", "GCG", "GCT"), G = c("GGA", "GGC", "GGG", "GGT"),
    V = c("GTA", "GTC", "GTG", "GTT"), R = c("CGA", "CGC", "CGG", "CGT"),
    L4 = c("CTA", "CTC", "CTG", "CTT"),
    S4 = c("AGC", "AGT", "TCA", "TCC", "TCG", "TCT")
  )
  eps <- 1e-8
  F_values <- sapply(codon_table, function(codons) {
    counts <- total_counts[codons]
    counts[is.na(counts)] <- 0
    counts <- counts + pseudocount
    n <- sum(counts)
    if (n <= 1) return(NA_real_)
    p <- counts / n
    sum_f2 <- sum(p * p)
    Fcf <- (n * sum_f2 - 1) / (n - 1)
    Fcf <- max(min(Fcf, 1 - eps), eps)
    Fcf
  })
  list(
    K2 = mean(F_values[c("K","N","Q","H","D","E","Y","C","F","L2")], na.rm = TRUE),
    K3 = mean(F_values[c("I")], na.rm = TRUE),
    K4 = mean(F_values[c("T","P","A","G","V","R","L4","S4")], na.rm = TRUE)
  )
}

# Calculate gene-level EnC
calculate_enc_with_fallback <- function(seq, sample_Ks, pseudocount = 1.0) {
  codon_counts <- uco(seq, index = "eff", frame = 0)
  names(codon_counts) <- toupper(names(codon_counts))
  codon_table <- list(
    K = c("AAA", "AAG"), N = c("AAC", "AAT"), Q = c("CAA", "CAG"),
    H = c("CAC", "CAT"), D = c("GAC", "GAT"), E = c("GAA", "GAG"),
    Y = c("TAC", "TAT"), C = c("TGC", "TGT"), F = c("TTC", "TTT"), L2 = c("TTA", "TTG"),
    I = c("ATA", "ATC", "ATT"),
    T = c("ACA", "ACC", "ACG", "ACT"), P = c("CCA", "CCC", "CCG", "CCT"),
    A = c("GCA", "GCC", "GCG", "GCT"), G = c("GGA", "GGC", "GGG", "GGT"),
    V = c("GTA", "GTC", "GTG", "GTT"), R = c("CGA", "CGC", "CGG", "CGT"),
    L4 = c("CTA", "CTC", "CTG", "CTT"),
    S4 = c("AGC", "AGT", "TCA", "TCC", "TCG", "TCT")
  )
  eps <- 1e-8
  F_values <- sapply(codon_table, function(codons) {
    raw <- codon_counts[codons]
    if (all(is.na(raw)) || sum(raw, na.rm = TRUE) == 0) return(NA_real_)
    counts <- raw
    counts[is.na(counts)] <- 0
    counts <- counts + pseudocount
    n <- sum(counts)
    if (n <= 1) return(NA_real_)
    p <- counts / n
    sum_f2 <- sum(p * p)
    Fcf <- (n * sum_f2 - 1) / (n - 1)
    max(min(Fcf, 1 - eps), eps)
  })
  K2 <- mean(F_values[c("K","N","Q","H","D","E","Y","C","F","L2")], na.rm = TRUE)
  K3 <- mean(F_values[c("I")], na.rm = TRUE)
  K4 <- mean(F_values[c("T","P","A","G","V","R","L4","S4")], na.rm = TRUE)
  if (is.na(K2)) K2 <- sample_Ks$K2
  if (is.na(K3)) K3 <- sample_Ks$K3
  if (is.na(K4)) K4 <- sample_Ks$K4
  K2 <- ifelse(is.na(K2), NA_real_, pmax(K2, eps))
  K3 <- ifelse(is.na(K3), NA_real_, pmax(K3, eps))
  K4 <- ifelse(is.na(K4), NA_real_, pmax(K4, eps))
  if (any(is.na(c(K2,K3,K4)))) return(NA_real_)
  enc <- 2 + 9 / K2 + 1 / K3 + 5 / K4
  enc <- min(max(enc, 20), 61)
  enc
}

# ========================================
# Helper function for EnC/GC3 Calculation
# ========================================

compute_enc_gc3_for_one_sample <- function(sample_row, min_length_nt = 200, pseudocount = 0.5) {
  full_path <- file.path(base_dir, sample_row$file_path)
  genome <- safe_read_fasta(full_path)
  if (length(genome) == 0) {
    warning("No sequences found: ", full_path)
    return(tibble())
  }
  is_valid <- function(x) {
    len <- length(x)
    len >= min_length_nt && (len %% 3 == 0) && !any(!toupper(x) %in% c("A","C","G","T"))
  }
  genome <- genome[vapply(genome, is_valid, logical(1))]
  if (length(genome) == 0) {
    warning("All sequences filtered out: ", full_path)
    return(tibble())
  }
  sample_Ks <- compute_sample_Ks(genome, pseudocount = pseudocount)
  enc_gc3_list <- lapply(genome, function(seq) {
    enc <- calculate_enc_with_fallback(seq, sample_Ks, pseudocount = pseudocount)
    gc3 <- GC3(seq)
    c(enc, gc3)
  })
  df <- as.data.frame(do.call(rbind, enc_gc3_list))
  colnames(df) <- c("EnC","GC3")
  df$gene_id      <- names(genome)
  df$sample_id    <- sample_row$sample_id
  df$sample_type  <- sample_row$sample_type
  df$host_genus   <- sample_row$host_genus
  df$host_species <- sample_row$host_species
  df$lifestyle    <- sample_row$lifestyle
  df$trna_group   <- sample_row$trna_group
  df
}

# Match each phage to its host - species or genus level
match_host_for_phage <- function(phage_row, host_gene_tbl) {
  sp <- phage_row$host_species
  gn <- phage_row$host_genus
  by_sp <- host_gene_tbl %>% filter(host_species == sp)
  if (nrow(by_sp) > 0) return(by_sp)
  by_gn <- host_gene_tbl %>% filter(host_genus == gn)
  if (nrow(by_gn) > 0) return(by_gn)
  tibble()
}

# ============================
# Load and normalize metadata
# ============================

required_cols <- c("sample_id","file_path","sample_type","host_genus","host_species","lifestyle","trna_group")

md_raw <- readxl::read_excel(metadata_xlsx, sheet = 1)

names(md_raw) <- names(md_raw) |>
  enc2utf8() |>
  sub("^\ufeff", "", x = _) |>
  gsub("^\\s+|\\s+$", "", x = _) |>
  gsub("\\s+", "_", x = _) |>
  gsub("[^A-Za-z0-9_]+", "_", x = _) |>
  tolower()

missing <- setdiff(required_cols, names(md_raw))
if (length(missing)) stop("Metadata is missing required columns: ", paste(missing, collapse = ", "))

metadata <- md_raw %>%
  mutate(across(where(is.character), standardize_text)) %>%
  mutate(
    sample_type  = factor(sample_type, levels = c("host","phage")),
    host_genus   = as.character(host_genus),
    host_species = as.character(host_species),
    lifestyle    = ifelse(is.na(lifestyle) | lifestyle == "", NA, as.character(lifestyle)),
    trna_group   = ifelse(is.na(trna_group) | trna_group == "", NA, as.character(trna_group))
  )

# check output (optional)
print(head(metadata, 3))

# ===========================
# Compute gene-level EnC/GC3
# ===========================

# computes both EnC/GC3 correlation data as well as EnC values
all_gene <- metadata %>%
  split(.$sample_id) %>%
  purrr::map_dfr(~ compute_enc_gc3_for_one_sample(.x[1, , drop = FALSE]))

all_gene %>%
  split(.$host_genus) %>%
  purrr::iwalk(function(df, g) {
    readr::write_csv(df, file.path(out_dir, paste0("gene_level_EnC_GC3_", g, ".csv")))
  })

readr::write_csv(all_gene, file.path(out_dir, "gene_level_EnC_GC3_all_samples.csv"))

# ============================================================
# Statistical summaries - medians of all phages within a host
# ============================================================

enc_medians_by_sample <- all_gene %>%
  group_by(sample_id, sample_type, host_genus) %>%
  summarise(median_EnC = median(EnC, na.rm = TRUE), .groups = "drop")

host_enc_summary <- enc_medians_by_sample %>%
  filter(sample_type == "host") %>%
  group_by(host_genus) %>%
  summarise(host_median_EnC = median(median_EnC, na.rm = TRUE), .groups = "drop")

phage_enc_summary <- enc_medians_by_sample %>%
  filter(sample_type == "phage") %>%
  group_by(host_genus) %>%
  summarise(
    phages_combined_median_EnC = median(median_EnC, na.rm = TRUE),
    sd_of_phage_medians        = sd(median_EnC, na.rm = TRUE),
    n_phages                   = dplyr::n(),
    .groups = "drop"
  )

enc_genus_summary <- host_enc_summary %>%
  full_join(phage_enc_summary, by = "host_genus") %>%
  arrange(host_genus)

readr::write_csv(enc_genus_summary, file.path(out_sum, "EnC_genus_summary.csv"))

# =========================================
# Statistical summaries - per-sample means
# =========================================

sample_means <- all_gene %>%
  group_by(sample_id, sample_type, host_genus, host_species, lifestyle, trna_group) %>%
  summarize(
    mean_EnC = mean(EnC, na.rm = TRUE),
    mean_GC3 = mean(GC3, na.rm = TRUE),
    n_genes  = sum(!is.na(EnC)),
    .groups = "drop"
  )

readr::write_csv(sample_means, file.path(out_dir, "sample_means_EnC_GC3.csv"))

phage_genes <- all_gene %>% filter(sample_type == "phage")
host_genes  <- all_gene %>% filter(sample_type == "host")
phage_means <- sample_means %>% filter(sample_type == "phage")
host_means  <- sample_means %>% filter(sample_type == "host")

# ================================================
# Per-phage vs host (means + gene-level Wilcoxon)
# ================================================

# computes both within host genus as well as globally across all samples
per_phage_results <- phage_means %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  pmap_dfr(function(sample_id, sample_type, host_genus, host_species,
                    lifestyle, trna_group, mean_EnC, mean_GC3, n_genes) {
    pr <- tibble::tibble(
      sample_id    = sample_id,
      host_species = host_species,
      host_genus   = host_genus,
      lifestyle    = lifestyle,
      trna_group   = trna_group,
      mean_EnC     = mean_EnC
    )
    ph_g   <- phage_genes %>% filter(sample_id == !!sample_id)
    host_g <- match_host_for_phage(pr, host_genes)
    if (nrow(ph_g) < 2 || nrow(host_g) < 2) {
      return(tibble::tibble(
        phage_id       = sample_id,
        host_species   = host_species,
        host_genus     = host_genus,
        lifestyle      = lifestyle,
        trna_group     = trna_group,
        mean_EnC_phage = mean_EnC,
        mean_EnC_host  = NA_real_,
        delta_EnC      = NA_real_,
        abs_delta_EnC  = NA_real_,
        W              = NA_real_,
        p_value        = NA_real_
      ))
    }
    host_scope <- if (!is.na(host_species) && any(host_g$host_species == host_species)) {
      host_means %>% filter(host_species == !!host_species)
    } else {
      host_means %>% filter(host_genus == !!host_genus)
    }
    host_mean_enc <- mean(host_scope$mean_EnC, na.rm = TRUE)
    wt <- try(stats::wilcox.test(ph_g$EnC, host_g$EnC, exact = FALSE), silent = TRUE)
    if (inherits(wt, "try-error")) {
      Wv <- NA_real_; pv <- NA_real_
    } else {
      Wv <- unname(wt$statistic); pv <- wt$p.value
    }
    tibble::tibble(
      phage_id       = sample_id,
      host_species   = host_species,
      host_genus     = host_genus,
      lifestyle      = lifestyle,
      trna_group     = trna_group,
      mean_EnC_phage = mean_EnC,
      mean_EnC_host  = host_mean_enc,
      delta_EnC      = mean_EnC - host_mean_enc,
      abs_delta_EnC  = abs(mean_EnC - host_mean_enc),
      W              = Wv,
      p_value        = pv
    )
  }) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(host_genus, host_species, phage_id)

readr::write_csv(per_phage_results, file.path(out_sum, "per_phage_vs_host_EnC_means_genelevel_wilcoxon.csv"))
readr::write_csv(per_phage_results, file.path(out_sum, "global_per_phage_EnC_vs_host_deltas.csv"))

# ======================================================
# Per-genus Stats Functions for tRNA group or lifestyle
# ======================================================

# Compares tRNA groups: High, Low, None; with/without tRNAs, and lifestyle: Virulent, Temperate

desired_pairs_for_genus <- function(genus) {
  c("High|Low", "High|None", "Low|None")}

.safe_dunn_pairs <- function(df, value_col, group_col) {
  res <- try(FSA::dunnTest(stats::reformulate(group_col, value_col),
                           data = df, method = "bh"),
             silent = TRUE)
  if (inherits(res, "try-error")) {
    return(tibble::tibble(
      group1 = character(), group2 = character(),
      Z = numeric(), p_raw = numeric(), p_adj_BH_dunn = numeric()
    ))
  }
  out <- tibble::tibble(res$res)
  if (!nrow(out)) {
    return(tibble::tibble(
      group1 = character(), group2 = character(),
      Z = numeric(), p_raw = numeric(), p_adj_BH_dunn = numeric()
    ))
  }
  out %>%
    tidyr::separate(Comparison, into = c("group1","group2"), sep = " - ", remove = FALSE) %>%
    dplyr::mutate(group1 = trimws(group1), group2 = trimws(group2)) %>%
    dplyr::transmute(group1, group2, Z = Z, p_raw = P.unadj, p_adj_BH_dunn = P.adj)
}

.compute_pair_wilcox <- function(df, value_col, group_col, g1, g2) {
  x <- df %>% dplyr::filter(.data[[group_col]] == g1) %>% dplyr::pull(.data[[value_col]])
  y <- df %>% dplyr::filter(.data[[group_col]] == g2) %>% dplyr::pull(.data[[value_col]])
  n1 <- sum(is.finite(x)); n2 <- sum(is.finite(y))
  if (n1 == 0 || n2 == 0) return(NULL)
  wt <- try(stats::wilcox.test(x, y, exact = FALSE), silent = TRUE)
  if (inherits(wt,"try-error")) return(NULL)
  m  <- outer(x, y, "-")
  cd <- (sum(m > 0) - sum(m < 0)) / (length(x)*length(y))
  tibble::tibble(
    group1 = as.character(g1), group2 = as.character(g2),
    n1 = n1, n2 = n2,
    median1 = stats::median(x, na.rm = TRUE), median2 = stats::median(y, na.rm = TRUE),
    p_raw = wt$p.value,
    cliffs_delta = cd
  )
}

.build_complete_rows_one_genus <- function(df, genus, value_col, group_col = "trna_group") {
  df <- df %>% dplyr::mutate(trna_group = factor(.data[[group_col]], levels = c("None","Low","High")))
  want <- desired_pairs_for_genus(genus)
  present <- stats::na.omit(unique(as.character(df[[group_col]])))
  k <- length(present)
  dunn_rows <- if (k >= 3) .safe_dunn_pairs(df, value_col, group_col) else tibble::tibble(
    group1 = character(), group2 = character(),
    Z = numeric(), p_raw = numeric(), p_adj_BH_dunn = numeric()
  )
  rows <- tibble::tibble()
  for (kp in want) {
    parts <- strsplit(kp, "\\|")[[1]]
    g1 <- parts[1]; g2 <- parts[2]
    if (g1 %in% present && g2 %in% present) {
      dr <- if (nrow(dunn_rows) > 0) dplyr::filter(dunn_rows, group1 == g1, group2 == g2) else tibble::tibble()
      add <- .compute_pair_wilcox(df, value_col, group_col, g1, g2)
      if (nrow(dr) == 1) {
        row <- tibble::tibble(
          host_genus = genus,
          comparison = paste0(g1, " - ", g2),
          group1 = g1, group2 = g2,
          n1 = add$n1 %||% NA_integer_, n2 = add$n2 %||% NA_integer_,
          median1 = add$median1 %||% NA_real_, median2 = add$median2 %||% NA_real_,
          Z = dr$Z,
          p_raw = dr$p_raw,
          p_adj_BH = dr$p_adj_BH_dunn,
          cliffs_delta = add$cliffs_delta %||% NA_real_,
          note = "Dunn"
        )
      } else {
        row <- tibble::tibble(
          host_genus = genus,
          comparison = paste0(g1, " - ", g2),
          group1 = g1, group2 = g2,
          n1 = add$n1 %||% NA_integer_, n2 = add$n2 %||% NA_integer_,
          median1 = add$median1 %||% NA_real_, median2 = add$median2 %||% NA_real_,
          Z = NA_real_,
          p_raw = add$p_raw %||% NA_real_,
          p_adj_BH = NA_real_,
          cliffs_delta = add$cliffs_delta %||% NA_real_,
          note = "Wilcoxon"
        )
      }
    } else {
      row <- tibble::tibble(
        host_genus = genus,
        comparison = paste0(g1, " - ", g2),
        group1 = g1, group2 = g2,
        n1 = 0L, n2 = 0L,
        median1 = NA_real_, median2 = NA_real_,
        Z = NA_real_, p_raw = NA_real_, p_adj_BH = NA_real_,
        cliffs_delta = NA_real_, note = "group_missing"
      )
    }
    rows <- dplyr::bind_rows(rows, row)
  }
  idx <- which(is.finite(rows$p_raw) & is.na(rows$p_adj_BH))
  if (length(idx) > 0) rows$p_adj_BH[idx] <- p.adjust(rows$p_raw[idx], method = "BH")
  rows
}

# ==========================================================
# Per-genus pairwise on absolute ∆EnC comparing tRNA groups
# ==========================================================

per_phage_trna <- per_phage_results %>%
  filter(!is.na(abs_delta_EnC), !is.na(trna_group)) %>%
  mutate(trna_group = factor(trna_group, levels = c("None","Low","High")))

pairwise_trna_per_genus_absdelta_complete <- per_phage_trna %>%
  split(.$host_genus) %>%
  purrr::imap_dfr(~ .build_complete_rows_one_genus(.x, .y, value_col = "abs_delta_EnC"))

write_csv(pairwise_trna_per_genus_absdelta_complete,
          file.path(out_sum, "PER_GENUS_trna_PAIRWISE_complete_on_absDeltaEnC.csv"))

# ==========================================
# Per-genus comparing with to without tRNAs
# ==========================================

with_vs_none_per_genus <- per_phage_results %>%
  filter(!is.na(trna_group), !is.na(abs_delta_EnC)) %>%
  mutate(tRNA_presence = ifelse(trna_group == "None", "No_tRNA", "With_tRNA")) %>%
  group_by(host_genus) %>%
  group_modify(~{
    df <- .x
    if (n_distinct(df$tRNA_presence) < 2) {
      return(tibble(test="Wilcoxon", W=NA_real_, p_value=NA_real_,
                    n_With_tRNA=sum(df$tRNA_presence=="With_tRNA"),
                    n_No_tRNA=sum(df$tRNA_presence=="No_tRNA"),
                    cliffs_delta=NA_real_))
    }
    wt <- try(wilcox.test(abs_delta_EnC ~ tRNA_presence, data = df, exact = FALSE), silent = TRUE)
    if (inherits(wt, "try-error")) {
      tibble(test="Wilcoxon", W=NA_real_, p_value=NA_real_,
             n_With_tRNA=sum(df$tRNA_presence=="With_tRNA"),
             n_No_tRNA=sum(df$tRNA_presence=="No_tRNA"),
             cliffs_delta=NA_real_)
    } else {
      tibble(test="Wilcoxon",
             W=unname(wt$statistic),
             p_value=wt$p.value,
             n_With_tRNA=sum(df$tRNA_presence=="With_tRNA"),
             n_No_tRNA=sum(df$tRNA_presence=="No_tRNA"),
             cliffs_delta=cliffs_delta(df$abs_delta_EnC, df$tRNA_presence))
    }
  }) %>% ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method="BH"))

write_csv(with_vs_none_per_genus,
          file.path(out_sum, "per_genus_trna_WithVsNone_Wilcoxon_absDeltaEnC.csv"))

# ========================================
# Per-genus pairwise lifestyle comparison
# ========================================

lifestyle_within_genus_delta <- phage_means %>%
  dplyr::filter(!is.na(lifestyle), !is.na(delta_EnC), !is.na(host_genus)) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::group_modify(~{
    df <- .x
    if (dplyr::n_distinct(stats::na.omit(df$lifestyle)) < 2) {
      return(tibble::tibble(
        W = NA_real_, p_value = NA_real_,
        n_virulent  = sum(df$lifestyle == "Virulent",  na.rm = TRUE),
        n_temperate = sum(df$lifestyle == "Temperate", na.rm = TRUE),
        cliffs_delta = NA_real_,
        median_virulent  = stats::median(df$delta_EnC[df$lifestyle == "Virulent"],  na.rm = TRUE),
        median_temperate = stats::median(df$delta_EnC[df$lifestyle == "Temperate"], na.rm = TRUE)
      ))
    }
    wt <- try(stats::wilcox.test(delta_EnC ~ lifestyle, data = df, exact = FALSE), silent = TRUE)
    if (inherits(wt, "try-error")) {
      tibble::tibble(
        W = NA_real_, p_value = NA_real_,
        n_virulent  = sum(df$lifestyle == "Virulent",  na.rm = TRUE),
        n_temperate = sum(df$lifestyle == "Temperate", na.rm = TRUE),
        cliffs_delta = NA_real_,
        median_virulent  = stats::median(df$delta_EnC[df$lifestyle == "Virulent"],  na.rm = TRUE),
        median_temperate = stats::median(df$delta_EnC[df$lifestyle == "Temperate"], na.rm = TRUE)
      )
    } else {
      # Cliff's delta for Virulent vs Temperate on ΔEnC
      x <- df$delta_EnC[df$lifestyle == "Virulent"]
      y <- df$delta_EnC[df$lifestyle == "Temperate"]
      m <- outer(x, y, "-")
      cd <- (sum(m > 0, na.rm = TRUE) - sum(m < 0, na.rm = TRUE)) / (length(x) * length(y))
      
      tibble::tibble(
        W = unname(wt$statistic),
        p_value = wt$p.value,
        n_virulent  = sum(df$lifestyle == "Virulent",  na.rm = TRUE),
        n_temperate = sum(df$lifestyle == "Temperate", na.rm = TRUE),
        cliffs_delta = cd,
        median_virulent  = stats::median(x, na.rm = TRUE),
        median_temperate = stats::median(y, na.rm = TRUE)
      )
    }
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_BH = stats::p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(host_genus)

readr::write_csv(
  lifestyle_within_genus_delta,
  file.path(out_sum, "per_genus_lifestyle_Wilcoxon_deltaEnC.csv"))

# =========================================
# Global tests on tRNA group and lifestyle
# =========================================

# by tRNA group
trna_tbl <- per_phage_results %>%
  filter(!is.na(abs_delta_EnC), !is.na(trna_group)) %>%
  mutate(trna_group = factor(trna_group, levels = c("None","Low","High")))

kw_trna <- safe_kw(abs_delta_EnC ~ trna_group, trna_tbl)

kw_trna_tbl <- if (!is.null(kw_trna)) {
  tibble(
    test      = "Kruskal",
    statistic = unname(kw_trna$statistic),
    df        = unname(kw_trna$parameter),
    p_value   = kw_trna$p.value,
    n         = nrow(trna_tbl),
    eta2      = kw_eta2(kw_trna, nrow(trna_tbl))
  )
} else tibble(test="Kruskal", statistic=NA_real_, df=NA_real_, p_value=NA_real_, n=nrow(trna_tbl), eta2=NA_real_)

write_csv(kw_trna_tbl, file.path(out_sum, "GLOBAL_trna_KW_absDeltaEnC.csv"))

dunn_trna <- NULL
if (!is.null(kw_trna) && kw_trna$p.value < 0.05) {
  dunn_trna <- FSA::dunnTest(abs_delta_EnC ~ trna_group, data = trna_tbl, method = "bh")$res
  write_csv(dunn_trna, file.path(out_sum, "GLOBAL_trna_Dunn_absDeltaEnC.csv"))
}

# with tRNA vs. without tRNA
trna_tbl2 <- trna_tbl %>% mutate(tRNA_presence = ifelse(trna_group == "None", "No_tRNA", "With_tRNA"))
wilx_trna_presence <- if (n_distinct(trna_tbl2$tRNA_presence) == 2) {
  wilcox.test(abs_delta_EnC ~ tRNA_presence, data = trna_tbl2, exact = FALSE)
} else NULL

wilx_trna_presence_tbl <- if (!is.null(wilx_trna_presence)) {
  tibble(
    test          = "Wilcoxon",
    W             = unname(wilx_trna_presence$statistic),
    p_value       = wilx_trna_presence$p.value,
    n_With_tRNA   = sum(trna_tbl2$tRNA_presence == "With_tRNA"),
    n_No_tRNA     = sum(trna_tbl2$tRNA_presence == "No_tRNA"),
    cliffs_delta  = cliffs_delta(trna_tbl2$abs_delta_EnC, trna_tbl2$tRNA_presence)
  )
} else tibble(test="Wilcoxon", W=NA_real_, p_value=NA_real_, n_With_tRNA=NA_integer_, n_No_tRNA=NA_integer_, cliffs_delta=NA_real_)

write_csv(wilx_trna_presence_tbl, file.path(out_sum, "GLOBAL_trna_WithVsNone_Wilcoxon_absDeltaEnC.csv"))

# lifestyle
life_tbl <- per_phage_results %>%
  filter(!is.na(lifestyle), !is.na(abs_delta_EnC)) %>%
  mutate(lifestyle = factor(lifestyle, levels = c("Virulent","Temperate")))

wilx_life <- if (n_distinct(life_tbl$lifestyle) == 2) {
  wilcox.test(abs_delta_EnC ~ lifestyle, data = life_tbl, exact = FALSE)
} else NULL

wilx_life_tbl <- if (!is.null(wilx_life)) {
  tibble(
    test        = "Wilcoxon",
    W           = unname(wilx_life$statistic),
    p_value     = wilx_life$p.value,
    n_virulent  = sum(life_tbl$lifestyle == "Virulent", na.rm = TRUE),
    n_temperate = sum(life_tbl$lifestyle == "Temperate", na.rm = TRUE),
    cliffs_delta= cliffs_delta(life_tbl$abs_delta_EnC, life_tbl$lifestyle)
  )
} else tibble(test="Wilcoxon", W=NA_real_, p_value=NA_real_, n_virulent=NA_integer_, n_temperate=NA_integer_, cliffs_delta=NA_real_)

write_csv(wilx_life_tbl, file.path(out_sum, "GLOBAL_lifestyle_Wilcoxon_absDeltaEnC.csv"))

# ============================================
# Optional: Excel workbook with all summaries
# ============================================
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
wb <- openxlsx::createWorkbook()

add_sheet <- function(name, df) {
  if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, df)
  }
}

add_sheet("EnC_genus_summary", enc_genus_summary)
add_sheet("PerPhage_vsHost", per_phage_results)
add_sheet("PerGenus_trna_pairs_absD", pairwise_trna_per_genus_absdelta_complete)
add_sheet("PerGenus_WithVsNone_absD", with_vs_none_per_genus)
add_sheet("PerGenus_Lifestyle_meanEnC", lifestyle_within_genus)
add_sheet("GLOBAL_KW_absD", kw_trna_tbl)
add_sheet("GLOBAL_Dunn_absD", if (!is.null(dunn_trna)) dunn_trna else NULL)
add_sheet("GLOBAL_WithVsNone_absD", wilx_trna_presence_tbl)
add_sheet("GLOBAL_Lifestyle_absD", wilx_life_tbl)

openxlsx::saveWorkbook(wb, file.path(out_sum, "EnC_tests_summary.xlsx"), overwrite = TRUE)

# =========================================
# Nc plots (faceted by host genus) EnC/GC3
# =========================================

# set desired host order
host_order <- c("Host1","Host2","...","HostN")

# clean gene table for plotting
gene_clean <- all_gene %>%
  filter(is.finite(EnC), is.finite(GC3), GC3 >= 0, GC3 <= 1) %>%
  mutate(
    sample_type = factor(as.character(sample_type), levels = c("host","phage")),
    host_genus  = factor(as.character(host_genus), levels = host_order)
  )

phage_only <- gene_clean %>% filter(sample_type == "phage")

# Expected Nc curve (Wright 1990
nc_expected <- function(gc3) 2 + gc3 + 29 / (gc3^2 + (1 - gc3)^2)
curve_df <- data.frame(GC3 = seq(0, 1, length.out = 500)) %>%
  mutate(Nc_expected = nc_expected(GC3))

# =========================
# Hosts vs Phages by genus
# =========================
p_hosts_phages_by_genus <- ggplot(gene_clean, aes(x = GC3, y = EnC)) +
  geom_point(aes(color = sample_type, shape = sample_type),
             alpha = 0.35, size = 1.4) +
  geom_line(data = curve_df,
            aes(x = GC3, y = Nc_expected),
            inherit.aes = FALSE,
            linewidth = 0.8) +
  scale_color_manual(values = c(host = "#2C7BB6", phage = "#D7191C")) +
  scale_shape_manual(values = c(host = 16, phage = 17)) +
  labs(x = "GC3 (GC at 3rd codon position)", y = "EnC (Nc)",
       color = "", shape = "") +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(8, "pt"),
    panel.grid = element_blank()   # <<< REMOVE GRIDLINES
  ) +
  facet_wrap(~ host_genus, scales = "free", drop = TRUE)

ggsave(file.path(out_dir, "Ncplot_BY_GENUS_hosts_vs_phages.png"),
       p_hosts_phages_by_genus, width = 12, height = 9, dpi = 300)

# ==========================
# Phage genes by tRNA group
# ==========================
ph_trna <- phage_only %>%
  filter(!is.na(trna_group)) %>%
  mutate(trna_group = factor(trna_group, levels = c("None","Low","High")))

if (nrow(ph_trna) > 0) {
  p_trna_by_genus <- ggplot(ph_trna, aes(x = GC3, y = EnC, color = trna_group)) +
    geom_point(alpha = 0.45, size = 1.3) +
    geom_line(data = curve_df,
              aes(x = GC3, y = Nc_expected),
              inherit.aes = FALSE,
              linewidth = 0.8) +
    scale_color_manual(values = c(
      None = "#5e3c99",
      Low  = "#b2abd2",
      High = "#fdb863"
    )) +
    labs(x = "GC3", y = "EnC (Nc)", color = "tRNA group") +
    theme_bw(base_size = 16) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(8, "pt"),
      panel.grid = element_blank()   # <<< REMOVE GRIDLINES
    ) +
    facet_wrap(~ host_genus, scales = "free", drop = TRUE)
  
  ggsave(file.path(out_dir, "Ncplot_PHAGE_by_tRNAgroup_BY_GENUS.png"),
         p_trna_by_genus, width = 12, height = 9, dpi = 300)
}

# =========================
# Phage genes by lifestyle
# =========================
ph_life <- phage_only %>%
  filter(!is.na(lifestyle)) %>%
  mutate(lifestyle = factor(lifestyle, levels = c("Virulent","Temperate")))

if (nrow(ph_life) > 0) {
  p_life_by_genus <- ggplot(ph_life, aes(x = GC3, y = EnC, color = lifestyle)) +
    geom_point(alpha = 0.45, size = 1.3) +
    geom_line(data = curve_df,
              aes(x = GC3, y = Nc_expected),
              inherit.aes = FALSE,
              linewidth = 0.8) +
    scale_color_manual(values = c(
      Virulent  = "#1f78b4",
      Temperate = "#33a02c"
    )) +
    labs(x = "GC3", y = "EnC (Nc)", color = "Lifestyle") +
    theme_bw(base_size = 16) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(8, "pt"),
      panel.grid = element_blank()   # <<< REMOVE GRIDLINES
    ) +
    facet_wrap(~ host_genus, scales = "free", drop = TRUE)
  
  ggsave(file.path(out_dir, "Ncplot_PHAGE_by_lifestyle_BY_GENUS.png"),
         p_life_by_genus, width = 12, height = 9, dpi = 300)
}

# =============================
# Spearman correlation EnC/GC3
# =============================

spearman_grouped <- function(df, group_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    group_modify(~{
      dd <- .x
      if (nrow(dd) < 3) {
        tibble(n = nrow(dd), rho = NA_real_, p_value = NA_real_)
      } else {
        ct <- suppressWarnings(cor.test(dd$GC3, dd$EnC, method = "spearman"))
        tibble(n = nrow(dd), rho = unname(ct$estimate), p_value = ct$p.value)
      }
    }) %>%
    ungroup() %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
}

corr_by_genus_host_vs_phage <- gene_clean %>%
  spearman_grouped(group_cols = c("host_genus", "sample_type"))
write_csv(corr_by_genus_host_vs_phage, file.path(out_sum, "CORR_byGenus_ALL_host_vs_phage.csv"))

corr_by_genus_phage_only <- gene_clean %>%
  filter(sample_type == "phage") %>%
  spearman_grouped(group_cols = c("host_genus"))
write_csv(corr_by_genus_phage_only, file.path(out_sum, "CORR_byGenus_PHAGE_only.csv"))

corr_by_genus_trna <- gene_clean %>%
  filter(sample_type == "phage", !is.na(trna_group)) %>%
  mutate(trna_group = factor(trna_group, levels = c("None","Low","High"))) %>%
  spearman_grouped(group_cols = c("host_genus", "trna_group"))
write_csv(corr_by_genus_trna, file.path(out_sum, "CORR_byGenus_PHAGE_by_tRNAgroup.csv"))

corr_by_genus_life <- gene_clean %>%
  filter(sample_type == "phage", !is.na(lifestyle)) %>%
  mutate(lifestyle = factor(lifestyle, levels = c("Virulent","Temperate"))) %>%
  spearman_grouped(group_cols = c("host_genus", "lifestyle"))
write_csv(corr_by_genus_life, file.path(out_sum, "CORR_byGenus_PHAGE_by_lifestyle.csv"))

# ================================
# Optional: Quick console summary
# ================================
cat("\n=== SUMMARY ===\n")
cat("Samples in metadata:", nrow(metadata), "\n")
cat("Gene rows computed:", nrow(all_gene), "\n")
cat("Phages evaluated (per_phage_results):", nrow(per_phage_results), "\n")
if (!is.null(kw_trna)) cat("GLOBAL KW(trna_group) p =", signif(kw_trna$p.value, 3),
                           " eta2=", signif(kw_eta2(kw_trna, nrow(trna_tbl)), 3), "\n")
cat("Outputs in:", out_dir, "\n")
