# ===========================
# CALCULATING RSCU and ∆RSCU
# ===========================

# This script is designed with bacterial hosts and bacteriophages in mind
# it can be applied to other virus/host pairs :)

# Set base and output directories
# Load metadata - should have columns "sample_id" (organism name), "file_path" (file path for the specific FASTA sequence),
#                 "sample_type" (i.e. phage or host), "host_genus", "host_species", "lifestyle" (Virulent or Temperate),
#                 and "trna_group" (High, Low, or None)

base_dir      <- "/your/file/path/FASTAs"   # assumes FASTA files are all saved to a single parent folder
metadata_xlsx <- "/your/file/path/metadata.xlsx" 
out_dir       <- "/your/file/path/Results/RSCU"

# Load required packages

req_pkgs <- c("dplyr","tidyr","purrr","readr", "seqinr","stringr","forcats",
              "ggplot2","pheatmap","broom","vegan","coin","effsize","FSA")

not_installed <- req_pkgs[!req_pkgs %in% installed.packages()[,"Package"]]
if(length(not_installed)) install.packages(not_installed, dependencies = TRUE)
invisible(lapply(req_pkgs, library, character.only = TRUE))

theme_set(theme_minimal(base_size = 14))

# ================================
# SETTINGS - can tweak as desired
# ================================

# Factor orders
lifestyle_order <- c("virulent","temperate")
trna_order      <- c("none","low","high")

# Palettes for plotting
palette_lifestyle <- c(virulent = "#1f78b4", temperate = "#33a02c")
palette_trna      <- c(none = "#5e3c99", low = "#b2abd2", high = "#fdb863")

# RSCU options
include_stop_codons <- FALSE
rscu_method <- "by_counts"  # "by_counts" or "by_gene_mean"

# Codon order
sense_codons <- c(
  "ttt","ttc","tta","ttg","ctt","ctc","cta","ctg",
  "att","atc","ata","atg",
  "gtt","gtc","gta","gtg",
  "tct","tcc","tca","tcg","agt","agc",
  "cct","ccc","cca","ccg",
  "act","acc","aca","acg",
  "gct","gcc","gca","gcg",
  "tat","tac","cat","cac","caa","cag","aat","aac",
  "aaa","aag","gat","gac","gaa","gag",
  "tgt","tgc","tgg",
  "cgt","cgc","cga","cgg","aga","agg",
  "ggt","ggc","gga","ggg"
)
stop_codons <- c("taa","tag","tga")

# Optional drop Met/Trp
excluded_codons <- c("atg","tgg")
filtered_codons <- setdiff(sense_codons, excluded_codons)

# =================================================
# Helper functions for RSCU and ∆RSCU Calculations
# =================================================

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste(...), "\n"))

is_absolute_path <- function(p) grepl("^/", p) || grepl("^[A-Za-z]:[/\\]", p)
resolve_path <- function(p){
  if(is.na(p) || p == "") return(p)
  if(is_absolute_path(p) || base_dir == "") return(p)
  file.path(base_dir, p)
}
check_paths <- function(df){
  missing <- df |> dplyr::filter(!file.exists(.data$file_path))
  if(nrow(missing)){
    warning("These file paths do not exist:\n", paste(missing$file_path, collapse = "\n"))
  }
  invisible(df)
}

uco_to_tidy <- function(uco_named_vec, include_stop = FALSE){
  v <- uco_named_vec
  nm <- tolower(names(v)); names(v) <- nm
  full_codons <- if(include_stop) c(sense_codons, stop_codons) else sense_codons
  tibble::tibble(codon = full_codons) |>
    dplyr::left_join(tibble::tibble(codon = nm, rscu = as.numeric(v)), by = "codon") |>
    dplyr::mutate(rscu = dplyr::coalesce(.data$rscu, 0))
}

compute_rscu_by_counts <- function(fasta_list){
  if(length(fasta_list) == 0) return(uco_to_tidy(rep(NA_real_, 64)))
  nuc <- unlist(fasta_list, use.names = FALSE)
  nuc <- nuc[nuc %in% c("a","t","g","c","A","T","G","C")]
  r <- seqinr::uco(nuc, index = "rscu")
  uco_to_tidy(r, include_stop = include_stop_codons)
}

compute_rscu_by_gene_mean <- function(fasta_list){
  if(length(fasta_list) == 0) return(uco_to_tidy(rep(NA_real_, 64)))
  per_gene <- purrr::map(fasta_list, function(seq){
    seq <- seq[seq %in% c("a","t","g","c","A","T","G","C")]
    if(length(seq) < 3) return(uco_to_tidy(rep(NA_real_, 64)))
    r <- seqinr::uco(seq, index = "rscu")
    uco_to_tidy(r, include_stop = include_stop_codons)
  })
  dplyr::bind_rows(per_gene, .id = "gene") |>
    dplyr::group_by(.data$codon) |>
    dplyr::summarise(rscu = mean(.data$rscu, na.rm = TRUE), .groups = "drop")
}

compute_rscu_for_file <- function(file){
  fasta <- seqinr::read.fasta(file = file, as.string = FALSE, seqtype = "DNA")
  if(identical(rscu_method, "by_gene_mean")) compute_rscu_by_gene_mean(fasta) else compute_rscu_by_counts(fasta)
}

# Clean label helpers
cap_title <- function(x) tools::toTitleCase(gsub("_", " ", x))
cap_first <- function(x) paste0(toupper(substr(x,1,1)), substr(x,2,nchar(x)))

# ===============================
# Load metadata and compute RSCU
# ===============================

msg("Loading metadata …")
meta <- readr::read_csv(metadata_csv, show_col_types = FALSE) |>
  dplyr::mutate(file_path = vapply(.data$file_path, resolve_path, character(1))) |>
  dplyr::mutate(
    sample_type = tolower(.data$sample_type),
    lifestyle   = dplyr::if_else(is.na(.data$lifestyle), NA_character_,
                                 stringr::str_trim(stringr::str_to_lower(as.character(.data$lifestyle)))),
    trna_group  = dplyr::if_else(is.na(.data$trna_group), NA_character_,
                                 stringr::str_trim(stringr::str_to_lower(as.character(.data$trna_group)))),
    lifestyle   = forcats::fct_relevel(.data$lifestyle, lifestyle_order),
    trna_group  = forcats::fct_relevel(.data$trna_group, trna_order)
  )

if(!"sample_name" %in% names(meta)) meta$sample_name <- meta$sample_id
meta <- meta |> dplyr::mutate(sample_name = dplyr::coalesce(.data$sample_name, .data$sample_id))

check_paths(meta)

msg("Computing RSCU for ", nrow(meta), " samples …")
rscu_list <- purrr::map2(meta$file_path, meta$sample_id, function(fp, sid){
  tibble::tibble(sample_id = sid) |>
    dplyr::bind_cols(compute_rscu_for_file(fp))
})

rscu_tidy <- dplyr::bind_rows(rscu_list) |>
  dplyr::left_join(meta, by = "sample_id") |>
  dplyr::rename(RSCU = rscu)

readr::write_csv(rscu_tidy, file.path(out_dir, "rscu_tidy_per_sample.csv"))
msg("Saved: rscu_tidy_per_sample.csv")

# Wide matrix
rscu_wide <- rscu_tidy |>
  dplyr::select(sample_id, codon, RSCU) |>
  tidyr::pivot_wider(names_from = codon, values_from = RSCU, values_fill = list(RSCU = 0)) |>
  dplyr::arrange(sample_id)
readr::write_csv(rscu_wide, file.path(out_dir, "rscu_wide_per_sample.csv"))

# Filtered codons to drop ATG/TGG
rscu_filt <- rscu_tidy |> dplyr::filter(.data$codon %in% filtered_codons)

# ===============================================================
# Calculating ∆RSCU (phage - host) (calculates host species mean)
# ===============================================================

msg("Computing host species means and ΔRSCU …")

host_means_species <- rscu_tidy |>
  dplyr::filter(.data$sample_type == "host") |>
  dplyr::group_by(.data$host_genus, .data$host_species, .data$codon) |>
  dplyr::summarise(host_RSCU = mean(.data$RSCU, na.rm = TRUE), .groups = "drop")

phage_delta <- rscu_tidy |>
  dplyr::filter(.data$sample_type == "phage") |>
  dplyr::left_join(host_means_species, by = c("host_genus","host_species","codon")) |>
  dplyr::mutate(Delta_RSCU = .data$RSCU - .data$host_RSCU)

readr::write_csv(phage_delta, file.path(out_dir, "delta_rscu_per_phage_per_codon.csv"))

# Filtered codons to drop ATG/TGG
phage_delta_filt <- phage_delta |> dplyr::filter(.data$codon %in% filtered_codons)

# ==============================================
# Per-phage ∆RSCU summaries - MAE and Euclidian
# ==============================================

phage_delta_scores <- phage_delta_filt |>
  dplyr::group_by(.data$sample_id, .data$sample_name, .data$host_genus, .data$host_species,
                  .data$lifestyle, .data$trna_group) |>
  dplyr::summarise(
    delta_MAE = mean(abs(.data$Delta_RSCU), na.rm = TRUE),
    delta_EUC = sqrt(sum((.data$Delta_RSCU)^2, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    lifestyle  = forcats::fct_relevel(.data$lifestyle, lifestyle_order),
    trna_group = forcats::fct_relevel(.data$trna_group, trna_order)
  )

readr::write_csv(phage_delta_scores, file.path(out_dir, "delta_rscu_scores_per_phage.csv"))

# ==========================================
# Statistics within genus and across genera
# ==========================================

# tRNA groups within genus
kw_within_genus <- phage_delta_scores |>
  dplyr::filter(!is.na(.data$trna_group)) |>
  dplyr::group_by(.data$host_genus) |>
  dplyr::group_modify(~{
    df <- .x
    if(length(unique(na.omit(df$trna_group))) < 2) return(tibble::tibble(p.value = NA_real_, n = nrow(df)))
    kt <- tryCatch(stats::kruskal.test(delta_MAE ~ trna_group, data = df), error = function(e) NULL)
    if(is.null(kt)) tibble::tibble(p.value = NA_real_, n = nrow(df))
    else broom::tidy(kt) |> dplyr::transmute(p.value, n = nrow(df))
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(p.adj = p.adjust(.data$p.value, method = "BH"))

readr::write_csv(kw_within_genus, file.path(out_dir, "stats_within_genus_kw_trnagroup.csv"))

# lifestyle within genus
wilcox_within_genus <- phage_delta_scores %>%
  dplyr::filter(!is.na(lifestyle)) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::group_modify(~{
    df <- .x
    if (length(unique(na.omit(df$lifestyle))) < 2) {
      return(tibble::tibble(p.value = NA_real_, method = NA_character_, n = nrow(df), cliffs_delta = NA_real_))
    }
    tab <- table(df$lifestyle)
    uniq_rate <- dplyr::n_distinct(df$delta_MAE) / nrow(df)
    use_perm <- (min(tab) < 10) || (max(tab) / min(tab) > 3) || (uniq_rate < 0.8)
    if (use_perm) {
      set.seed(123)
      wt <- coin::wilcox_test(delta_MAE ~ lifestyle, data = df,
                              distribution = coin::approximate(nresample = 10000))
      p  <- as.numeric(coin::pvalue(wt))
      method_used <- "Permutation (coin, 10000)"
    } else {
      wt <- stats::wilcox.test(delta_MAE ~ lifestyle, data = df, exact = FALSE, correct = TRUE)
      p  <- wt$p.value
      method_used <- "Normal Approx."
    }
    cd <- suppressWarnings(
      tryCatch(effsize::cliff.delta(delta_MAE ~ lifestyle, data = df)$estimate,
               error = function(e) NA_real_)
    )
    tibble::tibble(p.value = p, method = method_used, n = nrow(df), cliffs_delta = cd)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

readr::write_csv(wilcox_within_genus, file.path(out_dir, "stats_within_genus_wilcox_lifestyle.csv"))

# across genera tests on tRNA groups and lifestyle
kw_across <- phage_delta_scores |>
  dplyr::filter(!is.na(.data$trna_group)) |>
  (\(df) broom::tidy(stats::kruskal.test(delta_MAE ~ trna_group, data = df)))()
readr::write_csv(kw_across, file.path(out_dir, "stats_across_kw_trnagroup.csv"))

wilcox_across <- phage_delta_scores |>
dplyr::filter(!is.na(.data$lifestyle)) |>
  (\(df) broom::tidy(stats::wilcox.test(delta_MAE ~ lifestyle, data = df, exact = FALSE, correct = TRUE)))()
readr::write_csv(wilcox_across, file.path(out_dir, "stats_across_wilcox_lifestyle.csv"))

# =================================
# Per-codon Helpers and Statistics
# =================================

# use filtered codons for all per-codon outputs
stopifnot(exists("phage_delta_filt"))

percodon_raw <- phage_delta_filt |>
  dplyr::filter(sample_type == "phage") |>
  dplyr::select(sample_id, sample_name, host_genus, host_species,
                lifestyle, trna_group, codon, RSCU, host_RSCU, Delta_RSCU) |>
  dplyr::arrange(host_genus, host_species, sample_id, codon)

readr::write_csv(percodon_raw, file.path(out_dir, "percodon_phage_vs_host_RAW.csv"))

safe_sign_p <- function(x, nresample = 10000) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (all(x == 0))    return(1)
  
  df <- data.frame(v = x)
  set.seed(123)
  p <- tryCatch(
    {
      wt <- coin::wilcoxsign_test(
        v ~ 1, data = df,
        distribution = coin::approximate(nresample = nresample)
      )
      as.numeric(coin::pvalue(wt))
    },
    error = function(e) {
      suppressWarnings(stats::wilcox.test(x, mu = 0, exact = FALSE, correct = TRUE)$p.value)
    }
  )
  p
}

# Per-codon phage vs host within genus
percodon_phage_vs_host_within_genus <- phage_delta_filt |>
  dplyr::filter(sample_type == "phage") |>
  dplyr::group_by(host_genus, codon) |>
  dplyr::summarise(
    n            = sum(is.finite(Delta_RSCU)),
    median_delta = stats::median(Delta_RSCU, na.rm = TRUE),
    mean_delta   = mean(Delta_RSCU, na.rm = TRUE),
    mad_delta    = stats::mad(Delta_RSCU, na.rm = TRUE),
    p_value      = safe_sign_p(Delta_RSCU, nresample = 10000),
    .groups = "drop"
  ) |>
  dplyr::group_by(host_genus) |>
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::arrange(host_genus, codon)

readr::write_csv(
  percodon_phage_vs_host_within_genus,
  file.path(out_dir, "percodon_phage_vs_host_WITHIN_GENUS_significance.csv")
)

# Per-codon within genus - tRNA group
percodon_trna <- phage_delta_filt |>
  dplyr::filter(sample_type == "phage", !is.na(trna_group)) |>
  dplyr::group_by(host_genus, codon) |>
  dplyr::group_modify(~{
    df <- .x
    if (length(unique(na.omit(df$trna_group))) < 2) return(tibble::tibble(
      n_none = sum(df$trna_group == "none", na.rm = TRUE),
      n_low  = sum(df$trna_group == "low",  na.rm = TRUE),
      n_high = sum(df$trna_group == "high", na.rm = TRUE),
      med_none = NA_real_, med_low = NA_real_, med_high = NA_real_,
      eta2H = NA_real_, p.value = NA_real_
    ))
    
    kt <- tryCatch(stats::kruskal.test(Delta_RSCU ~ trna_group, data = df),
                   error = function(e) NULL)
    if (is.null(kt)) return(tibble::tibble(
      n_none = sum(df$trna_group == "none", na.rm = TRUE),
      n_low  = sum(df$trna_group == "low",  na.rm = TRUE),
      n_high = sum(df$trna_group == "high", na.rm = TRUE),
      med_none = stats::median(df$Delta_RSCU[df$trna_group=="none"], na.rm = TRUE),
      med_low  = stats::median(df$Delta_RSCU[df$trna_group=="low"],  na.rm = TRUE),
      med_high = stats::median(df$Delta_RSCU[df$trna_group=="high"], na.rm = TRUE),
      eta2H = NA_real_, p.value = NA_real_
    ))
    
    H <- unname(kt$statistic)
    N <- sum(stats::complete.cases(df$Delta_RSCU, df$trna_group))
    eta2H <- as.numeric(H / (N - 1))
    
    tibble::tibble(
      n_none = sum(df$trna_group == "none", na.rm = TRUE),
      n_low  = sum(df$trna_group == "low",  na.rm = TRUE),
      n_high = sum(df$trna_group == "high", na.rm = TRUE),
      med_none = stats::median(df$Delta_RSCU[df$trna_group=="none"], na.rm = TRUE),
      med_low  = stats::median(df$Delta_RSCU[df$trna_group=="low"],  na.rm = TRUE),
      med_high = stats::median(df$Delta_RSCU[df$trna_group=="high"], na.rm = TRUE),
      eta2H = eta2H,
      p.value = broom::tidy(kt)$p.value
    )
  }) |>
  dplyr::group_by(host_genus) |>
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::arrange(host_genus, codon)

readr::write_csv(
  percodon_trna,
  file.path(out_dir, "percodon_kw_within_genus_trnagroup_enhanced.csv")
)

# Per-codon within genus - lifestyle
percodon_lifestyle <- phage_delta_filt |>
  dplyr::filter(sample_type == "phage", !is.na(lifestyle)) |>
  dplyr::group_by(host_genus, codon) |>
  dplyr::group_modify(~{
    df <- .x
    lvls <- na.omit(df$lifestyle)
    if (length(unique(lvls)) < 2) return(tibble::tibble(
      n_virulent = sum(df$lifestyle == "virulent", na.rm = TRUE),
      n_temperate = sum(df$lifestyle == "temperate", na.rm = TRUE),
      med_virulent = NA_real_, med_temperate = NA_real_,
      cliffs_delta = NA_real_, p.value = NA_real_, p_perm = NA_real_
    ))
    
    v  <- df$Delta_RSCU[df$lifestyle == "virulent"]
    t  <- df$Delta_RSCU[df$lifestyle == "temperate"]
    med_v <- stats::median(v, na.rm = TRUE)
    med_t <- stats::median(t, na.rm = TRUE)
    
    wt <- tryCatch(stats::wilcox.test(Delta_RSCU ~ lifestyle, data = df, exact = FALSE),
                   error = function(e) NULL)
    p_w <- if (is.null(wt)) NA_real_ else broom::tidy(wt)$p.value
    
    set.seed(123)
    p_perm <- tryCatch({
      as.numeric(coin::pvalue(
        coin::wilcox_test(Delta_RSCU ~ lifestyle, data = df,
                          distribution = coin::approximate(nresample = 10000))
      ))
    }, error = function(e) NA_real_)
    
    cd <- tryCatch(effsize::cliff.delta(v, t)$estimate, error = function(e) NA_real_)
    
    tibble::tibble(
      n_virulent = sum(is.finite(v)),
      n_temperate = sum(is.finite(t)),
      med_virulent = med_v, med_temperate = med_t,
      cliffs_delta = cd,
      p.value = p_w,
      p_perm = p_perm
    )
  }) |>
  dplyr::group_by(host_genus) |>
  dplyr::mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    p_perm.adj = p.adjust(p_perm, method = "BH")
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(host_genus, codon)

readr::write_csv(
  percodon_lifestyle,
  file.path(out_dir, "percodon_wilcox_within_genus_lifestyle_enhanced.csv")
)

# Per-codon global - tRNA group
percodon_trna_across <- phage_delta_filt |>
  dplyr::filter(sample_type == "phage", !is.na(trna_group)) |>
  dplyr::group_by(codon) |>
  dplyr::summarise(
    n_none = sum(trna_group=="none" & is.finite(Delta_RSCU)),
    n_low  = sum(trna_group=="low"  & is.finite(Delta_RSCU)),
    n_high = sum(trna_group=="high" & is.finite(Delta_RSCU)),
    med_none = stats::median(Delta_RSCU[trna_group=="none"], na.rm = TRUE),
    med_low  = stats::median(Delta_RSCU[trna_group=="low"],  na.rm = TRUE),
    med_high = stats::median(Delta_RSCU[trna_group=="high"], na.rm = TRUE),
    KW_p = tryCatch(stats::kruskal.test(Delta_RSCU ~ trna_group)$p.value,
                    error = function(e) NA_real_),
    .groups = "drop"
  ) |>
  dplyr::mutate(p_adj_BH = p.adjust(KW_p, method = "BH")) |>
  dplyr::arrange(codon)

readr::write_csv(
  percodon_trna_across,
  file.path(out_dir, "percodon_tRNA_ACROSS_GENERA_kruskal.csv")
)

# Per-codon global - lifestyle
percodon_lifestyle_across <- phage_delta_filt |>
  dplyr::filter(sample_type == "phage", !is.na(lifestyle)) |>
  dplyr::group_by(codon) |>
  dplyr::summarise(
    n_vir = sum(lifestyle == "virulent"  & is.finite(Delta_RSCU)),
    n_tmp = sum(lifestyle == "temperate" & is.finite(Delta_RSCU)),
    med_vir = stats::median(Delta_RSCU[lifestyle=="virulent"],  na.rm = TRUE),
    med_tmp = stats::median(Delta_RSCU[lifestyle=="temperate"], na.rm = TRUE),
    diff_med = med_vir - med_tmp,
    W_p = tryCatch(stats::wilcox.test(Delta_RSCU ~ lifestyle, exact = FALSE)$p.value,
                   error = function(e) NA_real_),
    .groups = "drop"
  ) |>
  dplyr::mutate(p_adj_BH = p.adjust(W_p, method = "BH")) |>
  dplyr::arrange(codon)

readr::write_csv(
  percodon_lifestyle_across,
  file.path(out_dir, "percodon_lifestyle_ACROSS_GENERA_wilcox.csv")
)

# ==================================================================================
# Post-hoc tRNA tests per codon: Dunn (within & across) + binary (has tRNA vs none)
# ==================================================================================

# Helper functions
.safe_dunn <- function(df, adjust = "bh"){
  df <- df %>% dplyr::filter(is.finite(Delta_RSCU), !is.na(trna_group))
  if (dplyr::n_distinct(df$trna_group) < 2) {
    return(tibble::tibble(Comparison = character(), Z = numeric(), P.unadj = numeric(), P.adj = numeric()))
  }
  out <- tryCatch({
    FSA::dunnTest(Delta_RSCU ~ trna_group, data = df, method = adjust)$res
  }, error = function(e) {
    tibble::tibble(Comparison = character(), Z = numeric(), P.unadj = numeric(), P.adj = numeric())
  })
  tibble::as_tibble(out)
}

.counts_by_group <- function(df){
  df %>%
    dplyr::filter(is.finite(Delta_RSCU), !is.na(trna_group)) %>%
    dplyr::group_by(trna_group) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
}

.tidy_dunn_names <- function(dres){
  dres <- dplyr::as_tibble(dres)
  dres <- dres %>%
    dplyr::rename_with(~"Comparison", dplyr::matches("^comparison$", ignore.case = TRUE)) %>%
    dplyr::rename_with(~"Z",          dplyr::matches("^z($|\\b|\\.)", ignore.case = TRUE)) %>%
    dplyr::rename_with(~"P.unadj",    dplyr::matches("^p\\.?unadj",   ignore.case = TRUE)) %>%
    dplyr::rename_with(~"P.adj",      dplyr::matches("^p\\.?adj",     ignore.case = TRUE))
  dres
}

.split_comp_separate <- function(comp_vec){
  tidyr::separate(
    tibble::tibble(Comparison = as.character(comp_vec)),
    col = "Comparison",
    into = c("grp1","grp2"),
    sep = "\\s*-\\s*",
    remove = FALSE,
    fill = "right"
  )
}

.safe_binary <- function(df){
  df <- df %>%
    dplyr::filter(is.finite(Delta_RSCU), !is.na(trna_group)) %>%
    dplyr::mutate(trna_2 = dplyr::if_else(trna_group == "none", "none", "has_tRNA"))
  if (dplyr::n_distinct(df$trna_2) < 2) {
    return(tibble::tibble(n_none = sum(df$trna_2=="none"),
                          n_has  = sum(df$trna_2=="has_tRNA"),
                          W_stat = NA_real_, p_value = NA_real_))
  }
  wt <- tryCatch(stats::wilcox.test(Delta_RSCU ~ trna_2, data = df, exact = FALSE),
                 error = function(e) NULL)
  tibble::tibble(
    n_none = sum(df$trna_2=="none"),
    n_has  = sum(df$trna_2=="has_tRNA"),
    W_stat = if (is.null(wt)) NA_real_ else unname(wt$statistic),
    p_value= if (is.null(wt)) NA_real_ else wt$p.value
  )
}

# Dunn within genus
dunn_within <- phage_delta_filt %>%
  dplyr::filter(sample_type == "phage", !is.na(trna_group)) %>%
  dplyr::group_by(host_genus, codon) %>%
  dplyr::group_modify(~{
    df <- .x
    dres <- .safe_dunn(df, adjust = "bh")
    if (nrow(dres) == 0) return(tibble::tibble())
    
    dres <- .tidy_dunn_names(dres)
    split <- .split_comp_separate(dres$Comparison)
    counts <- .counts_by_group(df)
    
    dplyr::bind_cols(dres, split) %>%
      dplyr::left_join(counts %>% dplyr::rename(grp1 = trna_group, n1 = n), by = "grp1") %>%
      dplyr::left_join(counts %>% dplyr::rename(grp2 = trna_group, n2 = n), by = "grp2") %>%
      dplyr::select(grp1, n1, grp2, n2, Z, P.unadj, P.adj) %>%
      dplyr::rename(z_stat = Z, p_unadj = P.unadj, p_adj_BH = P.adj)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(host_genus, codon) %>%
  dplyr::arrange(host_genus, codon, grp1, grp2)

readr::write_csv(dunn_within, file.path(out_dir, "percodon_tRNA_WITHIN_GENUS_Dunn_pairwise.csv"))

# Dunn across genera (global)
dunn_across <- phage_delta_filt %>%
  dplyr::filter(sample_type == "phage", !is.na(trna_group)) %>%
  dplyr::group_by(codon) %>%
  dplyr::group_modify(~{
    df <- .x
    dres <- .safe_dunn(df, adjust = "bh")
    if (nrow(dres) == 0) return(tibble::tibble())
    
    dres <- .tidy_dunn_names(dres)
    split <- .split_comp_separate(dres$Comparison)
    counts <- .counts_by_group(df)
    
    dplyr::bind_cols(dres, split) %>%
      dplyr::mutate(scope = "ACROSS") %>%
      dplyr::left_join(counts %>% dplyr::rename(grp1 = trna_group, n1 = n), by = "grp1") %>%
      dplyr::left_join(counts %>% dplyr::rename(grp2 = trna_group, n2 = n), by = "grp2") %>%
      dplyr::select(scope, grp1, n1, grp2, n2, Z, P.unadj, P.adj) %>%
      dplyr::rename(z_stat = Z, p_unadj = P.unadj, p_adj_BH = P.adj)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(codon, .before = scope) %>%
  dplyr::arrange(codon, grp1, grp2)

readr::write_csv(dunn_across, file.path(out_dir, "percodon_tRNA_ACROSS_GENERA_Dunn_pairwise.csv"))

# Binary within genus - with vs. without tRNAs
bin_within <- phage_delta_filt %>%
  dplyr::filter(sample_type == "phage", !is.na(trna_group)) %>%
  dplyr::group_by(host_genus, codon) %>%
  dplyr::group_modify(~{
    .safe_binary(.x)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(host_genus) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(host_genus, codon) %>%
  dplyr::arrange(host_genus, codon)

readr::write_csv(bin_within,
                 file.path(out_dir, "percodon_tRNA_WITHIN_GENUS_binary_hasTRNA_vs_none.csv"))

# Binary across genera - with vs. without tRNAs
bin_across <- phage_delta_filt %>%
  dplyr::filter(sample_type == "phage", !is.na(trna_group)) %>%
  dplyr::group_by(codon) %>%
  dplyr::group_modify(~{
    .safe_binary(.x)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::relocate(codon) %>%
  dplyr::arrange(codon)

readr::write_csv(bin_across,
                 file.path(out_dir, "percodon_tRNA_ACROSS_GENERA_binary_hasTRNA_vs_none.csv"))

# ===============================
# ∆RSCU Summary tables per genus
# ===============================

codon_levels <- filtered_codons

# ΔRSCU by genus x codon x tRNA group
genus_delta_trna_long <- phage_delta_filt %>%
  dplyr::filter(sample_type == "phage", !is.na(trna_group)) %>%
  dplyr::mutate(trna_group = forcats::fct_relevel(trna_group, trna_order)) %>%
  dplyr::group_by(host_genus, codon, trna_group) %>%
  dplyr::summarise(
    n_phages     = sum(is.finite(Delta_RSCU)),
    mean_delta   = mean(Delta_RSCU, na.rm = TRUE),
    median_delta = stats::median(Delta_RSCU, na.rm = TRUE),
    mad_delta    = stats::mad(Delta_RSCU, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(codon = factor(codon, levels = codon_levels)) %>%
  dplyr::arrange(host_genus, codon, trna_group)

genus_delta_trna_mean_wide <- genus_delta_trna_long %>%
  dplyr::select(host_genus, codon, trna_group, mean_delta) %>%
  tidyr::pivot_wider(names_from = trna_group, values_from = mean_delta)

genus_delta_trna_median_wide <- genus_delta_trna_long %>%
  dplyr::select(host_genus, codon, trna_group, median_delta) %>%
  tidyr::pivot_wider(names_from = trna_group, values_from = median_delta)

readr::write_csv(genus_delta_trna_long,
                 file.path(out_dir, "deltaRSCU_BY_GENUS_trna_LONG.csv"))
readr::write_csv(genus_delta_trna_mean_wide,
                 file.path(out_dir, "deltaRSCU_BY_GENUS_trna_MEAN_WIDE.csv"))
readr::write_csv(genus_delta_trna_median_wide,
                 file.path(out_dir, "deltaRSCU_BY_GENUS_trna_MEDIAN_WIDE.csv"))

# ΔRSCU by genus x codon x lifestyle
genus_delta_life_long <- phage_delta_filt %>%
  dplyr::filter(sample_type == "phage", !is.na(lifestyle)) %>%
  dplyr::mutate(lifestyle = forcats::fct_relevel(lifestyle, lifestyle_order)) %>%
  dplyr::group_by(host_genus, codon, lifestyle) %>%
  dplyr::summarise(
    n_phages     = sum(is.finite(Delta_RSCU)),
    mean_delta   = mean(Delta_RSCU, na.rm = TRUE),
    median_delta = stats::median(Delta_RSCU, na.rm = TRUE),
    mad_delta    = stats::mad(Delta_RSCU, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(codon = factor(codon, levels = codon_levels)) %>%
  dplyr::arrange(host_genus, codon, lifestyle)

genus_delta_life_mean_wide <- genus_delta_life_long %>%
  dplyr::select(host_genus, codon, lifestyle, mean_delta) %>%
  tidyr::pivot_wider(names_from = lifestyle, values_from = mean_delta)

genus_delta_life_median_wide <- genus_delta_life_long %>%
  dplyr::select(host_genus, codon, lifestyle, median_delta) %>%
  tidyr::pivot_wider(names_from = lifestyle, values_from = median_delta)

readr::write_csv(genus_delta_life_long,
                 file.path(out_dir, "deltaRSCU_BY_GENUS_lifestyle_LONG.csv"))
readr::write_csv(genus_delta_life_mean_wide,
                 file.path(out_dir, "deltaRSCU_BY_GENUS_lifestyle_MEAN_WIDE.csv"))
readr::write_csv(genus_delta_life_median_wide,
                 file.path(out_dir, "deltaRSCU_BY_GENUS_lifestyle_MEDIAN_WIDE.csv"))

msg("Per-codon outputs written (ΔRSCU): RAW + within/across tests + Dunn/binary + genus summaries.")

# =========================
# PERMANOVA (ΔRSCU vectors)
# =========================

phage_delta_wide <- phage_delta_filt |>
  dplyr::select(sample_id, codon, Delta_RSCU, host_genus, host_species, lifestyle, trna_group) |>
  tidyr::pivot_wider(names_from = codon, values_from = Delta_RSCU, values_fill = 0)

phage_delta_mat <- phage_delta_wide %>%
  dplyr::select(-sample_id, -host_genus, -host_species, -lifestyle, -trna_group) %>%
  as.matrix()
rownames(phage_delta_mat) <- phage_delta_wide$sample_id

permanova_df <- phage_delta_wide %>%
  dplyr::select(sample_id, host_genus, host_species, lifestyle, trna_group)

dist_delta <- stats::dist(scale(phage_delta_mat, center = TRUE, scale = TRUE), method = "euclidean")

set.seed(123)
perm_across <- vegan::adonis2(
  dist_delta ~ lifestyle + trna_group + host_genus,
  data        = permanova_df,
  permutations = 4999,
  by          = "margin",
  strata      = permanova_df$host_genus
)
sink(file.path(out_dir, "permanova_across_marginal.txt")); print(perm_across); sink()

# Dispersion checks
bd_life <- vegan::betadisper(dist_delta, group = permanova_df$lifestyle)
bd_trna <- vegan::betadisper(dist_delta, group = permanova_df$trna_group)
bd_host <- vegan::betadisper(dist_delta, group = permanova_df$host_genus)

sink(file.path(out_dir, "permanova_dispersion.txt"))
cat("\n-- Dispersion lifestyle --\n"); print(anova(bd_life)); print(vegan::permutest(bd_life, permutations = 4999))
cat("\n-- Dispersion tRNA group --\n"); print(anova(bd_trna)); print(vegan::permutest(bd_trna, permutations = 4999))
cat("\n-- Dispersion host genus --\n"); print(anova(bd_host)); print(vegan::permutest(bd_host, permutations = 4999))
sink()

# ===============================================
# Visualizing with MAE violin plots and heatmaps
# ===============================================

# ∆RSCU MAE on tRNA group
p2 <- phage_delta_scores |>
  dplyr::filter(!is.na(trna_group)) |>
  ggplot(aes(x = trna_group, y = delta_MAE, fill = trna_group)) +
  geom_violin(trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~ host_genus, scales = "free_y", drop = TRUE) +
  scale_fill_manual(
    values = palette_trna,
    labels = c(none = "None", low = "Low", high = "High"),
    name   = "tRNA Group"
  ) +
  scale_x_discrete(labels = c(none = "None", low = "Low", high = "High")) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(
    title = "ΔRSCU (MAE) by tRNA Group within Host Genus",
    x = NULL, y = "Mean |ΔRSCU| across codons"
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )

ggsave(file.path(out_dir, "plot_deltaMAE_by_trna_within_genus.png"),
       p2, width = 12, height = 8, dpi = 300)

# ∆RSCU MAE on lifestyle
p1 <- phage_delta_scores |>
  dplyr::filter(!is.na(lifestyle)) |>
  ggplot(aes(x = lifestyle, y = delta_MAE, fill = lifestyle)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~ host_genus, scales = "free_y", drop = TRUE) +
  scale_fill_manual(
    values = palette_lifestyle,
    labels = c(virulent = "Virulent", temperate = "Temperate"),
    name   = "Lifestyle"
  ) +
  scale_x_discrete(labels = c(virulent = "Virulent", temperate = "Temperate")) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(
    title = "ΔRSCU (MAE) by Lifestyle within Host Genus",
    x = NULL, y = "Mean |ΔRSCU| across codons"
  ) +
  theme(
    legend.position = "top",
    legend.title = element_text(),
    panel.grid = element_blank()
  )

ggsave(file.path(out_dir, "plot_deltaMAE_by_lifestyle_within_genus.png"),
       p1, width = 12, height = 8, dpi = 300)

# Per-genus ∆RSCU heatmaps by tRNA group
heat_df_trna <- phage_delta_filt |>
  dplyr::filter(!is.na(.data$trna_group)) |>
  dplyr::group_by(.data$host_genus, .data$trna_group, .data$codon) |>
  dplyr::summarise(mean_delta = mean(.data$Delta_RSCU, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(codon = factor(.data$codon, levels = filtered_codons))

for(gn in unique(heat_df_trna$host_genus)){
  sub <- heat_df_trna |>
    dplyr::filter(.data$host_genus == gn) |>
    tidyr::pivot_wider(names_from = codon, values_from = mean_delta, values_fill = 0) |>
    dplyr::arrange(.data$trna_group)
  
  mat <- as.matrix(sub[, setdiff(colnames(sub), c("host_genus","trna_group"))])
  rownames(mat) <- tools::toTitleCase(as.character(sub$trna_group))
  
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE, cluster_cols = FALSE,
    main = paste0("Mean ΔRSCU by Codon — ", gn, " (tRNA Group)"),
    filename = file.path(out_dir, paste0("heatmap_mean_delta_by_trna_group_", gn, ".png")),
    width = 14, height = 4
  )
}

# Per-genus ∆RSCU heatmaps by lifestyle
heat_df_life <- phage_delta_filt |>
  dplyr::filter(!is.na(.data$lifestyle)) |>
  dplyr::group_by(.data$host_genus, .data$lifestyle, .data$codon) |>
  dplyr::summarise(mean_delta = mean(.data$Delta_RSCU, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(codon = factor(.data$codon, levels = filtered_codons))

for(gn in unique(heat_df_life$host_genus)){
  sub <- heat_df_life |>
    dplyr::filter(.data$host_genus == gn) |>
    tidyr::pivot_wider(names_from = codon, values_from = mean_delta, values_fill = 0) |>
    dplyr::arrange(.data$lifestyle)
  
  mat <- as.matrix(sub[, setdiff(colnames(sub), c("host_genus","lifestyle"))])
  rownames(mat) <- tools::toTitleCase(as.character(sub$lifestyle))
  
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE, cluster_cols = FALSE,
    main = paste0("Mean ΔRSCU by Codon — ", gn, " (Lifestyle)"),
    filename = file.path(out_dir, paste0("heatmap_mean_delta_by_lifestyle_", gn, ".png")),
    width = 14, height = 4
  )
}

# RSCU Host vs. Phages
host_means_genus <- rscu_filt |>
  dplyr::filter(sample_type=="host") |>
  dplyr::group_by(host_genus, codon) |>
  dplyr::summarise(RSCU = mean(RSCU, na.rm=TRUE), .groups="drop") |>
  dplyr::mutate(group = as.character(host_genus))

phage_means_genus <- rscu_filt |>
  dplyr::filter(sample_type=="phage") |>
  dplyr::group_by(host_genus, codon) |>
  dplyr::summarise(RSCU = mean(RSCU, na.rm=TRUE), .groups="drop") |>
  dplyr::mutate(group = paste0(as.character(host_genus), " phages"))

global_means <- dplyr::bind_rows(host_means_genus, phage_means_genus) %>%
  dplyr::select(codon, group, RSCU) %>%
  dplyr::mutate(codon = factor(as.character(codon), levels = filtered_codons)) %>%
  tidyr::pivot_wider(names_from = group, values_from = RSCU) %>%
  dplyr::arrange(codon) %>%
  dplyr::mutate(codon = as.character(codon))

wanted_cols <- c(
  "Host1","Phage1",
  "Host2","Phage2",
  "...","...",
  "HostN","PhageN"
)
existing_cols <- intersect(wanted_cols, names(global_means))
global_means <- global_means %>% dplyr::select(codon, dplyr::all_of(existing_cols))

mat_global <- global_means %>%
  tibble::column_to_rownames("codon") %>%
  as.matrix()
mat_global <- suppressWarnings(apply(mat_global, 2, as.numeric))
rownames(mat_global) <- global_means$codon

mat_global_t <- t(mat_global)

pheatmap::pheatmap(
  mat_global_t,
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "Global RSCU for Hosts and Phages (by Genus)",
  fontsize_row = 12, fontsize_col = 10,
  filename = file.path(out_dir, "heatmap_global_RSCU_codonsX.png"),
  width = 14, height = 4 + 0.25 * nrow(mat_global_t)
)

# Global ∆RSCU by tRNA group
delta_trna_global <- phage_delta_filt |>
  dplyr::filter(sample_type == "phage", !is.na(trna_group), is.finite(Delta_RSCU)) |>
  dplyr::mutate(
    trna_group = factor(tolower(as.character(trna_group)), levels = trna_order),
    codon      = factor(codon, levels = filtered_codons)
  ) |>
  dplyr::group_by(codon, trna_group) |>
  dplyr::summarise(mean_delta = mean(Delta_RSCU, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = trna_group, values_from = mean_delta)

readr::write_csv(delta_trna_global, file.path(out_dir, "global_deltaRSCU_by_tRNA_pooled.csv"))

mat_trna <- delta_trna_global %>%
  tibble::column_to_rownames("codon") %>%
  as.matrix()
mat_trna <- suppressWarnings(apply(mat_trna, 2, as.numeric))
mat_trna_t <- t(mat_trna)
rownames(mat_trna_t) <- tools::toTitleCase(rownames(mat_trna_t))

pheatmap::pheatmap(
  mat_trna_t,
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "Global ΔRSCU (Phages) by tRNA Group (Pooled)",
  filename = file.path(out_dir, "heatmap_global_deltaRSCU_tRNA_codonsX.png"),
  width = 10, height = 3.5 + 0.25 * nrow(mat_trna_t)
)

# Global ∆RSCU by lifestyle
delta_life_global <- phage_delta_filt |>
  dplyr::filter(sample_type == "phage", !is.na(lifestyle), is.finite(Delta_RSCU)) |>
  dplyr::mutate(
    lifestyle = factor(tolower(as.character(lifestyle)), levels = lifestyle_order),
    codon     = factor(codon, levels = filtered_codons)
  ) |>
  dplyr::group_by(codon, lifestyle) |>
  dplyr::summarise(mean_delta = mean(Delta_RSCU, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = lifestyle, values_from = mean_delta)

readr::write_csv(delta_life_global, file.path(out_dir, "global_deltaRSCU_by_lifestyle_pooled.csv"))

mat_life <- delta_life_global %>%
  tibble::column_to_rownames("codon") %>%
  as.matrix()
mat_life <- suppressWarnings(apply(mat_life, 2, as.numeric))
mat_life_t <- t(mat_life)
rownames(mat_life_t) <- tools::toTitleCase(rownames(mat_life_t))  # <-- fixed bug

pheatmap::pheatmap(
  mat_life_t,
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "Global ΔRSCU (Phages) by Lifestyle (Pooled)",
  filename = file.path(out_dir, "heatmap_global_deltaRSCU_lifestyle_codonsX.png"),
  width = 8, height = 3.5 + 0.25 * nrow(mat_life_t)
)

msg("All done. Outputs saved in: ", out_dir)


