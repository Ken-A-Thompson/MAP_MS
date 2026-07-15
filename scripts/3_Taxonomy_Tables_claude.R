# 3_Taxonomy_Tables_claude.R
# ──────────────────────────────────────────────────────────────────────────────
# Genus- and family-level benchmark tables using program-provided taxonomy.
# Datasets: MAP (ILL + ONT), mBRAVE (ILL), MetaWorks (ILL), SPCFY (ILL), QIIME2 (ILL).
#
# Prerequisites: run 1_read_clean_claude.R first (objects in environment):
#   PHAUS_MBC_ONT_MAP, PHAUS_MBC_ILL_MAP, PHAUS_ILL_mBRAVE,
#   PHAUS_MBC_ILL_MetaWorks, PHAUS_MBC_ILL_SPCFY, PHAUS_MBC_ILL_QIIME2,
#   GT, FACTOR_0001PCT, FACTOR_001PCT
# QIIME2 taxonomy from: data/benchmarking/QIIME2/2026-06-30/BOLDtaxonomy_60percID.tsv
#
# Outputs (figs_tables/tax_tables/):
#   PHAUS_benchmark_combined_gen.png   — genus-level combined table
#   PHAUS_benchmark_combined_fam.png   — family-level combined table
#   Ranks_Fig_For_Paul_gen.png         — genus-level ranks (ILL + ONT separate)
#   Ranks_Fig_For_Paul_fam.png         — family-level ranks (ILL + ONT separate)
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(flextable)
library(officer)
library(vegan)
library(DescTools)

dir.create("figs_tables/tax_tables", recursive = TRUE, showWarnings = FALSE)

PUMPKIN <- "#E67E22"
NAVY    <- "#1F4E79"

# ── 1. Data helpers ───────────────────────────────────────────────────────────

CONF_MIN <- 0.6   # minimum Bayesian posterior probability for MetaWorks taxonomy

# For MetaWorks: replace the BOLDistilled genus/family with MetaWorks' own
# classifier output, masking values below CONF_MIN to NA so agg_tax drops them.
prep_mw_tax <- function(data, tax_col) {
  bp_col  <- if (tax_col == "genus") "mw_gBP"    else "mw_fBP"
  mw_col  <- if (tax_col == "genus") "mw_genus"  else "mw_family"
  data %>%
    mutate(!!tax_col := if_else(
      !is.na(.data[[bp_col]]) & .data[[bp_col]] >= CONF_MIN,
      .data[[mw_col]], NA_character_
    ))
}

# For SPCFY: replace genus/family with the consensus taxonomy columns.
prep_spcfy_tax <- function(data, tax_col) {
  src_col <- if (tax_col == "genus") "spcfy_genus" else "spcfy_family"
  data %>% mutate(!!tax_col := .data[[src_col]])
}

# For QIIME2: join BOLD taxonomy file (60% ID threshold), parse family/genus
# from the semicolon-delimited Taxon string (e.g. "k__...;f__Oecophoridae;g__...").
qiime_tax_raw <- read_tsv(
  'data/benchmarking/QIIME2/2026-06-30/BOLDtaxonomy_60percID.tsv',
  show_col_types = FALSE
) %>%
  transmute(
    OTU_ID       = `Feature ID`,
    qiime_family = str_match(Taxon, "f__([^;]+)")[, 2] %>% na_if(""),
    qiime_genus  = str_match(Taxon, "g__([^;]+)")[, 2] %>% na_if("")
  )

PHAUS_MBC_ILL_QIIME2_tax <- PHAUS_MBC_ILL_QIIME2 %>%
  left_join(qiime_tax_raw, by = "OTU_ID")

prep_qiime_tax <- function(data, tax_col) {
  src_col <- if (tax_col == "genus") "qiime_genus" else "qiime_family"
  data %>% mutate(!!tax_col := .data[[src_col]])
}

# Collapse OTU/BIN-level MBC data to a single taxonomic rank.
# Handles comma-separated multi-taxon calls (e.g. mBRAVE "Nastra, Euphyes")
# by splitting into one row per taxon before aggregating.
agg_tax <- function(data, tax_col) {
  data %>%
    filter(!is.na(.data[[tax_col]]), .data[[tax_col]] != "") %>%
    mutate(taxon = str_split(.data[[tax_col]], ",\\s*")) %>%
    unnest(taxon) %>%
    mutate(taxon = str_trim(taxon)) %>%
    filter(taxon != "") %>%
    group_by(fieldid, taxon) %>%
    summarise(tot_reads = sum(tot_reads), replicates = max(replicates), .groups = "drop")
}

# Aggregate the ground truth to a single taxonomic rank.
# A taxon inherits record_type = "parent" if any underlying BIN is a parent.
agg_gt_tax <- function(gt, tax_col) {
  gt %>%
    filter(!is.na(.data[[tax_col]]), .data[[tax_col]] != "") %>%
    group_by(fieldid, .data[[tax_col]]) %>%
    summarise(
      record_type = if_else(any(record_type == "parent"), "parent", "non-parent"),
      .groups     = "drop"
    ) %>%
    rename(taxon = all_of(tax_col))
}

# ── 2. Bench computation ──────────────────────────────────────────────────────

assemble_bench <- function(metrics_df) {
  metrics_df %>%
    mutate(replicates = as.integer(as.character(replicates)),
           across(c(est, spearman, mntl, bc_mean, bc_sd, recall, f1), ~ signif(.x, 2))) %>%
    select(Software, Dataset, tot_reads, replicates, n_reads,
           est, spearman, mntl, bc_mean, bc_sd, recall, f1)
}

compute_tax_bench <- function(tax_col) {
  message("  Aggregating data to ", tax_col, " level...")

  tax_datasets <- list(
    list(data = agg_tax(PHAUS_MBC_ONT_MAP,                                        tax_col), Software = "MAP",       Dataset = "ONT", method = "MAP_ONT"),
    list(data = agg_tax(PHAUS_MBC_ILL_MAP,                                        tax_col), Software = "MAP",       Dataset = "ILL", method = "MAP_ILL"),
    list(data = agg_tax(prep_spcfy_tax(PHAUS_MBC_ILL_SPCFY,            tax_col), tax_col), Software = "spcfy.io",  Dataset = "ILL", method = "SPCFY"),
    list(data = agg_tax(PHAUS_ILL_mBRAVE,                                         tax_col), Software = "mBRAVE",    Dataset = "ILL", method = "mBRAVE"),
    list(data = agg_tax(prep_mw_tax(PHAUS_MBC_ILL_MetaWorks,           tax_col), tax_col), Software = "MetaWorks", Dataset = "ILL", method = "MetaWorks"),
    list(data = agg_tax(prep_qiime_tax(PHAUS_MBC_ILL_QIIME2_tax,       tax_col), tax_col), Software = "QIIME2",    Dataset = "ILL", method = "QIIME2")
  )

  gt_tax <- agg_gt_tax(GT, tax_col)

  run_at_factor <- function(factor, label) {
    message("  Computing metrics: ", label, "...")
    map_dfr(tax_datasets, function(d) {
      compute_MinFilt_Metrics(
        d$data, ground_truth    = gt_tax,
        threshold_factor = factor, sp_col = "taxon"
      ) %>%
        mutate(Software = d$Software, Dataset = d$Dataset)
    }) %>% assemble_bench()
  }

  list(
    unfilt         = run_at_factor(0,             "unfiltered"),
    filter_0001pct = run_at_factor(FACTOR_0001PCT, "0.0001% filter"),
    filter_001pct  = run_at_factor(FACTOR_001PCT,  "0.001% filter")
  )
}

# ── 3. Shared display helpers ─────────────────────────────────────────────────

# ── Format a value to 2 significant figures, preserving trailing zeros ────────
fmt_sig2 <- function(x) sub("\\.$", "", sprintf("%#.2g", x))

scale_col <- function(x, reverse = FALSE) {
  x_num <- suppressWarnings(as.numeric(x))
  rng   <- range(x_num, na.rm = TRUE)
  if (diff(rng) == 0) return(rep("#FFFFFF", length(x_num)))
  norm  <- ifelse(is.na(x_num), 0, (x_num - rng[1]) / diff(rng))
  if (reverse) norm <- 1 - norm
  r <- round(255 - norm * (255 - 34))
  g <- round(255 - norm * (255 - 139))
  b <- round(255 - norm * (255 - 34))
  sprintf("#%02X%02X%02X", r, g, b)
}

style_dataset_cells <- function(ft, dataset_vec) {
  ill_rows <- which(dataset_vec == "Illumina")
  ont_rows <- which(dataset_vec == "Nanopore")
  if (length(ill_rows) > 0) {
    ft <- color(ft, i = ill_rows, j = "Dataset", color = PUMPKIN)
    ft <- bold(ft,  i = ill_rows, j = "Dataset")
  }
  if (length(ont_rows) > 0) {
    ft <- color(ft, i = ont_rows, j = "Dataset", color = NAVY)
    ft <- bold(ft,  i = ont_rows, j = "Dataset")
  }
  ft
}

# ── 4. Combined benchmark table ───────────────────────────────────────────────

make_combined_ft <- function(bench_unfilt, bench_0001pct, bench_001pct, compact = FALSE) {

  bench_all <- dplyr::bind_rows(
    bench_unfilt  %>% mutate(Filter = "Unfiltered"),
    bench_0001pct %>% mutate(Filter = "Proportional filter\n(≥ 0.0001% of sample reads)"),
    bench_001pct  %>% mutate(Filter = "Proportional filter\n(≥ 0.001% of sample reads)")
  ) %>% select(Filter, everything())

  value_cols_hi <- c("est", "spearman", "mntl", "recall", "f1")
  group_indices <- split(seq_len(nrow(bench_all)), bench_all$Filter)
  filter_runs   <- rle(as.character(bench_all$Filter))
  group_breaks  <- head(cumsum(filter_runs$lengths), -1)

  bench_display <- bench_all %>%
    mutate(
      Dataset   = recode(Dataset, "ILL" = "Illumina", "ONT" = "Nanopore"),
      tot_reads = ifelse(Filter != "Unfiltered", paste0("≥", tot_reads), as.character(tot_reads)),
      n_reads   = paste0(round(n_reads / 1e6, 1), "M"),
      bc        = paste0(fmt_sig2(bc_mean), " ± ", fmt_sig2(bc_sd)),
      est       = fmt_sig2(est),
      spearman  = fmt_sig2(spearman),
      mntl      = fmt_sig2(mntl),
      recall    = as.character(round(recall * 100)),
      f1        = fmt_sig2(f1)
    ) %>%
    select(Filter, Software, Dataset, tot_reads, replicates, n_reads,
           est, spearman, mntl, bc, recall, f1)

  ft <- flextable(bench_display) %>%
    set_header_labels(
      Filter     = "Filter\nStrategy",
      tot_reads  = "Min.\nReads",
      replicates = "Min.\nReps",
      n_reads    = "# Reads\nin OTUs",
      est        = "CCC\n(α)",
      spearman   = "Spearman ρ\n(α)",
      mntl       = "Mantel R\n(β)",
      bc         = "BC Dist ± SD\n(β)",
      recall     = "% Target\n(γ)",
      f1         = "F1-score\n(γ)"
    ) %>%
    flextable::compose(part = "header", j = "mntl",     value = as_paragraph("Mantel ", as_i("R"), "\n(β)")) %>%
    flextable::compose(part = "header", j = "spearman", value = as_paragraph("Spearman ", as_i("ρ"), "\n(α)")) %>%
    flextable::compose(part = "header", j = "est",      value = as_paragraph("CCC\n(α)")) %>%
    flextable::compose(part = "header", j = "bc",       value = as_paragraph("BC Dist ± SD\n(β)")) %>%
    flextable::compose(part = "header", j = "recall",   value = as_paragraph("% Target\n(γ)")) %>%
    flextable::compose(part = "header", j = "f1",       value = as_paragraph("F1-score\n(γ)")) %>%
    merge_v(j = "Filter") %>%
    bold(part = "header") %>%
    bg(part = "header", bg = "#2C3E50") %>%
    color(part = "header", color = "white") %>%
    bg(part = "body", bg = "white") %>%
    bold(j = "Filter", part = "body") %>%
    bg(j = "Filter", bg = "#ECF0F1", part = "body") %>%
    align(align = "center", part = "all") %>%
    align(j  = "Filter", align = "center", part = "body") %>%
    valign(j = "Filter", valign = "center", part = "body") %>%
    border_outer(part = "all", border = fp_border(color = "#666666", width = 1.5)) %>%
    border_inner_h(part = "all", border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
    border_inner_v(part = "all", border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
    vline(j = "Filter", border = fp_border(color = "#2C3E50", width = 1.2), part = "body") %>%
    vline(j = "Filter", border = fp_border(color = "#2C3E50", width = 1.2), part = "header") %>%
    hline(i = group_breaks,      border = fp_border(color = "#2C3E50", width = 1.8)) %>%
    padding(i = group_breaks,     padding.bottom = if (compact) 5 else 9, part = "body") %>%
    padding(i = group_breaks + 1, padding.top    = if (compact) 5 else 9, part = "body") %>%
    padding(padding = if (compact) 3 else 8, part = "all") %>%
    fontsize(size = if (compact) 9 else 10, part = "all") %>%
    fontsize(size = if (compact) 9 else 11, j = "Filter", part = "body") %>%
    font(fontname = "Calibri", part = "all") %>%
    line_spacing(space = if (compact) 1.0 else 1.3, part = "all") %>%
    height(height = if (compact) 0.28 else 0.4, part = "header") %>%
    autofit() %>%
    flextable::width(j = "Filter", width = if (compact) 1.2 else 1.5)

  for (col in value_cols_hi) {
    for (grp in group_indices) {
      colors <- scale_col(bench_all[[col]][grp])
      for (k in seq_along(grp)) ft <- bg(ft, i = grp[k], j = col, bg = colors[k])
    }
  }
  for (grp in group_indices) {
    bc_colors <- scale_col(bench_all[["bc_mean"]][grp], reverse = TRUE)
    for (k in seq_along(grp)) ft <- bg(ft, i = grp[k], j = "bc", bg = bc_colors[k])
  }

  ft <- style_dataset_cells(ft, bench_display$Dataset)

  if (compact) {
    ft <- ft %>%
      height(height = 0.7 / 2.54, part = "body") %>%
      hrule(rule = "exact", part = "body")
  }

  ft
}

# ── 5. Ranks figure (ILL + ONT shown separately, ranked together) ─────────────

make_ranks_ft <- function(bench_unfilt, bench_0001pct, bench_001pct) {

  bench_all <- dplyr::bind_rows(
    bench_unfilt  %>% mutate(Filter = "Unfiltered"),
    bench_0001pct %>% mutate(Filter = "Proportional filter\n(≥ 0.0001% of sample reads)"),
    bench_001pct  %>% mutate(Filter = "Proportional filter\n(≥ 0.001% of sample reads)")
  ) %>% select(Filter, everything())

  group_indices <- split(seq_len(nrow(bench_all)), bench_all$Filter)
  filter_runs   <- rle(as.character(bench_all$Filter))
  group_breaks  <- head(cumsum(filter_runs$lengths), -1)

  # Within-group metric ranks (1 = best); bc_mean reversed (lower = better)
  for (col in c("est", "spearman", "mntl", "recall", "f1"))
    for (grp in group_indices)
      bench_all[[col]][grp] <- rank(-bench_all[[col]][grp], ties.method = "min")
  for (grp in group_indices)
    bench_all[["bc_mean"]][grp] <- rank(bench_all[["bc_mean"]][grp], ties.method = "min")

  bench_all <- bench_all %>%
    mutate(
      alpha_rank    = (est + spearman) / 2,
      beta_rank     = (mntl + bc_mean) / 2,
      gamma_rank    = (recall + f1) / 2,
      overall_score = (alpha_rank + beta_rank + gamma_rank) / 3
    )

  bench_all$overall_rank <- NA_real_
  for (grp in group_indices)
    bench_all[["overall_rank"]][grp] <- rank(bench_all[["overall_score"]][grp], ties.method = "min")

  bench_display <- bench_all %>%
    mutate(
      Dataset      = recode(Dataset, "ILL" = "Illumina", "ONT" = "Nanopore"),
      alpha_rank   = round(alpha_rank, 1),
      beta_rank    = round(beta_rank,  1),
      gamma_rank   = round(gamma_rank, 1),
      overall_rank = as.integer(overall_rank)
    ) %>%
    select(Filter, Software, Dataset, alpha_rank, beta_rank, gamma_rank, overall_rank)

  ft <- flextable(bench_display) %>%
    set_header_labels(
      Filter       = "Filter\nStrategy",
      alpha_rank   = "α:\nmean rank",
      beta_rank    = "β:\nmean rank",
      gamma_rank   = "γ:\nmean rank",
      overall_rank = "Overall\nRank"
    ) %>%
    flextable::compose(part = "header", j = "alpha_rank", value = as_paragraph("α:\nmean rank")) %>%
    flextable::compose(part = "header", j = "beta_rank",  value = as_paragraph("β:\nmean rank")) %>%
    flextable::compose(part = "header", j = "gamma_rank", value = as_paragraph("γ:\nmean rank")) %>%
    merge_v(j = "Filter") %>%
    bold(part = "header") %>%
    bg(part = "header", bg = "#2C3E50") %>%
    color(part = "header", color = "white") %>%
    bg(part = "body", bg = "white") %>%
    bold(j = "Filter", part = "body") %>%
    bg(j = "Filter", bg = "#ECF0F1", part = "body") %>%
    align(align = "center", part = "all") %>%
    align(j  = "Filter", align = "center", part = "body") %>%
    valign(j = "Filter", valign = "center", part = "body") %>%
    border_outer(part = "all", border = fp_border(color = "#666666", width = 1.5)) %>%
    border_inner_h(part = "all", border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
    border_inner_v(part = "all", border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
    vline(j = "Filter", border = fp_border(color = "#2C3E50", width = 1.2), part = "body") %>%
    vline(j = "Filter", border = fp_border(color = "#2C3E50", width = 1.2), part = "header") %>%
    hline(i = group_breaks,      border = fp_border(color = "#2C3E50", width = 1.8)) %>%
    padding(i = group_breaks,     padding.bottom = 9, part = "body") %>%
    padding(i = group_breaks + 1, padding.top    = 9, part = "body") %>%
    padding(padding = 8, part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    fontsize(size = 11, j = "Filter", part = "body") %>%
    font(fontname = "Calibri", part = "all") %>%
    line_spacing(space = 1.3, part = "all") %>%
    height(height = 0.4, part = "header") %>%
    add_footer_lines(paste0(
      "α: mean rank of CCC and Spearman ρ   |   ",
      "β: mean rank of Mantel R and Bray-Curtis distance   |   ",
      "γ: mean rank of % Target recall and F1-score"
    )) %>%
    fontsize(part = "footer", size = 9) %>%
    italic(part = "footer") %>%
    color(part = "footer", color = "#555555") %>%
    autofit() %>%
    flextable::width(j = "Filter", width = 1.5)

  for (col in c("alpha_rank", "beta_rank", "gamma_rank", "overall_rank")) {
    for (grp in group_indices) {
      colors <- scale_col(bench_all[[col]][grp], reverse = TRUE)
      for (k in seq_along(grp)) ft <- bg(ft, i = grp[k], j = col, bg = colors[k])
    }
  }

  style_dataset_cells(ft, bench_display$Dataset)
}

# ── 6. Main execution ─────────────────────────────────────────────────────────

for (tax_col in c("genus", "family")) {

  suffix    <- if (tax_col == "genus") "gen" else "fam"
  tax_label <- if (tax_col == "genus") "Genus" else "Family"
  message("\n── ", toupper(tax_col), " ──────────────────────────────────────────────")

  b <- compute_tax_bench(tax_col)

  combined_path <- paste0("figs_tables/tax_tables/PHAUS_benchmark_combined_", suffix, ".png")
  save_as_image(make_combined_ft(b$unfilt, b$filter_0001pct, b$filter_001pct),
                path = combined_path, zoom = 3, expand = 10)
  message("Saved: ", combined_path)

  docx_path <- paste0("figs_tables/tax_tables/PHAUS_benchmark_combined_", suffix, ".docx")
  doc <- read_docx() %>%
    body_set_default_section(
      prop_section(
        page_size    = page_size(orient = "landscape"),
        page_margins = page_mar(top = 0.4, bottom = 0.4, left = 0.4, right = 0.4, gutter = 0)
      )
    ) %>%
    body_add_par(
      paste0("Benchmark comparison — ", tax_label, " level: all filter strategies"),
      style = "heading 1"
    ) %>%
    body_add_flextable(make_combined_ft(b$unfilt, b$filter_0001pct, b$filter_001pct, compact = TRUE))
  print(doc, target = docx_path)
  message("Saved: ", docx_path)

  ranks_path <- paste0("figs_tables/tax_tables/Ranks_Fig_For_Paul_", suffix, ".png")
  save_as_image(make_ranks_ft(b$unfilt, b$filter_0001pct, b$filter_001pct),
                path = ranks_path, zoom = 3, expand = 10)
  message("Saved: ", ranks_path)
}

message("\nDone. All taxonomy-level figures saved to figs_tables/tax_tables/")
