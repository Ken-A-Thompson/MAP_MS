# run_benchmarking_tables.R
# Reproduces all PHAUS benchmarking tables and figures from pre-computed data.
#
# Usage:
#   1. Open this file in RStudio and click "Source", OR
#   2. From a terminal: Rscript run_benchmarking_tables.R
#
# Outputs (written to figs_tables/ in the same folder as this script):
#   PHAUS_benchmark_comparison.docx
#   PHAUS_benchmark_combined.docx / .png
#   PHAUS_benchmark_combined_TARGET_ONLY.docx / .png
#   PHAUS_benchmark_combined_RANKS.png
#   PHAUS_benchmark_combined_TARGET_ONLY_RANKS.png
#   Ranks_Fig_For_Paul.docx / .png

# ── Auto-install missing packages ─────────────────────────────────────────────
pkgs <- c("tidyverse", "flextable", "officer")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) {
  message("Installing missing packages: ", paste(new, collapse = ", "))
  install.packages(new, repos = "https://cloud.r-project.org")
}

library(tidyverse)
library(flextable)
library(officer)

# ── Locate script directory and set paths ─────────────────────────────────────
script_dir <- if (interactive()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))
  ][1], mustWork = FALSE) |> sub("^--file=", "", x = _))
}
if (is.na(script_dir) || script_dir == "") script_dir <- getwd()

data_file  <- file.path(script_dir, "bench_data.rda")
output_dir <- file.path(script_dir, "figs_tables")

if (!file.exists(data_file)) stop("bench_data.rda not found at: ", data_file)
dir.create(output_dir, showWarnings = FALSE)

load(data_file)
message("Loaded bench_data.rda")

# ── Dataset display colours ────────────────────────────────────────────────────
PUMPKIN <- "#E67E22"
NAVY    <- "#1F4E79"

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

# ── Individual table builder ───────────────────────────────────────────────────
make_bench_ft <- function(bench_table, prefix_ge = FALSE) {
  value_cols_hi  <- c("est", "spearman", "mntl", "recall", "f1")
  value_cols_lo  <- c("bc_mean")

  display_table <- bench_table %>%
    mutate(
      Dataset   = recode(Dataset, "ILL" = "Illumina", "ONT" = "Nanopore"),
      tot_reads = if (prefix_ge) paste0("≥", tot_reads) else as.character(tot_reads),
      n_reads   = paste0(round(n_reads / 1e6, 1), "M"),
      bc        = paste0(bc_mean, " ± ", bc_sd)
    ) %>%
    select(Software, Dataset, tot_reads, replicates, n_reads,
           est, spearman, mntl, bc, recall, f1)

  ft <- flextable(display_table) %>%
    set_header_labels(
      tot_reads  = "Min\nReads",
      replicates = "Min\nReps",
      n_reads    = "# Reads\nin OTUs",
      est        = "CCC\n(α)",
      spearman   = "Spearman ρ\n(α)",
      mntl       = "Mantel R\n(β)",
      bc         = "Bray-Curtis Dist ± SD\n(β)",
      recall     = "% Target\n(γ)",
      f1         = "F1-score\n(γ)"
    ) %>%
    flextable::compose(part = "header", j = "mntl",
      value = as_paragraph("Mantel ", as_i("R"), "\n(β)")) %>%
    flextable::compose(part = "header", j = "spearman",
      value = as_paragraph("Spearman ", as_i("ρ"), "\n(α)")) %>%
    flextable::compose(part = "header", j = "est",
      value = as_paragraph("CCC\n(α)")) %>%
    flextable::compose(part = "header", j = "bc",
      value = as_paragraph("Bray-Curtis Dist ± SD\n(β)")) %>%
    flextable::compose(part = "header", j = "recall",
      value = as_paragraph("% Target\n(γ)")) %>%
    flextable::compose(part = "header", j = "f1",
      value = as_paragraph("F1-score\n(γ)")) %>%
    bold(part = "header") %>%
    bg(part = "header", bg = "#2C3E50") %>%
    color(part = "header", color = "white") %>%
    align(align = "center", part = "all") %>%
    border_outer(part = "all", border = fp_border(color = "#AAAAAA", width = 1)) %>%
    border_inner(part = "all", border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    bg(part = "body", bg = "white") %>%
    padding(padding.top = 6, padding.bottom = 6, part = "body") %>%
    line_spacing(space = 1.3, part = "all") %>%
    autofit()

  for (col in value_cols_hi) {
    colors <- scale_col(bench_table[[col]])
    for (i in seq_along(colors)) ft <- bg(ft, i = i, j = col, bg = colors[i])
  }
  bc_colors <- scale_col(bench_table[["bc_mean"]], reverse = TRUE)
  for (i in seq_along(bc_colors)) ft <- bg(ft, i = i, j = "bc", bg = bc_colors[i])

  ft <- style_dataset_cells(ft, display_table$Dataset)
  ft
}

# ── Combined table builder ─────────────────────────────────────────────────────
make_combined_ft <- function(bench_unfilt, bench_0001pct, bench_001pct,
                             ranks = FALSE, compact = FALSE) {

  bench_all <- dplyr::bind_rows(
    bench_unfilt  %>% mutate(Filter = "Unfiltered"),
    bench_0001pct %>% mutate(Filter = "Proportional filter\n(≥ 0.0001% of sample reads)"),
    bench_001pct  %>% mutate(Filter = "Proportional filter\n(≥ 0.001% of sample reads)")
  ) %>%
    select(Filter, everything())

  value_cols_hi <- c("est", "spearman", "mntl", "recall", "f1")

  group_indices <- split(seq_len(nrow(bench_all)), bench_all$Filter)
  filter_runs   <- rle(as.character(bench_all$Filter))
  group_breaks  <- head(cumsum(filter_runs$lengths), -1)

  if (ranks) {
    for (col in value_cols_hi)
      for (grp in group_indices)
        bench_all[[col]][grp] <- rank(-bench_all[[col]][grp], ties.method = "min")
    for (grp in group_indices)
      bench_all[["bc_mean"]][grp] <- rank(bench_all[["bc_mean"]][grp], ties.method = "min")
  }

  bench_all_display <- bench_all %>%
    mutate(
      Dataset   = recode(Dataset, "ILL" = "Illumina", "ONT" = "Nanopore"),
      tot_reads = ifelse(Filter != "Unfiltered", paste0("≥", tot_reads), as.character(tot_reads)),
      n_reads   = paste0(round(n_reads / 1e6, 1), "M"),
      bc        = if (ranks) as.character(as.integer(bc_mean))
                  else paste0(bc_mean, " ± ", bc_sd)
    ) %>%
    select(Filter, Software, Dataset, tot_reads, replicates, n_reads,
           est, spearman, mntl, bc, recall, f1)

  bc_header <- if (ranks) "BC Dist\n(β)" else "BC Dist ± SD\n(β)"

  ft <- flextable(bench_all_display) %>%
    set_header_labels(
      Filter     = "Filter\nStrategy",
      tot_reads  = "Min.\nReads",
      replicates = "Min.\nReps",
      n_reads    = "# Reads\nin OTUs",
      est        = "CCC\n(α)",
      spearman   = "Spearman ρ\n(α)",
      mntl       = "Mantel R\n(β)",
      bc         = bc_header,
      recall     = "% Target\n(γ)",
      f1         = "F1-score\n(γ)"
    ) %>%
    flextable::compose(part = "header", j = "mntl",
      value = as_paragraph("Mantel ", as_i("R"), "\n(β)")) %>%
    flextable::compose(part = "header", j = "spearman",
      value = as_paragraph("Spearman ", as_i("ρ"), "\n(α)")) %>%
    flextable::compose(part = "header", j = "est",
      value = as_paragraph("CCC\n(α)")) %>%
    flextable::compose(part = "header", j = "bc",
      value = if (ranks) as_paragraph("BC Dist\n(β)") else as_paragraph("BC Dist ± SD\n(β)")) %>%
    flextable::compose(part = "header", j = "recall",
      value = as_paragraph("% Target\n(γ)")) %>%
    flextable::compose(part = "header", j = "f1",
      value = as_paragraph("F1-score\n(γ)")) %>%
    merge_v(j = "Filter") %>%
    bold(part = "header") %>%
    bg(part = "header", bg = "#2C3E50") %>%
    color(part = "header", color = "white") %>%
    bg(part = "body", bg = "white") %>%
    bold(j = "Filter", part = "body") %>%
    bg(j = "Filter", bg = "#ECF0F1", part = "body") %>%
    align(align = "center", part = "all") %>%
    align(j = "Filter", align = "center", part = "body") %>%
    valign(j = "Filter", valign = "center", part = "body") %>%
    border_outer(part = "all", border = fp_border(color = "#666666", width = 1.5)) %>%
    border_inner_h(part = "all", border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
    border_inner_v(part = "all", border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
    vline(j = "Filter", border = fp_border(color = "#2C3E50", width = 1.2), part = "body") %>%
    vline(j = "Filter", border = fp_border(color = "#2C3E50", width = 1.2), part = "header") %>%
    hline(i = group_breaks, border = fp_border(color = "#2C3E50", width = 1.8)) %>%
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
      cell_colors <- scale_col(bench_all[[col]][grp], reverse = ranks)
      for (k in seq_along(grp)) ft <- bg(ft, i = grp[k], j = col, bg = cell_colors[k])
    }
  }
  for (grp in group_indices) {
    bc_colors <- scale_col(bench_all[["bc_mean"]][grp], reverse = TRUE)
    for (k in seq_along(grp)) ft <- bg(ft, i = grp[k], j = "bc", bg = bc_colors[k])
  }

  ft <- style_dataset_cells(ft, bench_all_display$Dataset)

  if (compact) {
    ft <- ft %>%
      height(height = 0.7 / 2.54, part = "body") %>%
      hrule(rule = "exact", part = "body")
  }

  ft
}

# ── Illumina composite-rank summary ───────────────────────────────────────────
make_illumina_ranks_ft <- function(bench_unfilt, bench_0001pct, bench_001pct) {

  bench_all <- dplyr::bind_rows(
    bench_unfilt  %>% mutate(Filter = "Unfiltered"),
    bench_0001pct %>% mutate(Filter = "Proportional filter\n(≥ 0.0001% of sample reads)"),
    bench_001pct  %>% mutate(Filter = "Proportional filter\n(≥ 0.001% of sample reads)")
  ) %>%
    select(Filter, everything()) %>%
    filter(Dataset == "ILL")

  group_indices <- split(seq_len(nrow(bench_all)), bench_all$Filter)
  filter_runs   <- rle(as.character(bench_all$Filter))
  group_breaks  <- head(cumsum(filter_runs$lengths), -1)

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
      alpha_rank   = round(alpha_rank,   1),
      beta_rank    = round(beta_rank,    1),
      gamma_rank   = round(gamma_rank,   1),
      overall_rank = as.integer(overall_rank)
    ) %>%
    select(Filter, Software, alpha_rank, beta_rank, gamma_rank, overall_rank)

  ft <- flextable(bench_display) %>%
    set_header_labels(
      Filter       = "Filter\nStrategy",
      Software     = "Software",
      alpha_rank   = "a:\nmean rank",
      beta_rank    = "b:\nmean rank",
      gamma_rank   = "g:\nmean rank",
      overall_rank = "Overall\nRank"
    ) %>%
    flextable::compose(part = "header", j = "alpha_rank",
      value = as_paragraph("α:\nmean rank")) %>%
    flextable::compose(part = "header", j = "beta_rank",
      value = as_paragraph("β:\nmean rank")) %>%
    flextable::compose(part = "header", j = "gamma_rank",
      value = as_paragraph("γ:\nmean rank")) %>%
    merge_v(j = "Filter") %>%
    bold(part = "header") %>%
    bg(part = "header", bg = "#2C3E50") %>%
    color(part = "header", color = "white") %>%
    bg(part = "body", bg = "white") %>%
    bold(j = "Filter", part = "body") %>%
    bg(j = "Filter", bg = "#ECF0F1", part = "body") %>%
    align(align = "center", part = "all") %>%
    align(j = "Filter", align = "center", part = "body") %>%
    valign(j = "Filter", valign = "center", part = "body") %>%
    border_outer(part = "all", border = fp_border(color = "#666666", width = 1.5)) %>%
    border_inner_h(part = "all", border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
    border_inner_v(part = "all", border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
    vline(j = "Filter", border = fp_border(color = "#2C3E50", width = 1.2), part = "body") %>%
    vline(j = "Filter", border = fp_border(color = "#2C3E50", width = 1.2), part = "header") %>%
    hline(i = group_breaks, border = fp_border(color = "#2C3E50", width = 1.8)) %>%
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
      cell_colors <- scale_col(bench_all[[col]][grp], reverse = TRUE)
      for (k in seq_along(grp)) ft <- bg(ft, i = grp[k], j = col, bg = cell_colors[k])
    }
  }

  ft
}

# ── Build and save all outputs ─────────────────────────────────────────────────
message("Building tables...")

ft_unfilt  <- make_bench_ft(bench_unfilt)
ft_0001pct <- make_bench_ft(bench_0001pct, prefix_ge = TRUE)
ft_001pct  <- make_bench_ft(bench_001pct,  prefix_ge = TRUE)

doc <- read_docx() %>%
  body_add_par("Table 1: Unfiltered", style = "heading 1") %>%
  body_add_flextable(ft_unfilt) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Table 2: Proportional filter (≥ 0.0001% of sample reads)", style = "heading 1") %>%
  body_add_flextable(ft_0001pct) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Table 3: Proportional filter (≥ 0.001% of sample reads)", style = "heading 1") %>%
  body_add_flextable(ft_001pct)
print(doc, target = file.path(output_dir, "PHAUS_benchmark_comparison.docx"))
message("Saved: PHAUS_benchmark_comparison.docx")

ft_all         <- make_combined_ft(bench_unfilt, bench_0001pct, bench_001pct)
ft_all_compact <- make_combined_ft(bench_unfilt, bench_0001pct, bench_001pct, compact = TRUE)

doc_all <- read_docx() %>%
  body_set_default_section(prop_section(
    page_size    = page_size(orient = "landscape"),
    page_margins = page_mar(top = 0.4, bottom = 0.4, left = 0.4, right = 0.4, gutter = 0)
  )) %>%
  body_add_par("Benchmark comparison: all filter strategies", style = "heading 1") %>%
  body_add_flextable(ft_all_compact)
print(doc_all, target = file.path(output_dir, "PHAUS_benchmark_combined.docx"))
message("Saved: PHAUS_benchmark_combined.docx")

save_as_image(ft_all, path = file.path(output_dir, "PHAUS_benchmark_combined.png"), zoom = 3, expand = 10)
message("Saved: PHAUS_benchmark_combined.png")

ft_all_parent         <- make_combined_ft(bench_unfilt_parent, bench_0001pct_parent, bench_001pct_parent)
ft_all_parent_compact <- make_combined_ft(bench_unfilt_parent, bench_0001pct_parent, bench_001pct_parent, compact = TRUE)

save_as_image(ft_all_parent, path = file.path(output_dir, "PHAUS_benchmark_combined_TARGET_ONLY.png"), zoom = 3, expand = 10)
message("Saved: PHAUS_benchmark_combined_TARGET_ONLY.png")

doc_parent <- read_docx() %>%
  body_set_default_section(prop_section(
    page_size    = page_size(orient = "landscape"),
    page_margins = page_mar(top = 0.4, bottom = 0.4, left = 0.4, right = 0.4, gutter = 0)
  )) %>%
  body_add_par("Benchmark comparison: target BINs only — all filter strategies", style = "heading 1") %>%
  body_add_flextable(ft_all_parent_compact)
print(doc_parent, target = file.path(output_dir, "PHAUS_benchmark_combined_TARGET_ONLY.docx"))
message("Saved: PHAUS_benchmark_combined_TARGET_ONLY.docx")

ft_all_ranks <- make_combined_ft(bench_unfilt, bench_0001pct, bench_001pct, ranks = TRUE)
save_as_image(ft_all_ranks, path = file.path(output_dir, "PHAUS_benchmark_combined_RANKS.png"), zoom = 3, expand = 10)
message("Saved: PHAUS_benchmark_combined_RANKS.png")

ft_all_parent_ranks <- make_combined_ft(bench_unfilt_parent, bench_0001pct_parent, bench_001pct_parent, ranks = TRUE)
save_as_image(ft_all_parent_ranks, path = file.path(output_dir, "PHAUS_benchmark_combined_TARGET_ONLY_RANKS.png"), zoom = 3, expand = 10)
message("Saved: PHAUS_benchmark_combined_TARGET_ONLY_RANKS.png")

ft_paul <- make_illumina_ranks_ft(bench_unfilt, bench_0001pct, bench_001pct)
save_as_image(ft_paul, path = file.path(output_dir, "Ranks_Fig_For_Paul.png"), zoom = 3, expand = 10)
message("Saved: Ranks_Fig_For_Paul.png")

ft_paul_compact <- ft_paul %>%
  padding(padding.top = 3, padding.bottom = 3, part = "body") %>%
  padding(padding.top = 5, padding.bottom = 5, part = "header") %>%
  line_spacing(space = 1.0, part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  height(height = 0.7 / 2.54, part = "body") %>%
  hrule(rule = "exact", part = "body")

doc_paul <- read_docx() %>%
  body_set_default_section(prop_section(
    page_size    = page_size(orient = "portrait"),
    page_margins = page_mar(top = 0.4, bottom = 0.4, left = 0.4, right = 0.4, gutter = 0)
  )) %>%
  body_add_flextable(ft_paul_compact)
print(doc_paul, target = file.path(output_dir, "Ranks_Fig_For_Paul.docx"))
message("Saved: Ranks_Fig_For_Paul.docx")

message("\nDone. All outputs written to: ", output_dir)
