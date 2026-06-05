# 4_Supp_analyses_claude.R
# ──────────────────────────────────────────────────────────────────────────────
# Prerequisite: 1_read_clean_claude.R must be run first.
#   Provides: GT, SMART_FACTOR, taxon_filter, collections_sample_data,
#             euclid_from_nmds, assemble_bench (and sources 0_functions.R)
#
# Output: figs_tables/PHAUS_benchmark_MAP_97pctID.docx / .png
#   Benchmarking table for MAP only (ONT + ILL) with reads filtered to
#   ≥97.7 % ID match to BIN
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(flextable)
library(officer)

ID_THRESHOLD <- 97.7

# ── Re-load MAP datasets keeping %ID_match_to_BIN, then filter ───────────────
MAP_ONT_97 <- read_tsv(
  'data/MAP_output/PHAUS_ONT.MAP2026_03_24/2-TSV Versions of Results/Metabarcoding_Results_PHAUS_COI-5P_658_BySample.tsv'
) %>%
  rename(fieldid = Sample, bin_uri = BIN_Hit) %>%
  filter(`%ID_match_to_BIN` >= ID_THRESHOLD) %>%
  filter(fieldid %in% collections_sample_data$Lot_fieldID) %>%
  rename_with(tolower, Kingdom:Species) %>%
  rename(tot_reads = Reads, replicates = Replicates) %>%
  select(fieldid, tot_reads, ASV_ID, replicates, kingdom:bin_uri) %>%
  drop_na(bin_uri) %>%
  taxon_filter()

MAP_ILL_97 <- read_tsv(
  'data/MAP_output/PHAUS_ILL.MAP_2026-05-22/Illumina_PHAUS_contamFILT_20260522.tsv'
) %>%
  filter(`%ID_match_to_BIN` >= ID_THRESHOLD) %>%
  mutate(fieldid = str_replace_all(fieldid, "-", "#")) %>%
  filter(fieldid %in% collections_sample_data$Lot_fieldID) %>%
  rename(replicates = Replicates) %>%
  select(fieldid, tot_reads, OTU_ID, replicates, kingdom:bin_uri) %>%
  taxon_filter()

# ── Dataset registry (MAP only) ───────────────────────────────────────────────
datasets_map97 <- list(
  list(data = MAP_ONT_97, Software = "MAP", Dataset = "ONT", method = "MAP_ONT"),
  list(data = MAP_ILL_97, Software = "MAP", Dataset = "ILL", method = "MAP_ILL")
)

# ── Compute bench objects ─────────────────────────────────────────────────────
message("Computing metrics: unfiltered (≥97.7% ID)...")

unfilt_metrics_97 <- map_dfr(datasets_map97, function(d) {
  compute_MinFilt_Metrics(d$data, ground_truth = GT, threshold_factor = 0) %>%
    mutate(Software = d$Software, Dataset = d$Dataset, method = d$method)
})

unfilt_euclid_97 <- euclid_from_nmds(
  set_names(map(datasets_map97, "data"), map_chr(datasets_map97, "method")),
  GT
)

bench_unfilt_97 <- assemble_bench(unfilt_metrics_97, unfilt_euclid_97)

message("Computing metrics: best filter (≥97.7% ID)...")

best_metrics_97 <- map_dfr(datasets_map97, function(d) {
  abg <- compute_a_b_g(MBC_data = d$data, GT_data = GT)
  slice(abg, which.min(rank_sum)) %>%
    select(tot_reads, replicates, est, mntl, recall, f1) %>%
    mutate(Software = d$Software, Dataset = d$Dataset, method = d$method)
})

best_filtered_list_97 <- set_names(
  map(datasets_map97, function(d) {
    br <- filter(best_metrics_97, method == d$method)
    d$data %>%
      filter(tot_reads  >= br$tot_reads,
             replicates >= as.integer(as.character(br$replicates)))
  }),
  map_chr(datasets_map97, "method")
)

best_n_reads_97  <- imap_dfr(best_filtered_list_97,
                               ~ tibble(method = .y, n_reads = sum(.x$tot_reads)))
best_metrics_97  <- best_metrics_97 %>% left_join(best_n_reads_97, by = "method")
best_euclid_97   <- euclid_from_nmds(best_filtered_list_97, GT)
bench_best_97    <- assemble_bench(best_metrics_97, best_euclid_97)

message("Computing metrics: smart filter (≥97.7% ID)...")

smart_metrics_97 <- map_dfr(datasets_map97, function(d) {
  compute_MinFilt_Metrics(d$data, ground_truth = GT, threshold_factor = SMART_FACTOR) %>%
    mutate(Software = d$Software, Dataset = d$Dataset, method = d$method)
})

smart_filtered_list_97 <- set_names(
  map(datasets_map97, function(d) {
    d$data %>%
      group_by(fieldid) %>%
      mutate(threshold = SMART_FACTOR * sum(tot_reads)) %>%
      filter(tot_reads >= threshold) %>%
      ungroup()
  }),
  map_chr(datasets_map97, "method")
)

smart_euclid_97 <- euclid_from_nmds(smart_filtered_list_97, GT)
bench_smart_97  <- assemble_bench(smart_metrics_97, smart_euclid_97)

message("Metrics done.")

# ── Table helpers (mirrors 2_Benchmarking_Table_claude.R) ────────────────────
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

# ── Build combined table ──────────────────────────────────────────────────────
value_cols <- c("est", "mntl", "recall", "f1")

bench_all_97 <- dplyr::bind_rows(
  bench_unfilt_97 %>% mutate(Filter = "Unfiltered"),
  bench_best_97   %>% mutate(Filter = "Best filtering parameters"),
  bench_smart_97  %>% mutate(Filter = "Smart filter\n(≥ 0.001% reads in sample)")
) %>%
  select(Filter, everything())

bench_all_display_97 <- bench_all_97 %>%
  mutate(
    Software  = recode(Software, "SPCFY" = "spcfy.io"),
    Dataset   = recode(Dataset, "ILL" = "Illumina", "ONT" = "Nanopore"),
    tot_reads = ifelse(Filter == "Smart filter\n(≥ 0.001% reads in sample)",
                       paste0("≥", tot_reads),
                       as.character(tot_reads)),
    n_reads   = paste0(round(n_reads / 1e6, 1), "M")
  )

group_indices_97 <- split(seq_len(nrow(bench_all_97)), bench_all_97$Filter)
filter_runs_97   <- rle(as.character(bench_all_97$Filter))
group_breaks_97  <- head(cumsum(filter_runs_97$lengths), -1)

ft_all_97 <- flextable(bench_all_display_97 %>% select(-mean_Euc)) %>%
  set_header_labels(
    Filter     = "Filter\nStrategy",
    tot_reads  = "Min.\nReads",
    replicates = "Min.\nReps",
    n_reads    = "# Reads\nin OTUs",
    est        = "CCC\n(α)",
    mntl       = "Mantel R\n(β)",
    mean_sd    = "Mantel Distance\n(β; ± SD)",
    recall     = "% Target\n(γ)",
    f1         = "F1-score\n(γ)"
  ) %>%
  flextable::compose(
    part  = "header", j = "mntl",
    value = as_paragraph("Mantel ", as_i("R"), "\n(β)")
  ) %>%
  merge_v(j = "Filter") %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#2C3E50") %>%
  color(part = "header", color = "white") %>%
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
  hline(i = group_breaks_97,
        border = fp_border(color = "#2C3E50", width = 1.8)) %>%
  padding(i = group_breaks_97,     padding.bottom = 9, part = "body") %>%
  padding(i = group_breaks_97 + 1, padding.top    = 9, part = "body") %>%
  padding(padding = 6, part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  fontsize(size = 11, j = "Filter", part = "body") %>%
  font(fontname = "Calibri", part = "all") %>%
  height(height = 0.4, part = "header") %>%
  autofit() %>%
  width(j = "Filter", width = 1.1)

# Within-group colour scaling
for (col in value_cols) {
  for (grp in group_indices_97) {
    cell_colors <- scale_col(bench_all_97[[col]][grp])
    for (k in seq_along(grp)) {
      ft_all_97 <- bg(ft_all_97, i = grp[k], j = col, bg = cell_colors[k])
    }
  }
}

for (grp in group_indices_97) {
  euc_colors_g <- scale_col(bench_all_97$mean_Euc[grp], reverse = TRUE)
  for (k in seq_along(grp)) {
    ft_all_97 <- bg(ft_all_97, i = grp[k], j = "mean_sd", bg = euc_colors_g[k])
  }
}

ft_all_97 <- style_dataset_cells(ft_all_97, bench_all_display_97$Dataset)

# ── Save outputs ──────────────────────────────────────────────────────────────
doc_97 <- read_docx() %>%
  body_add_par(
    paste0("MAP benchmarking: ≥", ID_THRESHOLD, "% ID match to BIN"),
    style = "heading 1"
  ) %>%
  body_add_par(
    "MAP (ONT and Illumina) only. Reads filtered to ≥97.7% ID match to BIN prior to all analyses.",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(ft_all_97)

print(doc_97, target = "figs_tables/PHAUS_benchmark_MAP_97pctID.docx")
message("Saved: figs_tables/PHAUS_benchmark_MAP_97pctID.docx")

save_as_image(
  ft_all_97,
  path   = "figs_tables/PHAUS_benchmark_MAP_97pctID.png",
  zoom   = 3,
  expand = 10
)
message("Saved: figs_tables/PHAUS_benchmark_MAP_97pctID.png")
