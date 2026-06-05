library(dplyr)
library(flextable)

# ── 1. Combine BEST tibbles and join Euclidean/Mantel-distance stats ─────────
best_table <- bind_rows(
  PHAUS_MAP_ONT_BEST    %>% mutate(Software = "MAP",    Dataset = "ONT"),
  PHAUS_MAP_ILL_BEST    %>% mutate(Software = "MAP",    Dataset = "ILL"),
  PHAUS_ILL_SPCFY_BEST  %>% mutate(Software = "SPCFY",  Dataset = "ILL"),
  PHAUS_ILL_mBRAVE_BEST %>% mutate(Software = "mBRAVE", Dataset = "ILL")
) %>%
  mutate(replicates = as.integer(as.character(replicates))) %>%
  mutate(across(where(is.double), ~ round(.x, 3))) %>%
  mutate(method = case_when(Software == "MAP" ~ paste0("MAP_", Dataset), TRUE ~ Software)) %>%
  left_join(EUCLID_MEAN_SD %>% mutate(across(where(is.double), ~ round(.x, 3))),
            by = "method") %>%
  mutate(mean_sd = paste0(mean_Euc, " ± ", sd_Euc)) %>%
  # NEW column order: mean_sd between mntl and recall
  select(Software, Dataset, tot_reads, replicates, est, mntl, mean_sd, mean_Euc, recall, f1)

# ── 2. Colour scale (green → white, higher = better) ─────────────────────────
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

value_cols <- c("est", "mntl", "recall", "f1")

# ── 3. Build flextable (drop mean_Euc from display, keep for colouring) ──────
ft <- flextable(best_table %>% select(-mean_Euc)) %>%
  set_header_labels(
    tot_reads  = "Min.\nReads",
    replicates = "Min.\nReps",
    est        = "CCC\n(\u03B1)",            # α
    mntl       = "Mantel R\n(\u03B2)",       # β  (R italicised below)
    mean_sd    = "Mantel Distance\n(\u03B2; \u00B1 SD)",
    recall     = "% BINs\n(\u03B3)",         # γ
    f1         = "F1-score\n(\u03B3)"        # γ
  ) %>%
  # italicise the "R" in the Mantel R header
  compose(part = "header", j = "mntl",
          value = as_paragraph("Mantel ", as_i("R"), "\n(\u03B2)")) %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#2C3E50") %>%
  color(part = "header", color = "white") %>%
  align(align = "center", part = "all") %>%
  border_outer(part = "all", border = officer::fp_border(color = "#AAAAAA", width = 1)) %>%
  border_inner(part = "all", border = officer::fp_border(color = "#DDDDDD", width = 0.5)) %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Calibri", part = "all") %>%
  autofit()

# Apply shading to higher-is-better numeric columns
for (col in value_cols) {
  colors <- scale_col(best_table[[col]])
  for (i in seq_along(colors)) {
    ft <- bg(ft, i = i, j = col, bg = colors[i])
  }
}

# Apply shading to Mantel Distance — lower is better, so reverse = TRUE
euc_colors <- scale_col(best_table$mean_Euc, reverse = TRUE)
for (i in seq_along(euc_colors)) {
  ft <- bg(ft, i = i, j = "mean_sd", bg = euc_colors[i])
}

ft

# ── 4. Save to Word ───────────────────────────────────────────────────────────
save_as_docx(ft, path = "figs_tables/PHAUS_benchmark_table_noFilt.docx")
