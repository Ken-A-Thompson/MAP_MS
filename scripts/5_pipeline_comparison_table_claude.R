# 5_pipeline_comparison_table_claude.R
# ──────────────────────────────────────────────────────────────────────────────
# Recreates the pipeline feature-comparison table as a flextable.
# "Cloud" column renamed to "Local Compute" with checks/X's reversed.
#
# Output: figs_tables/pipeline_comparison_table.docx
#         figs_tables/pipeline_comparison_table.png
# ──────────────────────────────────────────────────────────────────────────────

library(flextable)
library(officer)
library(dplyr)

# ── Data (TRUE = green check, FALSE = red X) ──────────────────────────────────
bool_mat <- data.frame(
  Pipeline        = c("MAP", "spcfy.io", "mBRAVE", "MetaWorks", "OptimOTU", "DADA2"),
  `Long-Reads`    = c( TRUE,  FALSE,  FALSE,  FALSE,  FALSE,  FALSE),
  `Short-Reads`   = c( TRUE,   TRUE,   TRUE,   TRUE,   TRUE,   TRUE),
  `Local Compute` = c( TRUE,  FALSE,  FALSE,   TRUE,   TRUE,   TRUE),
  Demux           = c( TRUE,   TRUE,  FALSE,  FALSE,  FALSE,  FALSE),
  `BIN-match`     = c( TRUE,   TRUE,   TRUE,  FALSE,  FALSE,  FALSE),
  Free            = c( TRUE,  FALSE,   TRUE,   TRUE,   TRUE,   TRUE),
  check.names     = FALSE
)

feat_cols <- setdiff(names(bool_mat), "Pipeline")

# ── Display data frame (convert logicals → ✓ / ✗) ────────────────────────────
CHECK <- "✓"
CROSS <- "✗"

display_df <- bool_mat %>%
  mutate(across(all_of(feat_cols), ~ ifelse(.x, CHECK, CROSS)))

# ── Colours ───────────────────────────────────────────────────────────────────
GREEN_BG  <- "#C8E6C9"
RED_BG    <- "#FFCDD2"
GREEN_FG  <- "#2E7D32"
RED_FG    <- "#C62828"
HEADER_BG <- "#37474F"

# ── Build flextable ───────────────────────────────────────────────────────────
ft <- flextable(display_df) %>%
  set_header_labels(
    `Long-Reads`    = "Long-\nReads",
    `Short-Reads`   = "Short-\nReads",
    `Local Compute` = "Local\nCompute",
    `BIN-match`     = "BIN-\nmatch"
  ) %>%
  bold(part = "header") %>%
  bg(part = "header", bg = HEADER_BG) %>%
  color(part = "header", color = "white") %>%
  align(align = "center", part = "all") %>%
  align(j = "Pipeline", align = "left", part = "body") %>%
  bg(part = "body", bg = "white") %>%
  border_outer(part = "all", border = fp_border(color = "#888888", width = 1.2)) %>%
  border_inner(part = "all", border = fp_border(color = "#CCCCCC", width = 0.5)) %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Calibri", part = "all") %>%
  line_spacing(space = 1.2, part = "all") %>%
  autofit()

# Apply per-cell background + text colour
for (col in feat_cols) {
  true_rows  <- which(bool_mat[[col]])
  false_rows <- which(!bool_mat[[col]])
  if (length(true_rows)  > 0) {
    ft <- bg(ft,    i = true_rows,  j = col, bg = GREEN_BG)
    ft <- color(ft, i = true_rows,  j = col, color = GREEN_FG)
    ft <- bold(ft,  i = true_rows,  j = col)
  }
  if (length(false_rows) > 0) {
    ft <- bg(ft,    i = false_rows, j = col, bg = RED_BG)
    ft <- color(ft, i = false_rows, j = col, color = RED_FG)
    ft <- bold(ft,  i = false_rows, j = col)
  }
}

# ── Save outputs ──────────────────────────────────────────────────────────────
save_as_docx(ft, path = "figs_tables/pipeline_comparison_table.docx")
message("Saved: figs_tables/pipeline_comparison_table.docx")

save_as_image(ft, path = "figs_tables/pipeline_comparison_table.png",
              zoom = 3, expand = 10)
message("Saved: figs_tables/pipeline_comparison_table.png")
