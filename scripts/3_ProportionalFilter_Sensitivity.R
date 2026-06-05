# 3_ProportionalFilter_Sensitivity.R
# ──────────────────────────────────────────────────────────────────────────────
# Sweeps ~20 proportional filter thresholds (log-spaced around current 5e-6)
# across all four methods and plots how α/β/γ metrics respond.
#
# Requires in environment (run 1_read_clean_claude.R first):
#   PHAUS_MBC_ONT_MAP, PHAUS_MBC_ILL_MAP, PHAUS_MBC_ILL_SPCFY,
#   PHAUS_ILL_mBRAVE, PHAUS_BOLD_Clean_NTS
#
# Output: figs_tables/claude_figs/Fig_sensitivity.png
# ──────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)
library(ggtext)

source("scripts/0_functions.R")
source("scripts/1_read_clean_claude.R")

dir.create("figs_tables/claude_figs", recursive = TRUE, showWarnings = FALSE)

CURRENT_THRESHOLD <- 5e-06   # 0.0005% — the value used in the paper

# 20 log-spaced thresholds from 1e-7 to 1e-3
# (current value sits at position ~8 in this sequence)
thresholds <- 10^seq(-7, -3, length.out = 20)

# ── Method registry ───────────────────────────────────────────────────────────
methods <- list(
  list(name = "MAP (ONT)", data = PHAUS_MBC_ONT_MAP),
  list(name = "MAP (ILL)", data = PHAUS_MBC_ILL_MAP),
  list(name = "spcfy.io",  data = PHAUS_MBC_ILL_SPCFY),
  list(name = "mBRAVE",   data = PHAUS_ILL_mBRAVE)
)

# ── Sweep thresholds × methods ────────────────────────────────────────────────
# mantel_permutations = 99 keeps runtime reasonable for 20 × 4 = 80 calls;
# use 999 for final production figures.
message("Sweeping ", length(thresholds), " thresholds × ",
        length(methods), " methods (mantel perms = 99)...")

results <- map_dfr(methods, function(m) {
  message("  ", m$name, "...")
  map_dfr(thresholds, function(thr) {
    compute_MinFilt_Metrics(
      MBC_data            = m$data,
      ground_truth        = PHAUS_BOLD_Clean_NTS,
      threshold_factor    = thr,
      mantel_permutations = 99
    ) %>%
      mutate(method = m$name, threshold = thr)
  })
})

# ── Reshape for faceting ──────────────────────────────────────────────────────
# Ordered vector drives panel order; bc_mean listed last (lower = better)
metric_levels <- c(
  "CCC (&alpha;)",
  "Spearman &rho; (&alpha;)",
  "Mantel R (&beta;)",
  "Bray-Curtis dist (&beta;)",
  "% BINs recalled (&gamma;)",
  "F1-score (&gamma;)"
)

metric_labels <- c(
  est      = "CCC (&alpha;)",
  spearman = "Spearman &rho; (&alpha;)",
  mntl     = "Mantel R (&beta;)",
  bc_mean  = "Bray-Curtis dist (&beta;)",
  recall   = "% BINs recalled (&gamma;)",
  f1       = "F1-score (&gamma;)"
)

results_long <- results %>%
  select(method, threshold, est, spearman, mntl, bc_mean, recall, f1) %>%
  pivot_longer(c(est, spearman, mntl, bc_mean, recall, f1),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(
    metric        = factor(recode(metric, !!!metric_labels), levels = metric_levels),
    threshold_pct = threshold * 100
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
method_colours <- c(
  "MAP (ONT)" = "#1F4E79",
  "MAP (ILL)" = "#2ECC71",
  "spcfy.io"  = "#E67E22",
  "mBRAVE"    = "#C0392B"
)

p <- ggplot(results_long,
            aes(x = threshold_pct, y = value,
                colour = method, group = method)) +

  geom_line(linewidth = 0.85, alpha = 0.9) +
  geom_point(size = 1.8, alpha = 0.9) +

  # Current threshold marker
  geom_vline(xintercept = CURRENT_THRESHOLD * 100,
             linetype = "dashed", colour = "black", linewidth = 0.65) +

  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  facetted_pos_scales(y = list(
    NULL,              # CCC
    NULL,              # Spearman rho
    NULL,              # Mantel R
    scale_y_reverse(), # Bray-Curtis dist (lower = better → top)
    NULL,              # % BINs recalled
    NULL               # F1-score
  )) +

  scale_x_log10(
    name   = "Proportional filter threshold (% of sample reads, log₁₀ scale)",
    labels = function(x) {
      sapply(x, function(v) {
        if (is.na(v))        return("")
        if      (v >= 0.01)  sprintf("%.2f%%", v)
        else if (v >= 0.001) sprintf("%.3f%%", v)
        else                 sprintf("%.4f%%", v)
      })
    }
  ) +

  scale_colour_manual(values = method_colours, name = NULL) +

  labs(
    title    = "Metric sensitivity across proportional filter thresholds",
    subtitle = paste0("Dashed line = current threshold (5e-6 = 0.0005% of sample reads)  |  ",
                      "20 log-spaced values from 1e-7 to 1e-3"),
    y        = "Metric value"
  ) +

  theme_bw(base_size = 12) +
  theme(
    legend.position      = "bottom",
    legend.key.width     = unit(1.5, "cm"),
    panel.grid.minor     = element_blank(),
    strip.background     = element_rect(fill = "grey92", colour = "grey70"),
    strip.text           = element_markdown(face = "bold", size = 10),
    plot.title           = element_text(face = "bold", size = 13),
    plot.subtitle        = element_text(size = 10, colour = "grey30"),
    axis.text.x          = element_text(angle = 35, hjust = 1, size = 9)
  )

out_path <- "figs_tables/claude_figs/Fig_sensitivity.png"
ggsave(out_path, p, width = 340, height = 200, units = "mm", dpi = 300, bg = "white")

message("Saved: ", out_path)
