# 4_varying_prop_filter_claude.R
# ──────────────────────────────────────────────────────────────────────────────
# Sweeps 20 log-spaced proportional filter thresholds across all four methods
# and plots ALL metrics reported in the benchmark tables:
#   α  CCC                       (compute_MinFilt_Metrics → est)
#   β  Mantel R                  (compute_MinFilt_Metrics → mntl)
#   β  Mean paired Bray-Curtis   (bc_dist_per_method      → mean_Euc)
#   γ  % BINs recalled           (compute_MinFilt_Metrics → recall)
#   γ  F1-score                  (compute_MinFilt_Metrics → f1)
#
# Requires in environment (run 1_read_clean_claude.R first):
#   PHAUS_MBC_ONT_MAP, PHAUS_MBC_ILL_MAP, PHAUS_MBC_ILL_SPCFY,
#   PHAUS_ILL_mBRAVE, PHAUS_BOLD_Clean_NTS
#   + helpers: compute_MinFilt_Metrics, bc_dist_per_method, make_comm_matrix
#
# Outputs (figs_tables/):
#   Fig_sensitivity_CCC.png
#   Fig_sensitivity_MantelR.png
#   Fig_sensitivity_MantelDist.png
#   Fig_sensitivity_Recall.png
#   Fig_sensitivity_F1.png
#   Fig_sensitivity_all.png   ← combined 5-panel
# ──────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

source("scripts/0_functions.R")

CURRENT_THRESHOLD <- 5e-06   # 0.0005% — current paper value

# 20 log-spaced thresholds from 1e-7 to 1e-3
thresholds <- 10^seq(-7, -3, length.out = 20)

# ── Method registry ───────────────────────────────────────────────────────────
methods <- list(
  list(name = "MAP (ONT)", short = "MAP_ONT", data = PHAUS_MBC_ONT_MAP),
  list(name = "MAP (ILL)", short = "MAP_ILL", data = PHAUS_MBC_ILL_MAP),
  list(name = "spcfy.io",  short = "SPCFY",   data = PHAUS_MBC_ILL_SPCFY),
  list(name = "mBRAVE",   short = "mBRAVE",  data = PHAUS_ILL_mBRAVE)
)

method_colours <- c(
  "MAP (ONT)" = "#1F4E79",
  "MAP (ILL)" = "#2ECC71",
  "spcfy.io"  = "#E67E22",
  "mBRAVE"    = "#C0392B"
)

# ── Sweep: CCC / Mantel R / recall / F1 (all from compute_MinFilt_Metrics) ───
message("Sweeping ", length(thresholds), " thresholds x ", length(methods),
        " methods (Mantel perms = 999)...")

abgf_results <- map_dfr(methods, function(m) {
  message("  ", m$name, "...")
  map_dfr(thresholds, function(thr) {
    compute_MinFilt_Metrics(
      MBC_data            = m$data,
      ground_truth        = PHAUS_BOLD_Clean_NTS,
      threshold_factor    = thr,
      mantel_permutations = 999
    ) %>%
      mutate(method = m$name, threshold = thr)
  })
})

# ── Sweep: paired Bray-Curtis distance (bc_dist_per_method, one call/threshold)
message("Computing paired Bray-Curtis distances across thresholds...")

bc_results <- map_dfr(thresholds, function(thr) {
  filtered_list <- set_names(
    map(methods, function(m) {
      m$data %>%
        group_by(fieldid) %>%
        filter(tot_reads >= thr * sum(tot_reads)) %>%
        ungroup()
    }),
    map_chr(methods, "name")
  )
  bc_dist_per_method(filtered_list, PHAUS_BOLD_Clean_NTS) %>%
    mutate(threshold = thr)
})

# ── Combine all metrics ───────────────────────────────────────────────────────
all_results <- abgf_results %>%
  select(method, threshold, est, mntl, recall, f1) %>%
  left_join(bc_results %>% select(method, threshold, mean_Euc),
            by = c("method", "threshold")) %>%
  mutate(threshold_pct = threshold * 100)

message("All metrics computed.")

# ── Shared plot elements ──────────────────────────────────────────────────────
vline_thr <- CURRENT_THRESHOLD * 100

x_scale <- scale_x_log10(
  name   = "Proportional filter threshold (% of sample reads, log scale)",
  labels = function(x) {
    sapply(x, function(v) {
      if (is.na(v))        return("")
      if (v >= 0.01)       sprintf("%.2f%%", v)
      else if (v >= 0.001) sprintf("%.3f%%", v)
      else                 sprintf("%.4f%%", v)
    })
  }
)

base_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.key.width = unit(1.5, "cm"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, colour = "grey30"),
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 9)
  )

colour_scale <- scale_colour_manual(values = method_colours, name = NULL)

subtitle_txt <- "Dashed line = current threshold (5e-6 = 0.0005% of sample reads)"

# ── Plot builder ──────────────────────────────────────────────────────────────
make_plot <- function(data, y_var, y_lab, title, reverse_better = FALSE) {
  ggplot(data, aes(x = threshold_pct, y = .data[[y_var]],
                   colour = method, group = method)) +
    geom_line(linewidth = 0.85, alpha = 0.9) +
    geom_point(size = 2, alpha = 0.9) +
    geom_vline(xintercept = vline_thr,
               linetype = "dashed", colour = "black", linewidth = 0.65) +
    x_scale +
    colour_scale +
    labs(title = title, subtitle = subtitle_txt, y = y_lab) +
    base_theme
}

# ── Individual figures ────────────────────────────────────────────────────────
p_ccc <- make_plot(all_results, "est",      "CCC",
                   "alpha-diversity: CCC across proportional filter thresholds")

p_mntl <- make_plot(all_results, "mntl",    "Mantel R",
                    "beta-diversity: Mantel R across proportional filter thresholds")

p_bc   <- make_plot(all_results, "mean_Euc","Mean paired Bray-Curtis distance",
                    "beta-diversity: Bray-Curtis distance across proportional filter thresholds")

p_rec  <- make_plot(all_results, "recall",  "% BINs recalled",
                    "gamma-diversity: % BINs recalled across proportional filter thresholds")

p_f1   <- make_plot(all_results, "f1",      "F1-score",
                    "gamma-diversity: F1-score across proportional filter thresholds")

# ── Combined 5-panel figure ───────────────────────────────────────────────────
# Reuse same data in long format for faceting
panel_labels <- c(
  est     = "CCC (alpha)",
  mntl    = "Mantel R (beta)",
  mean_Euc= "Bray-Curtis dist (beta)",
  recall  = "% BINs recalled (gamma)",
  f1      = "F1-score (gamma)"
)

long_all <- all_results %>%
  pivot_longer(c(est, mntl, mean_Euc, recall, f1),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(recode(metric, !!!panel_labels),
                         levels = panel_labels))

p_all <- ggplot(long_all,
                aes(x = threshold_pct, y = value,
                    colour = method, group = method)) +
  geom_line(linewidth = 0.75, alpha = 0.9) +
  geom_point(size = 1.6, alpha = 0.9) +
  geom_vline(xintercept = vline_thr,
             linetype = "dashed", colour = "black", linewidth = 0.55) +
  x_scale +
  colour_scale +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    title    = "All benchmark metrics across proportional filter thresholds",
    subtitle = subtitle_txt,
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position      = "bottom",
    legend.key.width     = unit(1.4, "cm"),
    panel.grid.minor     = element_blank(),
    strip.background     = element_rect(fill = "#2C3E50"),
    strip.text           = element_text(colour = "white", face = "bold", size = 10),
    plot.title           = element_text(face = "bold", size = 13),
    plot.subtitle        = element_text(size = 10, colour = "grey30"),
    axis.text.x          = element_text(angle = 35, hjust = 1, size = 8)
  )

# ── Save all ──────────────────────────────────────────────────────────────────
out <- list(
  "Fig_sensitivity_CCC.png"        = list(p = p_ccc,  w = 9,  h = 6),
  "Fig_sensitivity_MantelR.png"    = list(p = p_mntl, w = 9,  h = 6),
  "Fig_sensitivity_MantelDist.png" = list(p = p_bc,   w = 9,  h = 6),
  "Fig_sensitivity_Recall.png"     = list(p = p_rec,  w = 9,  h = 6),
  "Fig_sensitivity_F1.png"         = list(p = p_f1,   w = 9,  h = 6),
  "Fig_sensitivity_all.png"        = list(p = p_all,  w = 11, h = 13)
)

for (fname in names(out)) {
  path <- file.path("figs_tables", fname)
  ggsave(path, out[[fname]]$p,
         width = out[[fname]]$w, height = out[[fname]]$h,
         dpi = 300, bg = "white")
  message("Saved: ", path)
}

message("Done. 6 figures written to figs_tables/")
