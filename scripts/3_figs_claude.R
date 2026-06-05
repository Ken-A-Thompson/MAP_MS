# 3_figs.R
# ──────────────────────────────────────────────────────────────────────────────
# Three-panel figures: (A) α-diversity CCC, (B) β-diversity PCoA, (C) BIN Venn
# Produced for every dataset × filter-regime combination.
# Output → figs_tables/claude_figs/Fig_<method>_<regime>.png
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(patchwork)
library(eulerr)
library(vegan)
library(DescTools)
library(ggtext)

source("scripts/0_functions.R")

# Prerequisite: 1_read_clean_claude.R and 2_Benchmarking_Table_claude.R must
# be run first.
# Requires in environment: GT, datasets, best_filtered_list,
#                          filter_0001pct_list, filter_001pct_list,
#                          bench_unfilt, bench_0001pct, bench_001pct, bench_best

# ── Output directory ──────────────────────────────────────────────────────────
dir.create("figs_tables/claude_figs", recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Panel builders
# ══════════════════════════════════════════════════════════════════════════════

# ── (A) α-diversity: CCC scatterplot, GT on x ─────────────────────────────────
make_alpha_panel <- function(filt_data, ground_truth) {

  gt_rich <- ground_truth %>%
    select(fieldid, bin_uri) %>% distinct() %>%
    count(fieldid, name = "GT")

  mbc_rich <- filt_data %>%
    select(fieldid, bin_uri) %>% distinct() %>%
    count(fieldid, name = "MBC")

  joined <- left_join(gt_rich, mbc_rich, by = "fieldid") %>%
    mutate(MBC = replace_na(MBC, 0L))

  ccc_val <- tryCatch(
    round(DescTools::CCC(joined$GT, joined$MBC)$rho.c[[1]], 3),
    error = function(e) NA_real_
  )

  ggplot(joined, aes(x = GT, y = MBC)) +
    geom_abline(slope = 1, intercept = 0,
                colour = "red", linetype = "dashed", linewidth = 0.6) +
    geom_smooth(method = "lm", colour = "steelblue",
                se = TRUE, alpha = 0.2, linewidth = 0.8) +
    geom_point(size = 2.5, alpha = 0.85) +
    coord_fixed() +
    labs(
      x        = "BOLD BIN richness",
      y        = "META BIN richness",
      title    = "(A) α-diversity",
      subtitle = paste0("CCC = ", ccc_val)
    ) +
    theme_bw(base_size = 11)
}

# ── (B) β-diversity: PCoA with segments connecting matched samples ─────────────
make_beta_panel <- function(filt_data, ground_truth, mantel_r = NULL) {

  gt_long <- ground_truth %>%
    mutate(mf = paste0("BOLD.", fieldid)) %>%
    select(mf, bin_uri) %>% distinct()

  focal_long <- filt_data %>%
    mutate(mf = paste0("META.", fieldid)) %>%
    select(mf, bin_uri) %>% distinct() %>% na.omit()

  comm <- make_comm_matrix(bind_rows(gt_long, focal_long), "mf", "bin_uri")
  mds  <- cmdscale(vegan::vegdist(comm, method = "bray"), k = 2)

  sc <- as.data.frame(mds) %>%
    setNames(c("MDS1", "MDS2")) %>%
    rownames_to_column("mf") %>%
    separate(mf, into = c("method", "fieldid"), sep = "\\.", extra = "merge")

  segs <- sc %>%
    select(method, fieldid, MDS1, MDS2) %>%
    pivot_wider(names_from = method, values_from = c(MDS1, MDS2)) %>%
    transmute(fieldid,
              x    = MDS1_BOLD, y    = MDS2_BOLD,
              xend = MDS1_META, yend = MDS2_META)

  ggplot() +
    geom_segment(data = segs,
                 aes(x = x, y = y, xend = xend, yend = yend, colour = fieldid),
                 linewidth = 0.4, alpha = 0.5, show.legend = FALSE) +
    geom_point(data = sc,
               aes(MDS1, MDS2, fill = fieldid, shape = method),
               colour = "black", size = 2.5, alpha = 0.9, stroke = 0.4) +
    scale_shape_manual(values = c(BOLD = 21, META = 22),
                       labels = c(BOLD = "BOLD", META = "META")) +
    scale_fill_viridis_d(option = "turbo") +
    guides(
      fill  = "none",
      shape = guide_legend(title = NULL,
                           override.aes = list(size = 3, fill = "grey60"))
    ) +
    coord_fixed() +
    labs(x = "MDS1", y = "MDS2", title = "(B) β-diversity",
         subtitle = if (!is.null(mantel_r)) paste0("Mantel R = ", mantel_r) else NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

# ── (C) Venn: Euler diagram of BIN overlap (ggplot, labels outside) ───────────
make_venn_panel <- function(filt_data, ground_truth, f1_score = NULL) {

  gt_bins   <- ground_truth %>% pull(bin_uri) %>% na.omit() %>% unique()
  meta_bins <- filt_data    %>% pull(bin_uri) %>% na.omit() %>% unique()

  n_bold_only <- length(setdiff(gt_bins,   meta_bins))
  n_overlap   <- length(intersect(gt_bins,  meta_bins))
  n_meta_only <- length(setdiff(meta_bins,  gt_bins))

  fit <- euler(list(BOLD = gt_bins, META = meta_bins))
  el  <- fit$ellipses

  # General ellipse path (handles non-circular fits)
  ellipse_path <- function(h, k, a, b, phi, n = 300) {
    theta <- seq(0, 2 * pi, length.out = n + 1)
    data.frame(
      x = h + a * cos(theta) * cos(phi) - b * sin(theta) * sin(phi),
      y = k + a * cos(theta) * sin(phi) + b * sin(theta) * cos(phi)
    )
  }

  bold_path <- ellipse_path(el["BOLD", "h"], el["BOLD", "k"],
                             el["BOLD", "a"], el["BOLD", "b"],
                             el["BOLD", "phi"]) %>% mutate(grp = "BOLD")
  meta_path <- ellipse_path(el["META", "h"], el["META", "k"],
                             el["META", "a"], el["META", "b"],
                             el["META", "phi"]) %>% mutate(grp = "META")
  circles <- bind_rows(bold_path, meta_path)

  # Overlap centre
  x_overlap <- (el["BOLD", "h"] + el["META", "h"]) / 2
  y_overlap <- (el["BOLD", "k"] + el["META", "k"]) / 2

  # Plot extent — keep tight so circles fill the panel
  x_rng <- range(circles$x);  y_rng <- range(circles$y)
  x_lo  <- x_rng[1] - diff(x_rng) * 0.03
  x_hi  <- x_rng[2] + diff(x_rng) * 0.03
  y_lo  <- y_rng[1] - diff(y_rng) * 0.03
  y_hi  <- y_rng[2] + diff(y_rng) * 0.09   # just enough for corner labels

  # Vertical gap between name and count in each corner
  lbl_gap <- diff(y_rng) * 0.05

  ggplot() +
    geom_polygon(data = circles,
                 aes(x = x, y = y, group = grp, fill = grp),
                 alpha = 0.65, colour = NA) +
    scale_fill_manual(values = c(BOLD = "skyblue", META = "#E5A840")) +
    # Overlap count (inside)
    annotate("text", x = x_overlap, y = y_overlap,
             label = format(n_overlap, big.mark = ","),
             size = 3.8, colour = "grey20") +
    # BOLD: top-left corner — name then count
    annotate("text", x = x_lo, y = y_hi,
             label = "BOLD", size = 3.2, fontface = "bold",
             hjust = 0, vjust = 1) +
    annotate("text", x = x_lo, y = y_hi - lbl_gap,
             label = format(n_bold_only, big.mark = ","),
             size = 3.2, hjust = 0, vjust = 1) +
    # META: top-right corner — name then count
    annotate("text", x = x_hi, y = y_hi,
             label = "META", size = 3.2, fontface = "bold",
             hjust = 1, vjust = 1) +
    annotate("text", x = x_hi, y = y_hi - lbl_gap,
             label = format(n_meta_only, big.mark = ","),
             size = 3.2, hjust = 1, vjust = 1) +
    labs(
      title    = "(C) BIN overlap",
      subtitle = if (!is.null(f1_score)) paste0("F1 = ", f1_score) else NULL
    ) +
    coord_equal(clip = "off") +
    expand_limits(x = c(x_lo, x_hi), y = c(y_lo, y_hi)) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      axis.title      = element_blank(),
      legend.position = "none",
      plot.clip       = "off",
      plot.title      = element_text(face = "plain", margin = margin(b = 2)),
      plot.subtitle   = element_text(colour = "grey30", margin = margin(b = 4))
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# Filtered data pool
# Re-uses already-computed lists from 1_read_clean_claude.R — nothing re-runs.
# ══════════════════════════════════════════════════════════════════════════════
filt_pool <- list(
  unfilt         = set_names(map(datasets, "data"), map_chr(datasets, "method")),
  filter_0001pct = filter_0001pct_list,
  filter_001pct  = filter_001pct_list,
  best           = best_filtered_list
)

# ══════════════════════════════════════════════════════════════════════════════
# Main loop — one PNG per dataset × regime
# Unfiltered / 0.0001% / 0.001% = main figures
# Best = supplementary figures
# ══════════════════════════════════════════════════════════════════════════════
regimes <- c("unfilt", "filter_0001pct", "filter_001pct", "best")

regime_labels <- c(
  unfilt         = "Unfiltered",
  filter_0001pct = "0.0001% filter",
  filter_001pct  = "0.001% filter",
  best           = "Best filter (supplementary)"
)

walk(datasets, function(d) {
  walk(regimes, function(regime) {
    message("Building: ", d$method, " | ", regime)

    filt <- filt_pool[[regime]][[d$method]]

    if (nrow(filt) < 3) {
      message("  Skipping — too few rows after filtering.")
      return(invisible(NULL))
    }

    # Pull pre-computed metrics from bench objects
    bench_obj <- switch(regime,
      unfilt         = bench_unfilt,
      filter_0001pct = bench_0001pct,
      filter_001pct  = bench_001pct,
      best           = bench_best
    )
    bench_row <- bench_obj %>%
      filter(Software == d$Software, Dataset == d$Dataset)
    mntl_val <- bench_row$mntl
    f1_val   <- bench_row$f1

    p_alpha <- make_alpha_panel(filt, GT)
    p_beta  <- make_beta_panel(filt, GT, mantel_r = mntl_val)
    p_venn  <- make_venn_panel(filt, GT, f1_score = f1_val)

    regime_lbl <- regime_labels[[regime]]

    combined <- (p_alpha | p_beta | p_venn) +
      plot_layout(widths = c(1, 1, 1.1)) +
      plot_annotation(
        title = paste0("**Pipeline:** ", d$Software,
                       ";  **Data:** ", d$Dataset,
                       ";  **Filter strategy:** ", regime_lbl),
        theme = theme(plot.title = element_markdown(
          size = 13, hjust = 0))
      )

    fname <- file.path(
      "figs_tables/claude_figs",
      paste0("Fig_", d$method, "_", regime, ".png")
    )

    ggsave(fname, combined,
           width = 300, height = 110, units = "mm", dpi = 150)
    message("  → saved: ", fname)
  })
})

message("Done. Per-dataset figures in figs_tables/claude_figs/")

# ══════════════════════════════════════════════════════════════════════════════
# Joint MDS figures — all platforms in a single shared ordination, faceted
# One PNG per filtering regime (unfilt / 0.0001% / 0.001% / best)
# ══════════════════════════════════════════════════════════════════════════════

make_joint_mds_figure <- function(filt_list, ground_truth, regime_label) {

  # ── Build long data: GT + every MBC dataset ──────────────────────────────
  gt_long <- ground_truth %>%
    mutate(mf = paste0("BOLD.", fieldid)) %>%
    select(mf, bin_uri) %>% distinct()

  mbc_long <- imap_dfr(filt_list, function(dat, method_label) {
    dat %>%
      mutate(mf = paste0(method_label, ".", fieldid)) %>%
      select(mf, bin_uri) %>% distinct() %>% na.omit()
  })

  # ── Single joint ordination ───────────────────────────────────────────────
  comm <- make_comm_matrix(bind_rows(gt_long, mbc_long), "mf", "bin_uri")
  mds  <- cmdscale(vegan::vegdist(comm, method = "bray"), k = 2)

  sc <- as.data.frame(mds) %>%
    setNames(c("MDS1", "MDS2")) %>%
    rownames_to_column("mf") %>%
    separate(mf, into = c("method", "fieldid"), sep = "\\.", extra = "merge")

  bold_sc     <- sc %>% filter(method == "BOLD")
  mbc_sc      <- sc %>% filter(method != "BOLD")
  mbc_methods <- unique(mbc_sc$method)

  # ── Duplicate BOLD into every facet so it appears as reference ───────────
  bold_faceted <- map_dfr(mbc_methods, function(m) {
    bold_sc %>% mutate(facet = m, point_type = "BOLD")
  })

  mbc_faceted <- mbc_sc %>%
    mutate(facet = method, point_type = "META")

  all_points <- bind_rows(bold_faceted, mbc_faceted)

  # ── Segments: BOLD → MBC within each facet ───────────────────────────────
  segs <- mbc_sc %>%
    left_join(
      bold_sc %>% select(fieldid, MDS1_b = MDS1, MDS2_b = MDS2),
      by = "fieldid"
    ) %>%
    transmute(facet = method, fieldid,
              x = MDS1_b, y = MDS2_b,
              xend = MDS1, yend = MDS2) %>%
    drop_na()

  # ── Friendly facet labels ─────────────────────────────────────────────────
  method_labels <- c(
    MAP_ONT = "MAP (ONT)",
    MAP_ILL = "MAP (ILL)",
    SPCFY   = "SPCFY",
    mBRAVE  = "mBRAVE"
  )
  lvl_labels <- method_labels[mbc_methods]

  all_points <- all_points %>%
    mutate(facet = factor(facet, levels = mbc_methods, labels = lvl_labels))
  segs <- segs %>%
    mutate(facet = factor(facet, levels = mbc_methods, labels = lvl_labels))

  # ── Plot ──────────────────────────────────────────────────────────────────
  ggplot() +
    geom_segment(data = segs,
                 aes(x = x, y = y, xend = xend, yend = yend,
                     colour = fieldid),
                 linewidth = 0.35, alpha = 0.45, show.legend = FALSE) +
    geom_point(data = all_points,
               aes(MDS1, MDS2, fill = fieldid, shape = point_type),
               colour = "black", size = 2.2, alpha = 0.9, stroke = 0.4) +
    scale_shape_manual(
      values = c(BOLD = 21, META = 22),
      labels = c(BOLD = "BOLD", META = "META")
    ) +
    scale_fill_viridis_d(option = "turbo") +
    guides(
      fill  = "none",
      shape = guide_legend(title = NULL,
                           override.aes = list(size = 3.5, fill = "grey60"))
    ) +
    facet_wrap(~ facet, nrow = 1) +
    coord_fixed() +
    labs(
      x     = "MDS1",
      y     = "MDS2",
      title = paste0("Joint β-diversity ordination — ", regime_label)
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92", colour = "grey70"),
      strip.text       = element_text(face = "bold", size = 11)
    )
}

# ── One joint MDS figure per regime ─────────────────────────────────────────
plots <- map(regimes, function(regime) {
  message("Joint MDS: ", regime)
  make_joint_mds_figure(filt_pool[[regime]], GT, regime_labels[[regime]])
})

combined <- wrap_plots(plots, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fname <- file.path("figs_tables/claude_figs", "Fig_joint_MDS_combined.png")
ggsave(fname, combined, width = 320, height = 330, units = "mm", dpi = 150)
message("→ saved: ", fname)
