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
# Requires in environment: GT, datasets,
#                          filter_0001pct_list, filter_001pct_list,
#                          bench_unfilt, bench_0001pct, bench_001pct

# ── Output directory ──────────────────────────────────────────────────────────
dir.create("figs_tables/claude_figs", recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Panel builders
# ══════════════════════════════════════════════════════════════════════════════

# ── (A) α-diversity: CCC scatterplot, GT on x ─────────────────────────────────
make_alpha_panel <- function(filt_data, ground_truth, spearman_val = NULL) {

  gt_rich <- ground_truth %>%
    select(fieldid, bin_uri) %>% distinct() %>%
    count(fieldid, name = "GT")

  mbc_rich <- filt_data %>%
    select(fieldid, bin_uri) %>% distinct() %>%
    count(fieldid, name = "MBC")

  joined <- left_join(gt_rich, mbc_rich, by = "fieldid") %>%
    mutate(MBC = replace_na(MBC, 0L))

  ccc_val <- tryCatch(
    signif(DescTools::CCC(joined$GT, joined$MBC)$rho.c[[1]], 2),
    error = function(e) NA_real_
  )

  subtitle_parts <- c(paste0("CCC = ", ccc_val))
  if (!is.null(spearman_val)) subtitle_parts <- c(subtitle_parts, paste0("Spearman ρ = ", spearman_val))

  ggplot(joined, aes(x = GT, y = MBC)) +
    geom_abline(slope = 1, intercept = 0,
                colour = "red", linetype = "dashed", linewidth = 0.6) +
    geom_smooth(method = "lm", colour = "steelblue",
                se = TRUE, alpha = 0.2, linewidth = 0.8) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(
      x        = "Barcode BIN richness",
      y        = "Metabarcode BIN richness",
      title    = "α-diversity",
      subtitle = paste(subtitle_parts, collapse = "  |  ")
    ) +
    theme_bw(base_size = 14) +
    theme(aspect.ratio = 1,
          axis.text    = element_text(size = 13),
          axis.title   = element_text(size = 14),
          plot.title   = element_text(size = 14, face = "bold"),
          plot.subtitle= element_text(size = 12, colour = "grey30"),
          legend.text  = element_text(size = 13),
          legend.title = element_text(size = 13))
}

# ── (B) β-diversity: PCoA with segments connecting matched samples ─────────────
make_beta_panel <- function(filt_data, ground_truth, mantel_r = NULL, bc_mean = NULL, bc_sd = NULL) {

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

  subtitle_parts <- c()
  if (!is.null(mantel_r)) subtitle_parts <- c(subtitle_parts, paste0("Mantel R = ", mantel_r))
  if (!is.null(bc_mean))  subtitle_parts <- c(subtitle_parts, paste0("BC Dist = ", bc_mean,
                                                                      if (!is.null(bc_sd)) paste0(" ± ", bc_sd) else ""))

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
    guides(fill = "none", shape = "none") +
    coord_fixed() +
    labs(x = "MDS1", y = "MDS2", title = "β-diversity",
         subtitle = if (length(subtitle_parts)) paste(subtitle_parts, collapse = "  |  ") else NULL) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          axis.text    = element_text(size = 13),
          axis.title   = element_text(size = 14),
          plot.title   = element_text(size = 14, face = "bold"),
          plot.subtitle= element_text(size = 12, colour = "grey30"))
}

# ── (C) Venn: Euler diagram of BIN overlap (ggplot, labels outside) ───────────
make_venn_panel <- function(filt_data, ground_truth, f1_score = NULL, pct_target = NULL) {

  gt_bins   <- ground_truth %>% pull(bin_uri) %>% na.omit() %>% unique()
  meta_bins <- filt_data    %>% pull(bin_uri) %>% na.omit() %>% unique()

  n_barcode_only <- length(setdiff(gt_bins,   meta_bins))
  n_overlap      <- length(intersect(gt_bins,  meta_bins))
  n_metabc_only  <- length(setdiff(meta_bins,  gt_bins))

  fit <- euler(list(BARCODE = gt_bins, METABARCODE = meta_bins))
  el  <- fit$ellipses

  # General ellipse path (handles non-circular fits)
  ellipse_path <- function(h, k, a, b, phi, n = 300) {
    theta <- seq(0, 2 * pi, length.out = n + 1)
    data.frame(
      x = h + a * cos(theta) * cos(phi) - b * sin(theta) * sin(phi),
      y = k + a * cos(theta) * sin(phi) + b * sin(theta) * cos(phi)
    )
  }

  barcode_path <- ellipse_path(el["BARCODE", "h"], el["BARCODE", "k"],
                                el["BARCODE", "a"], el["BARCODE", "b"],
                                el["BARCODE", "phi"]) %>% mutate(grp = "BARCODE")
  metabc_path  <- ellipse_path(el["METABARCODE", "h"], el["METABARCODE", "k"],
                                el["METABARCODE", "a"], el["METABARCODE", "b"],
                                el["METABARCODE", "phi"]) %>% mutate(grp = "METABARCODE")
  circles <- bind_rows(barcode_path, metabc_path)

  # Overlap centre
  x_overlap <- (el["BARCODE", "h"] + el["METABARCODE", "h"]) / 2
  y_overlap <- (el["BARCODE", "k"] + el["METABARCODE", "k"]) / 2

  # Plot extent — keep tight so circles fill the panel
  x_rng <- range(circles$x);  y_rng <- range(circles$y)
  x_lo  <- x_rng[1] - diff(x_rng) * 0.03
  x_hi  <- x_rng[2] + diff(x_rng) * 0.03
  y_lo  <- y_rng[1] - diff(y_rng) * 0.03
  y_hi  <- y_rng[2] + diff(y_rng) * 0.22   # room for corner name + count labels

  # Vertical gap between name and count in each corner
  lbl_gap <- diff(y_rng) * 0.12

  ggplot() +
    geom_polygon(data = circles,
                 aes(x = x, y = y, group = grp, fill = grp),
                 alpha = 0.65, colour = NA) +
    scale_fill_manual(values = c(BARCODE = "skyblue", METABARCODE = "#E5A840")) +
    # Overlap count (inside)
    annotate("text", x = x_overlap, y = y_overlap,
             label = format(n_overlap, big.mark = ","),
             size = 5, colour = "grey20") +
    # BARCODE: top-left corner — name then count
    annotate("text", x = x_lo, y = y_hi,
             label = "BARCODE", size = 4.2, fontface = "bold",
             hjust = 0, vjust = 1) +
    annotate("text", x = x_lo, y = y_hi - lbl_gap,
             label = format(n_barcode_only, big.mark = ","),
             size = 4.2, hjust = 0, vjust = 1) +
    # METABARCODE: top-right corner — name then count
    annotate("text", x = x_hi, y = y_hi,
             label = "METABARCODE", size = 4.2, fontface = "bold",
             hjust = 1, vjust = 1) +
    annotate("text", x = x_hi, y = y_hi - lbl_gap,
             label = format(n_metabc_only, big.mark = ","),
             size = 4.2, hjust = 1, vjust = 1) +
    labs(
      title    = "BIN overlap",
      subtitle = {
        parts <- c()
        if (!is.null(f1_score))  parts <- c(parts, paste0("F1 = ", f1_score))
        if (!is.null(pct_target)) parts <- c(parts, paste0("% Target = ", pct_target))
        if (length(parts)) paste(parts, collapse = "  |  ") else NULL
      }
    ) +
    coord_equal(clip = "off") +
    expand_limits(x = c(x_lo, x_hi), y = c(y_lo, y_hi)) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      axis.title      = element_blank(),
      legend.position = "none",
      plot.clip       = "off",
      plot.title      = element_text(size = 14, face = "plain", margin = margin(b = 2)),
      plot.subtitle   = element_text(size = 12, colour = "grey30", margin = margin(b = 4))
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# Filtered data pool
# Re-uses already-computed lists from 1_read_clean_claude.R — nothing re-runs.
# ══════════════════════════════════════════════════════════════════════════════
filt_pool <- list(
  unfilt         = set_names(map(datasets, "data"), map_chr(datasets, "method")),
  filter_0001pct = filter_0001pct_list,
  filter_001pct  = filter_001pct_list
)

# ══════════════════════════════════════════════════════════════════════════════
# Main loop — one PNG per dataset, three rows (one per filter regime)
# ══════════════════════════════════════════════════════════════════════════════
regimes <- c("unfilt", "filter_0001pct", "filter_001pct")

row_letters <- c("A", "B", "C")

regime_headings <- c(
  unfilt         = "Unfiltered",
  filter_0001pct = "Proportional filter (≥ 0.0001% of sample reads)",
  filter_001pct  = "Proportional filter (≥ 0.001% of sample reads)"
)

# Thin ggplot used as a row heading — reliably renders inside wrap_plots
make_row_header <- function(text) {
  ggplot() +
    annotate("text", x = 0.01, y = 0.5, label = text,
             hjust = 0, vjust = 0.5, size = 5, fontface = "bold") +
    theme_void() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")
}

walk(datasets, function(d) {
  message("Building combined figure: ", d$method)

  rows <- imap(regimes, function(regime, idx) {
    filt <- filt_pool[[regime]][[d$method]]

    if (is.null(filt) || nrow(filt) < 3) {
      message("  Skipping regime ", regime, " — too few rows.")
      return(NULL)
    }

    bench_obj <- switch(regime,
      unfilt         = bench_unfilt,
      filter_0001pct = bench_0001pct,
      filter_001pct  = bench_001pct
    )
    bench_row <- bench_obj %>%
      filter(Software == d$Software, Dataset == d$Dataset)

    pct_target <- if (length(bench_row$recall)) round(bench_row$recall * 100) else NULL

    p_alpha <- make_alpha_panel(filt, GT, spearman_val = bench_row$spearman)
    p_beta  <- make_beta_panel(filt, GT, mantel_r = bench_row$mntl,
                                bc_mean = bench_row$bc_mean, bc_sd = bench_row$bc_sd)
    p_venn  <- make_venn_panel(filt, GT, f1_score = bench_row$f1, pct_target = pct_target)

    heading <- paste0("(", row_letters[idx], ")  ", regime_headings[[regime]])
    header  <- make_row_header(heading)
    panels  <- (p_alpha | p_beta | p_venn) + plot_layout(widths = c(1, 1, 1.1))

    (header / panels) + plot_layout(heights = c(0.07, 1))
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(invisible(NULL))

  combined <- wrap_plots(rows, ncol = 1)

  fname <- file.path("figs_tables/claude_figs", paste0("Fig_", d$method, ".png"))
  ggsave(fname, combined, width = 300, height = 345, units = "mm", dpi = 150)
  message("  → saved: ", fname)
})

message("Done. Combined per-dataset figures in figs_tables/claude_figs/")

# ══════════════════════════════════════════════════════════════════════════════
# Joint MDS: single facet_grid — rows = filter regime, columns = method
# ══════════════════════════════════════════════════════════════════════════════

# Column order matches Table 2
mds_method_levels <- c("MAP_ONT", "MAP_ILL", "SPCFY", "mBRAVE", "MetaWorks", "QIIME2")
mds_method_labels <- c(
  MAP_ONT   = "MAP (ONT)",
  MAP_ILL   = "MAP (ILL)",
  SPCFY     = "spcfy.io",
  mBRAVE    = "mBRAVE",
  MetaWorks = "MetaWorks",
  QIIME2    = "QIIME2"
)

mds_regime_labels <- c(
  unfilt         = "Unfiltered",
  filter_0001pct = "Prop. filter ≥ 0.0001%",
  filter_001pct  = "Prop. filter ≥ 0.001%"
)

# Extract ordination data for one regime; tag rows with regime label
extract_mds_data <- function(filt_list, ground_truth, regime_label) {

  gt_long <- ground_truth %>%
    mutate(mf = paste0("BOLD.", fieldid)) %>%
    select(mf, bin_uri) %>% distinct()

  mbc_long <- imap_dfr(filt_list, function(dat, method_label) {
    dat %>%
      mutate(mf = paste0(method_label, ".", fieldid)) %>%
      select(mf, bin_uri) %>% distinct() %>% na.omit()
  })

  comm <- make_comm_matrix(bind_rows(gt_long, mbc_long), "mf", "bin_uri")
  mds  <- cmdscale(vegan::vegdist(comm, method = "bray"), k = 2)

  sc <- as.data.frame(mds) %>%
    setNames(c("MDS1", "MDS2")) %>%
    rownames_to_column("mf") %>%
    separate(mf, into = c("method", "fieldid"), sep = "\\.", extra = "merge")

  bold_sc  <- sc %>% filter(method == "BOLD")
  mbc_sc   <- sc %>% filter(method != "BOLD")
  mbc_methods <- unique(mbc_sc$method)

  bold_faceted <- map_dfr(mbc_methods, function(m) {
    bold_sc %>% mutate(facet = m, point_type = "BOLD")
  })
  mbc_faceted <- mbc_sc %>% mutate(facet = method, point_type = "META")

  points <- bind_rows(bold_faceted, mbc_faceted) %>% mutate(regime = regime_label)

  segs <- mbc_sc %>%
    left_join(bold_sc %>% select(fieldid, MDS1_b = MDS1, MDS2_b = MDS2),
              by = "fieldid") %>%
    transmute(facet = method, fieldid,
              x = MDS1_b, y = MDS2_b,
              xend = MDS1, yend = MDS2,
              regime = regime_label) %>%
    drop_na()

  list(points = points, segs = segs)
}

message("Building joint MDS ordinations...")
mds_data <- map(regimes, function(regime) {
  message("  ", regime)
  extract_mds_data(filt_pool[[regime]], GT, mds_regime_labels[[regime]])
})

all_points <- bind_rows(map(mds_data, "points")) %>%
  mutate(
    facet  = factor(facet,  levels = mds_method_levels, labels = mds_method_labels),
    regime = factor(regime, levels = mds_regime_labels)
  )
all_segs <- bind_rows(map(mds_data, "segs")) %>%
  mutate(
    facet  = factor(facet,  levels = mds_method_levels, labels = mds_method_labels),
    regime = factor(regime, levels = mds_regime_labels)
  )

p_joint <- ggplot() +
  geom_segment(data = all_segs,
               aes(x = x, y = y, xend = xend, yend = yend, colour = fieldid),
               linewidth = 0.35, alpha = 0.45, show.legend = FALSE) +
  geom_point(data = all_points,
             aes(MDS1, MDS2, fill = fieldid, shape = point_type),
             colour = "black", size = 2.0, alpha = 0.9, stroke = 0.4) +
  scale_shape_manual(values = c(BOLD = 21, META = 22),
                     labels = c(BOLD = "BOLD", META = "META")) +
  scale_fill_viridis_d(option = "turbo") +
  guides(fill  = "none",
         shape = guide_legend(title = NULL,
                              override.aes = list(size = 3.5, fill = "grey60"))) +
  facet_grid(regime ~ facet, scales = "free") +
  labs(x = "MDS1", y = "MDS2") +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", colour = "grey70"),
    strip.text.x     = element_text(face = "bold", size = 13),
    strip.text.y     = element_text(face = "bold", size = 12, angle = 0),
    axis.text        = element_text(size = 13),
    axis.title       = element_text(size = 14),
    legend.text      = element_text(size = 13),
    legend.title     = element_text(size = 13)
  )

fname <- file.path("figs_tables/claude_figs", "Fig_joint_MDS_combined.png")
ggsave(fname, p_joint, width = 420, height = 260, units = "mm", dpi = 150)
message("→ saved: ", fname)

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — MAP ONT, unfiltered only (single row: A / B / C)
# Output → figs_tables/Fig2.png
# ══════════════════════════════════════════════════════════════════════════════
message("Building Fig2 (MAP ONT, unfiltered)...")

ont_data  <- filt_pool[["unfilt"]][["MAP_ONT"]]
bench_row <- bench_unfilt %>% filter(Software == "MAP", Dataset == "ONT")
fig2_pct_target <- if (length(bench_row$recall)) round(bench_row$recall * 100) else NULL

fig2_alpha <- make_alpha_panel(ont_data, GT, spearman_val = bench_row$spearman) + labs(title = NULL)
fig2_beta  <- make_beta_panel( ont_data, GT, mantel_r = bench_row$mntl,
                                bc_mean = bench_row$bc_mean, bc_sd = bench_row$bc_sd) + labs(title = NULL)
fig2_venn  <- make_venn_panel( ont_data, GT, f1_score = bench_row$f1, pct_target = fig2_pct_target) + labs(title = NULL)

fig2 <- (fig2_alpha | fig2_beta | fig2_venn) +
  plot_layout(widths = c(1, 1, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("figs_tables/Fig2.png", fig2, width = 300, height = 108, units = "mm", dpi = 300)
message("→ saved: figs_tables/Fig2.png")
