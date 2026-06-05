# 6_ONTdownsample_claude.R
# ──────────────────────────────────────────────────────────────────────────────
# Compare three conditions (unfiltered only), all 6 metrics:
#   • ILL MAP          — single value (full Illumina dataset)
#   • ONT 47.1M        — violin: N_ITER multinomial resamples of ONT at 47.1M reads
#   • ONT Full         — single value (full ONT dataset ~110M reads)
# Output → figs_tables/claude_figs/Fig_ONT_downsample.png
# ──────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggtext)
library(ggh4x)
library(vegan)
library(DescTools)

source("scripts/0_functions.R")
source("scripts/1_read_clean_claude.R")   # PHAUS_MBC_ONT_MAP, PHAUS_MBC_ILL_MAP, GT

dir.create("figs_tables/claude_figs", recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Settings
# ══════════════════════════════════════════════════════════════════════════════
N_ITER       <- 20
MANTEL_PERMS <- 199
TARGET_M     <- 47.1
TARGET_READS <- round(TARGET_M * 1e6)

DS_LABEL <- paste0("ONT ", TARGET_M, "M\n(resampled)")

metric_meta <- tribble(
  ~key,      ~label,
  "est",      "CCC (&alpha;-diversity)",
  "spearman", "Spearman &rho; (&alpha;-diversity)",
  "mntl",     "Mantel r (&beta;-diversity)",
  "bc_mean",  "Bray-Curtis dist (&beta;-diversity)",
  "recall",   "Recall (&gamma;-diversity)",
  "f1",       "F1 (&gamma;-diversity)"
)

# ══════════════════════════════════════════════════════════════════════════════
# Helper: run compute_MinFilt_Metrics (unfiltered) and pivot to long
# ══════════════════════════════════════════════════════════════════════════════
metrics_long <- function(data, gt, label, seed = 142) {
  compute_MinFilt_Metrics(data, gt,
                          threshold_factor    = 0,
                          mantel_permutations = MANTEL_PERMS,
                          seed                = seed) %>%
    select(all_of(metric_meta$key)) %>%
    pivot_longer(everything(), names_to = "key", values_to = "value") %>%
    mutate(condition = label)
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. ONT full (unfiltered)
# ══════════════════════════════════════════════════════════════════════════════
message("Computing ONT Full (unfiltered)...")
ont_full_long <- metrics_long(PHAUS_MBC_ONT_MAP, GT, "ONT Full")

# ══════════════════════════════════════════════════════════════════════════════
# 2. ONT downsampled to 47.1M — N_ITER resamples
# ══════════════════════════════════════════════════════════════════════════════
message("Downsampling ONT to ", TARGET_M, "M (", N_ITER, " iterations)...")
total_reads <- sum(PHAUS_MBC_ONT_MAP$tot_reads)
probs       <- PHAUS_MBC_ONT_MAP$tot_reads / total_reads

set.seed(42)
iter_seeds <- sample.int(1e6, N_ITER)

ont_ds_long <- map_dfr(seq_len(N_ITER), function(i) {
  message("  iter ", i, "/", N_ITER)
  set.seed(iter_seeds[i])
  new_reads <- as.vector(rmultinom(1, size = TARGET_READS, prob = probs))
  ds <- PHAUS_MBC_ONT_MAP %>%
    mutate(tot_reads = new_reads) %>%
    filter(tot_reads > 0)

  tryCatch(
    metrics_long(ds, GT, DS_LABEL, seed = iter_seeds[i]) %>% mutate(iter = i),
    error = function(e) {
      message("    [error] ", conditionMessage(e))
      tibble(key = metric_meta$key, value = NA_real_,
             condition = DS_LABEL, iter = i)
    }
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# 3. ILL MAP (unfiltered)
# ══════════════════════════════════════════════════════════════════════════════
message("Computing ILL MAP (unfiltered)...")
ill_long <- metrics_long(PHAUS_MBC_ILL_MAP, GT, "ILL MAP")

# ══════════════════════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════════════════════
cond_levels <- c("ONT Full", DS_LABEL, "ILL MAP")

all_data <- bind_rows(
  ill_long,
  ont_ds_long %>% select(key, value, condition),
  ont_full_long
) %>%
  left_join(metric_meta, by = "key") %>%
  mutate(
    condition = factor(condition, levels = cond_levels),
    label     = factor(label,     levels = metric_meta$label)
  )

violin_df <- all_data %>% filter(condition == DS_LABEL)
point_df  <- all_data %>% filter(condition != DS_LABEL)

# ══════════════════════════════════════════════════════════════════════════════
# Plot
# ══════════════════════════════════════════════════════════════════════════════
cond_colours <- c(
  "ILL MAP"    = "#2ca02c",
  "ONT Full"   = "#d62728"
)

p <- ggplot() +
  geom_violin(data = violin_df,
              aes(x = condition, y = value),
              fill = "steelblue", colour = "steelblue4",
              alpha = 0.35, width = 0.65, trim = TRUE) +
  geom_jitter(data = violin_df,
              aes(x = condition, y = value),
              colour = "steelblue4", width = 0.07, size = 1.6, alpha = 0.75) +
  stat_summary(data = violin_df,
               aes(x = condition, y = value),
               fun = median, geom = "point",
               shape = 23, size = 4, fill = "black", colour = "black") +
  geom_point(data = point_df,
             aes(x = condition, y = value, colour = condition, fill = condition),
             shape = 23, size = 4, stroke = 0.6) +
  scale_colour_manual(values = cond_colours, guide = "none") +
  scale_fill_manual(  values = cond_colours, guide = "none") +
  facet_wrap(~ label, ncol = 2, scales = "free_y") +
  facetted_pos_scales(y = list(
    NULL,              # CCC
    NULL,              # Spearman rho
    NULL,              # Mantel r
    scale_y_reverse(), # Bray-Curtis dist (lower = better → top)
    NULL,              # Recall
    NULL               # F1
  )) +
  labs(
    x     = NULL,
    y     = "Metric value",
    title = "ONT vs Illumina MAP — unfiltered",
    subtitle = paste0(
      "Blue violin: ONT resampled to ", TARGET_M, "M reads (n = ", N_ITER, " iterations)  |  ",
      "Green diamond: ILL MAP  |  Red diamond: ONT Full"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = "grey70"),
    strip.text       = element_markdown(face = "bold", size = 10),
    axis.text.x      = element_text(size = 9)
  )

out_path <- "figs_tables/claude_figs/Fig_ONT_downsample.png"
ggsave(out_path, p, width = 280, height = 180, units = "mm", dpi = 150)
message("Done → ", out_path)
