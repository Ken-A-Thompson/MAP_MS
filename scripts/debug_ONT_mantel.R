# debug_ONT_mantel.R
# ──────────────────────────────────────────────────────────────────────────────
# Bare, transparent walk-through of the Mantel R calculation for MAP_ONT
# (unfiltered) vs. the BOLD ground truth. No helper functions from
# 0_functions.R are used for the actual math — every step is inline and
# printed so the full code path is visible.
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(vegan)

source("scripts/1_read_clean_claude.R")   # gives us PHAUS_MBC_ONT_MAP and GT, cleaned + filtered

cat("\n================ INPUT DATA ================\n")
cat("PHAUS_MBC_ONT_MAP: ", nrow(PHAUS_MBC_ONT_MAP), "rows,",
    n_distinct(PHAUS_MBC_ONT_MAP$fieldid), "samples,",
    n_distinct(PHAUS_MBC_ONT_MAP$bin_uri), "unique BINs\n")
cat("GT (BOLD):         ", nrow(GT), "rows,",
    n_distinct(GT$fieldid), "samples,",
    n_distinct(GT$bin_uri), "unique BINs\n")

cat("\nGT file source: data/BOLD_BCDM/PHAUS_parent_nts_BCDM_2.xlsx\n")
cat("MAP_ONT file source: data/MAP_output/PHAUS_ONT.MAP2026_07-02/",
    "Metabarcoding_Results_PHAUS_1K_COI-5P_658_BySample.tsv\n", sep = "")

cat("\nSample tot_reads range in MAP_ONT: ",
    min(PHAUS_MBC_ONT_MAP$tot_reads), "-", max(PHAUS_MBC_ONT_MAP$tot_reads), "\n")
cat("Unfiltered means: tot_reads >= 0, replicates >= 0 (no filtering applied here)\n")

# ── Step 1: label rows by source, no filtering (this is the "Unfiltered" regime) ──
gt_long <- GT %>%
  mutate(mf = paste0("BOLD.", fieldid)) %>%
  select(mf, bin_uri) %>%
  distinct()

map_long <- PHAUS_MBC_ONT_MAP %>%
  mutate(mf = paste0("MAP.", fieldid)) %>%
  select(mf, bin_uri) %>%
  distinct() %>%
  na.omit()

cat("\n================ SAMPLE MATCHING ================\n")
gt_samples  <- unique(sub("^BOLD\\.", "", gt_long$mf))
map_samples <- unique(sub("^MAP\\.",  "", map_long$mf))
cat("GT samples (n=", length(gt_samples), "): ", paste(sort(gt_samples), collapse = ", "), "\n", sep = "")
cat("MAP samples (n=", length(map_samples), "): ", paste(sort(map_samples), collapse = ", "), "\n", sep = "")
cat("Samples present in BOTH: ", length(intersect(gt_samples, map_samples)), "\n")
cat("Samples only in GT:      ", paste(setdiff(gt_samples, map_samples), collapse = ", "), "\n")
cat("Samples only in MAP:     ", paste(setdiff(map_samples, gt_samples), collapse = ", "), "\n")

# ── Step 2: build one combined community matrix (rows = BOLD.x / MAP.x, cols = BINs) ──
all_data <- bind_rows(gt_long, map_long)

comm_all <- all_data %>%
  count(mf, bin_uri) %>%
  pivot_wider(names_from = bin_uri, values_from = n, values_fill = 0) %>%
  column_to_rownames("mf") %>%
  as.matrix()

cat("\n================ COMMUNITY MATRIX ================\n")
cat("Combined matrix dims: ", nrow(comm_all), "rows x", ncol(comm_all), "BIN columns\n")

gt_rows  <- grep("^BOLD\\.", rownames(comm_all))
map_rows <- grep("^MAP\\.",  rownames(comm_all))
cat("GT rows in matrix: ", length(gt_rows), " | MAP rows in matrix: ", length(map_rows), "\n")

BOLD_comm <- comm_all[gt_rows, , drop = FALSE]
MAP_comm  <- comm_all[map_rows, , drop = FALSE]

# Reorder MAP_comm rows to match BOLD_comm sample order (same as beta_surface/compute_MinFilt_Metrics logic)
bold_ids <- sub("^BOLD\\.", "", rownames(BOLD_comm))
map_ids  <- sub("^MAP\\.",  "", rownames(MAP_comm))
cat("\nDo BOLD and MAP rows line up 1:1 in the same order?", identical(bold_ids, map_ids), "\n")

# ── Step 3: Bray-Curtis distance matrices ──────────────────────────────────────
BOLD_dist <- vegan::vegdist(BOLD_comm, method = "bray")
MAP_dist  <- vegan::vegdist(MAP_comm,  method = "bray")

cat("\n================ MANTEL TEST ================\n")
cat("BOLD_dist: distance matrix over", attr(BOLD_dist, "Size"), "samples\n")
cat("MAP_dist:  distance matrix over", attr(MAP_dist,  "Size"), "samples\n")

set.seed(142)
mantel_result <- vegan::mantel(
  BOLD_dist,
  MAP_dist,
  method       = "pearson",
  permutations = 999
)

print(mantel_result)

cat("\n================ RESULT ================\n")
cat("Mantel R (statistic) = ", round(mantel_result$statistic, 4), "\n")
cat("p-value               = ", mantel_result$signif, "\n")
