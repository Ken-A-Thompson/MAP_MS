# 1_read_clean_claude.R
# ──────────────────────────────────────────────────────────────────────────────
# Self-contained: loads all raw data, cleans it, then computes bench objects.
# Produces: bench_unfilt, bench_0001pct, bench_001pct, bench_best
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(vegan)
library(DescTools)
library(readxl)

source("scripts/0_functions.R")

# ── Filter Cassette ───────────────────────────────────────────────────────────
# Edit here to apply consistent taxonomy filters across all datasets
# Also potentially just edit to have general quality
taxon_filter <- function(df) {
  df %>%
    filter(phylum == "Arthropoda") %>%
    drop_na(bin_uri) %>% 
    drop_na(class)
}

# ── Load raw data ─────────────────────────────────────────────────────────────
collections_sample_data <- read_excel('data/collections_data/CBGMB01277 sample data.xlsx')

PHAUS_BOLD_Clean_NTS <- read_excel('data/BOLD_BCDM/PHAUS_parent_nts_BCDM_2.xlsx') %>%
  rename(bin_uri = bin) %>%
  # Need to filter by >85% ID to BIN
  filter(is.na(flagged) == TRUE) %>%
  taxon_filter() %>% 
  filter(Pct_ID >= 0.977)

PHAUS_BOLD_BINs <- unique(PHAUS_BOLD_Clean_NTS$bin_uri)

PHAUS_MBC_ONT_MAP <- read_tsv('data/MAP_output/PHAUS_ONT.MAP2026_03_24/2-TSV Versions of Results/Metabarcoding_Results_PHAUS_COI-5P_658_BySample.tsv') %>%
  rename(fieldid = Sample, bin_uri = BIN_Hit) %>%
  filter(fieldid %in% collections_sample_data$Lot_fieldID) %>%
  rename_with(tolower, Kingdom:Species) %>%
  rename(tot_reads = Reads, replicates = Replicates) %>%
  select(fieldid, tot_reads, ASV_ID, replicates, kingdom:bin_uri) %>%
  drop_na(bin_uri) %>%
  taxon_filter()
# 
# #
# PHAUS_MBC_ILL_MAP <- read_tsv('data/MAP_output/PHAUS_ILL.MAP2026-05-30/PHAUS_illuminaAlt_20260529.tsv') %>%
#   filter(`%ID_match_to_BIN` > 85) %>%
#     filter(Number_N <= 4) %>%
#     filter(between(Seq_Length, 409, 424)) %>%
#   rename(fieldid = Sample, bin_uri = `BIN_Hit`) %>%
#   mutate(fieldid = str_replace_all(fieldid, "-", "#")) %>%
#   filter(fieldid %in% collections_sample_data$Lot_fieldID) %>%
#   rename_with(tolower, Kingdom:Species) %>%
#   rename(tot_reads = Reads, replicates = Replicates) %>%
#   # filter(tot_reads > 1) %>%
#   select(fieldid, tot_reads, OTU_ID, replicates, kingdom:bin_uri) %>%
#   taxon_filter()

PHAUS_MBC_ILL_MAP <- read_tsv('data/MAP_output/PHAUS_ILL.MAP2026-06-02/Metabarcoding_Results_PHAUS_COI-5P_418_BySample 3.tsv') %>%
  # filter(`%ID_match_to_BIN` > 97.6) %>%
  filter(Number_N <= 4) %>%
  filter(between(Seq_Length, 409, 424)) %>%
  # filter(Reads > 1) %>%
  rename(fieldid = Sample, bin_uri = BIN_Hit, tot_reads = Reads) %>%
  mutate(fieldid = str_replace_all(fieldid, "-", "#")) %>%
  filter(fieldid %in% collections_sample_data$Lot_fieldID) %>%
  rename_with(tolower, Kingdom:Species) %>%
  rename(replicates = Replicates) %>%
  select(fieldid, tot_reads, OTU_ID, replicates, kingdom:bin_uri) %>%
  taxon_filter()

# Last good run
# PHAUS_MBC_ILL_MAP <- read_tsv('data/MAP_output/PHAUS_ILL.MAP_2026-05-22/Illumina_PHAUS_contamFILT_20260522.tsv') %>%
#   # filter(`%ID_match_to_BIN` > 97.6) %>%
#   # filter(Number_N <= 4) %>%
#   # filter(between(Seq_Length, 409, 424)) %>%
#   # filter(Reads > 1) %>%
#   # rename(fieldid = Sample, bin_uri = BIN_Hit, tot_reads = Reads) %>%
#   mutate(fieldid = str_replace_all(fieldid, "-", "#")) %>%
#   filter(fieldid %in% collections_sample_data$Lot_fieldID) %>%
#   # rename_with(tolower, Kingdom:Species) %>%
#   rename(replicates = Replicates) %>%
#   select(fieldid, tot_reads, OTU_ID, replicates, kingdom:bin_uri) %>%
#   taxon_filter()


PHAUS_MBC_ILL_spcfy <- readxl::read_xlsx('data/benchmarking/spcfy/PHAUS_Illumina_Benchmarking_SPARK.xlsx', skip = 1) %>%
  mutate(OTU_NUMBER = sub(";.*", "", OTU_ID)) 

PHAUS_SPCFY_BOLDistilled <- read_tsv('data/BOLDistilled_ID/BOLDistilled_ID_Spcfy_Jan2025/BIN_MATCH_20260415_095915.tsv') %>%
  select(Query, BIN, kingdom:species) %>%
  rename(OTU_NUMBER = Query, bin_uri = BIN) 

PHAUS_MBC_ILL_SPCFY_OTU_BOLDTax <- PHAUS_MBC_ILL_spcfy %>%
  select(OTU_ID, OTU_NUMBER) %>%
  left_join(PHAUS_SPCFY_BOLDistilled, by = "OTU_NUMBER")

PHAUS_MBC_ILL_SPCFY <- PHAUS_MBC_ILL_spcfy %>%
  select(BOLD_Process_ID:OTU_ID, GMP_58219_Rep1:GMP_58228_Rep8) %>%
  select(-`BOLD_Grade%ID`) %>%
  pivot_longer(GMP_58219_Rep1:GMP_58228_Rep8, names_to = "sample", values_to = "reads") %>%
  filter(reads != 0) %>%
  separate(sample, into = c("fieldid", "rep"), sep = "_Rep", remove = TRUE) %>%
  group_by(fieldid, OTU_ID) %>%
  summarise(replicates = n(), tot_reads = sum(reads), .groups = "drop") %>%
  left_join(PHAUS_MBC_ILL_SPCFY_OTU_BOLDTax, by = "OTU_ID") %>%
  mutate(fieldid = gsub("_", "#", fieldid)) %>%
  select(fieldid, replicates, tot_reads, bin_uri, phylum:species) %>%
  mutate(species = gsub("_", " ", species)) %>%
  drop_na(bin_uri) %>%
  taxon_filter()

PHAUS_ILL_mBRAVE <- read_tsv('data/benchmarking/mBRAVE/All_Sets_of_Data_Illumina.tsv') %>%
  select(sampleId, bin_uri, sequences, phylum:bin_uri) %>%
  separate(sampleId, into = c('prefix', 'fieldid', 'rep'), sep = "_", extra = "drop") %>%
  mutate(fieldid = gsub("-", "#", fieldid)) %>%
  filter(fieldid %in% collections_sample_data$Lot_fieldID) %>%
  select(-prefix) %>%
  group_by(fieldid, phylum, class, order, family, genus, species, bin_uri) %>%
  summarise(tot_reads = sum(sequences), replicates = length(unique(rep)), .groups = "drop") %>%
  mutate(method = "mBRAVE_ILL") %>%
  taxon_filter()

GT           <- PHAUS_BOLD_Clean_NTS
GT_parent    <- PHAUS_BOLD_Clean_NTS %>% filter(record_type == "parent")
FACTOR_0001PCT <- 1e-06   # 0.0001% of sample reads
FACTOR_001PCT  <- 1e-05   # 0.001%  of sample reads

# Apply a proportional per-sample filter (replicates logic in compute_MinFilt_Metrics)
apply_prop_filter <- function(data, factor) {
  data %>%
    group_by(fieldid) %>%
    mutate(sample_reads = sum(tot_reads),
           threshold    = factor * sample_reads) %>%
    filter(tot_reads >= threshold) %>%
    ungroup() %>%
    select(-sample_reads, -threshold)
}

# ── Dataset registry ──────────────────────────────────────────────────────────
datasets <- list(
  list(data = PHAUS_MBC_ONT_MAP,    Software = "MAP",    Dataset = "ONT", method = "MAP_ONT"),
  list(data = PHAUS_MBC_ILL_MAP,    Software = "MAP",    Dataset = "ILL", method = "MAP_ILL"),
  list(data = PHAUS_MBC_ILL_SPCFY,  Software = "SPCFY",  Dataset = "ILL", method = "SPCFY"),
  list(data = PHAUS_ILL_mBRAVE,     Software = "mBRAVE", Dataset = "ILL", method = "mBRAVE")
)

# ── Helper: assemble bench table ─────────────────────────────────────────────
assemble_bench <- function(metrics_df) {
  metrics_df %>%
    mutate(replicates = as.integer(as.character(replicates)),
           across(where(is.double), ~ round(.x, 3))) %>%
    select(Software, Dataset, tot_reads, replicates, n_reads, est, spearman, mntl, bc_mean, bc_sd, recall, f1)
}

# ── Core: compute all bench tables for a given ground truth ───────────────────
compute_all_bench <- function(datasets, GT_data, label = "") {

  pfx <- if (nchar(label)) paste0("[", label, "] ") else ""

  # Version 1: Unfiltered
  message(pfx, "Computing metrics: unfiltered...")
  unfilt_metrics <- map_dfr(datasets, function(d) {
    compute_MinFilt_Metrics(d$data, ground_truth = GT_data, threshold_factor = 0) %>%
      mutate(Software = d$Software, Dataset = d$Dataset, method = d$method)
  })
  bench_unfilt <- assemble_bench(unfilt_metrics)

  # Version 2: 0.0001% proportional filter
  message(pfx, "Computing metrics: 0.0001% proportional filter...")
  filter_0001pct_list <- set_names(
    map(datasets, function(d) apply_prop_filter(d$data, FACTOR_0001PCT)),
    map_chr(datasets, "method")
  )
  metrics_0001pct <- map_dfr(datasets, function(d) {
    compute_MinFilt_Metrics(d$data, ground_truth = GT_data, threshold_factor = FACTOR_0001PCT) %>%
      mutate(Software = d$Software, Dataset = d$Dataset, method = d$method)
  })
  bench_0001pct <- assemble_bench(metrics_0001pct)

  # Version 3: 0.001% proportional filter
  message(pfx, "Computing metrics: 0.001% proportional filter...")
  filter_001pct_list <- set_names(
    map(datasets, function(d) apply_prop_filter(d$data, FACTOR_001PCT)),
    map_chr(datasets, "method")
  )
  metrics_001pct <- map_dfr(datasets, function(d) {
    compute_MinFilt_Metrics(d$data, ground_truth = GT_data, threshold_factor = FACTOR_001PCT) %>%
      mutate(Software = d$Software, Dataset = d$Dataset, method = d$method)
  })
  bench_001pct <- assemble_bench(metrics_001pct)

  # Version 4: Best (optimised reads × replicates surface search)
  message(pfx, "Computing metrics: best filter (surface search)...")
  best_metrics <- map_dfr(datasets, function(d) {
    abg <- compute_a_b_g(MBC_data = d$data, GT_data = GT_data)
    slice(abg, which.min(rank_sum)) %>%
      select(tot_reads, replicates, est, spearman, mntl, bc_mean, bc_sd, recall, f1) %>%
      mutate(Software = d$Software, Dataset = d$Dataset, method = d$method)
  })
  best_filtered_list <- set_names(
    map(datasets, function(d) {
      br <- filter(best_metrics, method == d$method)
      d$data %>%
        filter(tot_reads  >= br$tot_reads,
               replicates >= as.integer(as.character(br$replicates)))
    }),
    map_chr(datasets, "method")
  )
  best_n_reads <- imap_dfr(best_filtered_list,
                           ~ tibble(method = .y, n_reads = sum(.x$tot_reads)))
  best_metrics <- best_metrics %>% left_join(best_n_reads, by = "method")
  bench_best   <- assemble_bench(best_metrics)

  list(
    unfilt              = bench_unfilt,
    filter_0001pct      = bench_0001pct,
    filter_001pct       = bench_001pct,
    best                = bench_best,
    filter_0001pct_list = filter_0001pct_list,
    filter_001pct_list  = filter_001pct_list,
    best_filtered_list  = best_filtered_list
  )
}

# ── Full ground truth (all record types) ─────────────────────────────────────
gt_all              <- compute_all_bench(datasets, GT, label = "all records")
bench_unfilt        <- gt_all$unfilt
bench_0001pct       <- gt_all$filter_0001pct
bench_001pct        <- gt_all$filter_001pct
bench_best          <- gt_all$best
filter_0001pct_list <- gt_all$filter_0001pct_list
filter_001pct_list  <- gt_all$filter_001pct_list
best_filtered_list  <- gt_all$best_filtered_list

# ── Parent-only ground truth (record_type == "parent") ───────────────────────
gt_par                    <- compute_all_bench(datasets, GT_parent, label = "parent only")
bench_unfilt_parent       <- gt_par$unfilt
bench_0001pct_parent      <- gt_par$filter_0001pct
bench_001pct_parent       <- gt_par$filter_001pct
bench_best_parent         <- gt_par$best

message("Done. Objects ready: bench_unfilt/0001pct/001pct/best and parent variants")
