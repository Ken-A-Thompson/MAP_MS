# MAPLE & RAMS Validation #
# Author: Ken A. Thompson
# Script 3 of 4: Main data analysis

library(tidyverse)
library(readxl)
library(BOLDconnectR)
library(visreg)

# For revision: Remove records with flags

# Functions #

# Compute & plot Beta Diversity with grouping toggle #
plot_PA_MDS <- function(data, plot_title, group_col = c("bin_uri", "genus", "family")) { 
  
  beta_div_out_dat <- NULL
  group_col <- match.arg(group_col)  # ensures valid choice
  
  data_site_info <- data %>% 
    select(unique_sample, fieldid) %>% 
    unique() 
  
  # Beta diversity #### 
  PA_df <- data %>% 
    group_by(unique_sample, .data[[group_col]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(presence = 1) %>% 
    select(-count) %>% 
    pivot_wider(
      names_from = all_of(group_col),
      values_from = presence,
      values_fill = 0
    ) %>%
    ungroup() %>% 
    na.omit()
  
  # Site names and community matrix
  PA_site_names <- PA_df$unique_sample
  PA_mat <- PA_df %>% select(-unique_sample) %>% as.matrix()
  
  # Bray–Curtis distance (binary)
  PA_Bray_DistMat <- vegdist(PA_mat, method = "jaccard")
  
  # MDS ordination
  PA_BC_result <- cmdscale(PA_Bray_DistMat, k = 2)
  
  PA_for_plot <- as_tibble(PA_BC_result) %>% 
    mutate(unique_sample = PA_site_names) %>% 
    left_join(data_site_info, by = "unique_sample") %>% 
    separate(unique_sample, into = c("method", "fieldid"), sep = "_")
  
  # Convex hulls by collection_notes
  PA_hulls <- PA_for_plot %>%
    group_by(fieldid) %>%
    slice(chull(V1, V2))
  
  shape_values <- c(
    "Barcode" = 22,   # square with border
    "Metabarcode" = 21   # circle with border
  )
  
  out_PA_plot <- ggplot(
    PA_for_plot,
    aes(x = V1, y = V2, fill = fieldid)
  ) +
    coord_fixed() +
    scale_colour_manual(values = my_colors) +
    # hull polygons
    geom_polygon(
      data = PA_hulls,
      aes(x = V1, y = V2, fill = fieldid, colour = fieldid),
      alpha = 0.7,
      inherit.aes = T
    ) +
    
    # points with shapes for each method
    geom_point(
      aes(shape = method),
      size = 3,
      colour = "black"   # black border
    ) +
    
    # assign square to Barcode, circle to Metabarcode
    scale_shape_manual(values = c(
      "Barcode" = 22,
      "Metabarcode" = 21
    )) +
    scale_fill_manual(values = my_colors) + 
    theme_bw() +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none")
  
  beta_div_out_dat <- NULL
  
  return(out_PA_plot)
  
}

# Read raw data ####
## Sample / Plate Info from Collections Unit ####
collections_plate_info <- read_excel('data/collections_data/CBGMB01277 sample data.xlsx', sheet = 1)
collections_sample_info <- read_excel('data/collections_data/CBGMB01277 sample data.xlsx', sheet = 2, skip = 1) %>% 
  mutate(duration = collection_end_date -collection_start_date)

which_samples_shortDur <-  collections_sample_info %>% 
  filter(duration == 8) %>% 
  pull(trackingsl_sample)

which_samples_longDur <-  collections_sample_info %>% 
  filter(duration == 17) %>% 
  pull(trackingsl_sample)


## BOLD ####
PHAUS_res_BOLD_noTax <- read_tsv('data/BOLD_downloads/PHAUS_res.tsv') %>% 
  select(processid, fieldid, processid_minted_date, collection_date_start)

## MAPLE ####
mbc_results_MAPLE <- read_excel('data/metabarcode_data/Metabarcoding Results - PHAUS_CBGMB01277_COI-5P_658.xlsx', sheet = 1) %>% 
  mutate(row_number = row_number()) %>% 
  rename(bin_uri = "BIN Hit")

mbc_results_MAPLE_noTax <- mbc_results_MAPLE %>% 
  select(row_number, Sample:Replicates) %>% 
  rename(Query = row_number)

## Traits ####
trait_df <- read_csv('data/trait_data/insect_traits_secondary_sources.csv') %>% 
  select(1:6)

## Geography ####
BINs_by_country <- read_tsv('data/derived_data/BINs_by_country.tsv')


# How many ASVs?
length(unique(mbc_results_MAPLE$ASV_Name))

# Read BIN Match data #####
## Barcode - BOLD ####
BC_full <- PHAUS_res_BOLDistilled %>% 
  rename_with(tolower, kingdom:species) %>% 
  filter(!kingdom == "Bacteria") %>% 
  rename(BIN = bin_uri)

## Metabarcode - MAPLE ####
# NOTE: If BOLDistilled is update can do this without running MAPLE again
# Only consider high-confidence ASVs for this analysis
BOLDistilled_comp_min_reps <- 7

MBC_full <- mbc_results_MAPLE %>% 
  rename_with(tolower, Kingdom:Species) %>% 
  filter(!kingdom == "Bacteria") %>% 
  rename(BIN_Match = `BIN Match Status`) %>% 
  filter(Replicates >= BOLDistilled_comp_min_reps) %>% 
  group_by(bin_uri) %>% 
  slice(1) %>%       # keep one ASV per bin_uri
  ungroup()

MBC_fasta_lines <- MBC_full  %>% 
  dplyr::mutate(fasta_entry = paste0(">", Sample, "|", ASV_ID, "\n", ASV_Consensus_Sequence)) %>% 
  dplyr::pull(fasta_entry)

# write_tsv(MBC_Info_FastaJoin, 'data/derived_data/MBC_Info_FastaJoin.tsv')  

# Write to a file
writeLines(MBC_fasta_lines, "BIN_Match/PHAUS_mbc.fasta")


# Trap info / Reads info
MBC_Info_FastaJoin <- read_tsv('data/derived_data/MBC_Info_FastaJoin.tsv')

MBC_noAus <- read_tsv('BIN_Match/BIN_match_results/NoAus_MBC/PHAUS_mbc_BIN_MATCH_RESULTS.tsv')  %>% 
  separate(Query,
           into = c("Sample", "ASV_ID"),
           sep = "\\|",
           remove = FALSE) %>% 
  rename(BIN_Match = `BIN Match?`,
         bin_uri = `Hit (BIN)`)

MBC_noNSW <- read_tsv('BIN_Match/BIN_match_results/NoNSW_MBC/PHAUS_mbc_BIN_MATCH_RESULTS.tsv')  %>% 
  separate(Query,
           into = c("Sample", "ASV_ID"),
           sep = "\\|",
           remove = FALSE) %>% 
  rename(BIN_Match = `BIN Match?`,
         bin_uri = `Hit (BIN)`)

MBC_noPark <- read_tsv('BIN_Match/BIN_match_results/NoPark_MBC/PHAUS_mbc_BIN_MATCH_RESULTS.tsv')  %>% 
  separate(Query,
           into = c("Sample", "ASV_ID"),
           sep = "\\|",
           remove = FALSE) %>% 
  rename(BIN_Match = `BIN Match?`,
         bin_uri = `Hit (BIN)`)

## BOLDistilled BIN Counts ####
count_BINs <- function(data) {length(unique(filter(data, `BIN_Match` == "BIN_MATCH")$bin_uri))}
pct_tax_ID <- function(data) {
  data_unique_ASVs <- data %>% 
    group_by(ASV_ID) %>% 
    slice(1)
  
  pct_no_fam <- data_unique_ASVs$Family %>% is.na(.) %>% mean(.)
  pct_no_gen <- data_unique_ASVs$Genus %>% is.na(.) %>% mean(.)
  pct_no_spec <- data_unique_ASVs$Species %>% is.na(.) %>% mean(.)

  print(paste0(round(100 * (1-pct_no_fam), 1), "% with Family ID"))
  print(paste0(round(100 * (1-pct_no_gen), 1), "% with Genus ID"))
  print(paste0(round(100 * (1-pct_no_spec), 1), "% with Species ID"))
  
  
  
  }


n_ASVs_total <- MBC_full %>% 
  select(ASV_ID) %>% 
  unique() %>% 
  nrow(.)

count_BINs(MBC_noAus)
count_BINs(MBC_noNSW)
count_BINs(MBC_noPark)
count_BINs(MBC_full)

1 - MBC_noAus$Species %>% is.na(.) %>% mean(.)
1 - MBC_noNSW$Species %>% is.na(.) %>% mean(.)
1 - MBC_noPark$Species %>% is.na(.) %>% mean(.)

# pct_tax_ID(MBC_noNSW)
# pct_tax_ID(MBC_noPark)
# pct_tax_ID(mbc_results_MAPLE)

# Frac with family IDs and Frac with Genus/Species IDs



# Data Analysis ####
## Biomass correlations ####
# Specimens
biomass_specimens_lm <- lm(specimen_count ~ sample_wet_weight, data = collections_sample_info)
summary(biomass_specimens_lm)
visreg::visreg(biomass_specimens_lm)
cor.test(x = collections_sample_info$sample_wet_weight, y = collections_sample_info$specimen_count, method = "p")
pspearman::spearman.test(collections_sample_info$sample_wet_weight, collections_sample_info$specimen_count)

# Short / Long Duration ####
collections_sample_info %>% 
  group_by(duration) %>% 
  summarise(mean_specimen_count = mean(specimen_count),
            sd_specimen_count = sd(specimen_count))

PHAUS_res_BOLDistilled %>% 
  mutate(duration = ifelse(fieldid %in% which_samples_longDur, "17 days", "8 days")) %>% 
  select(fieldid, duration, bin_uri) %>%
  unique() %>% 
  group_by(fieldid, duration) %>% 
  summarise(n_BINs = n()) %>% 
  group_by(duration) %>% 
  summarise(mean_BIN_count = mean(n_BINs),
            sd_BIN_count = sd(n_BINs))

# Filtering Parameters ####
## False Pos/Negs ####
BC_full_BINs_samples <- PHAUS_res_BOLDistilled %>% 
  filter(`BIN Match?` == "BIN_MATCH") %>% 
  select(fieldid, bin_uri) %>% 
  rename(sample = fieldid) %>% 
  unique() %>% 
  na.omit() %>% 
  mutate(in_BC = TRUE) %>% 
  mutate(sample_BIN = paste(sample, bin_uri, sep = "_"))

MBC_full_BINs_samples <- mbc_results_MAPLE %>% 
  filter(`BIN Match Status` == "BIN_MATCH") %>% 
  select(Sample, bin_uri) %>% 
  rename(sample = Sample) %>% 
  unique() %>% 
  na.omit() %>% 
  mutate(in_MBC = TRUE)

# If a BIN was detected by MBC but not by BC, it is a FALSE POSITIVE
# I.e., it was detected but shouldn't have been
MBC_FALSE_POS <- MBC_full_BINs_samples %>% 
  left_join(., BC_full_BINs_samples, by = c("sample", "bin_uri")) %>% 
  mutate(error_type = ifelse(is.na(in_BC)==T, "FALSE_POS", NA)) %>% 
  select(sample, bin_uri, error_type) %>% 
  na.omit()
  

# If a BIN was detected by BC but not by MBC, it is a FALSE NEGATIVE
# I.e., it was not detected but should have been
MBC_FALSE_NEG <- BC_full_BINs_samples %>% 
  left_join(., MBC_full_BINs_samples, by = c("sample", "bin_uri")) %>% 
  mutate(error_type = ifelse(is.na(in_MBC)==T, "FALSE_NEG", NA)) %>% 
  select(sample, bin_uri, error_type) %>% 
  na.omit()
  
FALSE_POS_NEGS <- bind_rows(MBC_FALSE_POS, MBC_FALSE_NEG) %>% 
  mutate(sample_BIN = paste(sample, bin_uri, sep = "_"))
  
FALSE_POS_SAMPLE_BIN_Vec <- filter(FALSE_POS_NEGS, error_type == "FALSE_POS") %>% 
  .$sample_BIN

length(FALSE_POS_SAMPLE_BIN_Vec)

FALSE_NEG_SAMPLE_BIN_Vec <- filter(FALSE_POS_NEGS, error_type == "FALSE_NEG") %>% 
  .$sample_BIN

length(FALSE_NEG_SAMPLE_BIN_Vec)

## Filtering parameters and false positive/negative ####
filt_reads <- 1:100
filt_reps <- 1:8

all_filt_combs <- expand_grid(reads = filt_reads,
                              reps = filt_reps)

# Make dataframes. We want to see false pos/neg and total BINs
NUM_TRUE_POS_Vec <- NULL
NUM_FALSE_NEG_Vec <- NULL
NUM_FALSE_POS_Vec <- NULL
ACCURACY_SCORE_Vec <- NULL
PCT_BINs_SHARED_Vec <- NULL

for (i in 1:nrow(all_filt_combs)) {
  
  MBC_full_filtered_inLoop <- mbc_results_MAPLE %>% 
    filter(Reads >= all_filt_combs$reads[i],
           Replicates >= all_filt_combs$reps[i]) %>% 
    select(Sample, bin_uri) %>% 
    na.omit() %>% 
    select(Sample, bin_uri) %>% 
    rename(sample = Sample) %>% 
    unique() %>% 
    na.omit() %>% 
    mutate(in_MBC = TRUE) %>% 
    mutate(sample_BIN = paste(sample, bin_uri, sep = "_"))

  # Calculate FALSE NEGATIVES
  MBC_FALSE_NEG_InLoop <- BC_full_BINs_samples %>% 
    left_join(., MBC_full_filtered_inLoop, by = c("sample", "bin_uri")) %>% 
    mutate(error_type = ifelse(is.na(in_MBC)==T, "FALSE_NEG", NA)) %>% 
    select(sample, bin_uri, error_type) %>% 
    na.omit()
  
  # Calculate FALSE POSITIVES
  MBC_FALSE_POS_InLoop <- MBC_full_filtered_inLoop %>% 
    left_join(., BC_full_BINs_samples, by = c("sample", "bin_uri")) %>% 
    mutate(error_type = ifelse(is.na(in_BC)==T, "FALSE_POS", NA)) %>% 
    select(sample, bin_uri, error_type) %>% 
    na.omit()
  
  NUM_TRU_POS <- sum(MBC_full_filtered_inLoop$sample_BIN %in% BC_full_BINs_samples$sample_BIN)
  NUM_FALSE_NEG <- nrow(MBC_FALSE_NEG_InLoop)
  NUM_FALSE_POS <- nrow(MBC_FALSE_POS_InLoop)
  ACCURACY_SCORE = NUM_TRU_POS - NUM_FALSE_NEG - NUM_FALSE_POS
  PCT_BINs_SHARED <- mean(unique(MBC_full_filtered_inLoop$bin_uri) %in% PHAUS_res$bin_uri)
    
  # Add data 
  NUM_TRUE_POS_Vec <- c(NUM_TRUE_POS_Vec, NUM_TRU_POS)
  ACCURACY_SCORE_Vec <- c(ACCURACY_SCORE_Vec, ACCURACY_SCORE)
  NUM_FALSE_NEG_Vec <- c(NUM_FALSE_NEG_Vec, NUM_FALSE_NEG)
  NUM_FALSE_POS_Vec <- c(NUM_FALSE_POS_Vec, NUM_FALSE_POS)
  PCT_BINs_SHARED_Vec <- c(PCT_BINs_SHARED_Vec, PCT_BINs_SHARED)
  
}

Data_For_Fig2_Plot <- all_filt_combs %>% 
  mutate(NUM_TRUE_POS = NUM_TRUE_POS_Vec,
         NUM_FALSE_NEG = NUM_FALSE_NEG_Vec,
         NUM_FALSE_POS = NUM_FALSE_POS_Vec,
         ACCURACY_SCORE = ACCURACY_SCORE_Vec,
         NUM_ERRORS = NUM_FALSE_NEG + NUM_FALSE_POS,
         PCT_BINs_SHARED = PCT_BINs_SHARED_Vec)

# PHAUS Barcode Diversity ####
# Number of BINs
BC_full %>% 
  filter(`BIN Match?` == "BIN_MATCH") %>% 
  select(BIN) %>% 
  unique %>% nrow(.)

# Abundnace of BINs
PHAUS_BIN_count_taxa <- BC_full %>%
  drop_na(BIN) %>%
  group_by(BIN, class, order, family, genus, species) %>% 
  mutate(n_records = n()) %>% 
  mutate(row_number = row_number()) %>% 
  filter(row_number == 1) %>% 
  select(-row_number) %>% 
  select(BIN, n_records, class, order, family, genus, species)

table(PHAUS_BIN_count_taxa$n_records)

## Metabarcoding Summary ####
# make options to add filters?
MBC_full_filt <- filter(MBC_full, Reads >= 25, Replicates >=2)
# Need to mess around with 

table(MBC_full_filt$`BIN Match?`)

# What queries to keep

unique_ASVs_full <- MBC_full %>% 
    filter(kingdom == "Animalia") %>% 
    filter(Replicates >= 7) %>% 
    select(ASV_ID, bin_uri, BIN_Match) %>% 
  group_by(ASV_ID, bin_uri) %>% 
  unique() %>% 
  # drop ASVs appearing multiple times
  group_by(ASV_ID) %>% 
  arrange(BIN_Match) %>% 
  mutate(num_in_group = row_number()) %>% 
  ungroup() %>% 
  filter(num_in_group == 1)

unique_ASVs_noAus <- MBC_noAus %>% 
  filter(ASV_ID %in% unique_ASVs_full$ASV_ID) 

unique_ASVs_noNSW <- MBC_noNSW %>% 
  filter(ASV_ID %in% unique_ASVs_full$ASV_ID) 

unique_ASVs_noMua <- MBC_noPark %>% 
  filter(ASV_ID %in% unique_ASVs_full$ASV_ID) 

full_BM_table <- table(unique_ASVs_full$BIN_Match)
noMua_BM_table <- table(unique_ASVs_noMua$BIN_Match)
noNSW_BM_table <- table(unique_ASVs_noNSW$BIN_Match)
noAus_BM_table <- table(unique_ASVs_noAus$BIN_Match)

ref_tot <- nrow(MBC_full)

# Small bug, fix...
sequential_dat <- bind_rows(full_BM_table, noMua_BM_table, noNSW_BM_table, noAus_BM_table) %>% 
  mutate(BD_lib = c("complete", "no park", "no NSW", "no Aus")) %>% 
  mutate(tot = BIN_MATCH + NO_MATCH) %>% 
  mutate(
    NO_MATCH = NO_MATCH + (ref_tot - tot),
    new_tot = ref_tot
  )

sequential_dat_improvement <- sequential_dat %>% 
  arrange(BIN_MATCH) %>% 
  mutate(pct_match = BIN_MATCH / (BIN_MATCH + NO_MATCH),
         additional_pct = abs(pct_match - lag(pct_match)))

# Taxonomic mismatch ####
## Get two datasets of all unique BINs ####
unique_tax_mbc <- mbc_results_MAPLE %>% 
  filter(`BIN Match Status` == "BIN_MATCH") %>% 
  rename(sample = Sample) %>% 
  drop_na(bin_uri) %>%
  drop_na(Phylum) %>% 
  # filter(Reads >= 25) %>% 
  select(bin_uri, Kingdom:Species) %>% 
  unique() 

in_MBC_df <- unique_tax_mbc %>% 
  select(bin_uri) %>% 
  mutate(det_MBC = T)

unique_tax_bc <- PHAUS_res_BOLDistilled %>% 
  rename(sample = fieldid) %>% 
  drop_na(bin_uri) %>%
  drop_na(phylum) %>% 
  select(bin_uri, kingdom:species) %>% 
  unique()

in_BC_df <- unique_tax_bc %>% 
  select(bin_uri) %>% 
  mutate(det_BC = T)

## Get ONLY BINs and Tax DF ####
BINs_det_by_method <- bind_rows(unique_tax_mbc, unique_tax_bc) %>% 
  unique() %>% 
  left_join(., in_BC_df, by = "bin_uri") %>% 
  left_join(., in_MBC_df, by = "bin_uri") %>% 
  relocate(det_BC, det_MBC, .after = bin_uri) %>% 
  mutate(
    det_BC = if_else(is.na(det_BC), FALSE, det_BC),
    det_MBC = if_else(is.na(det_MBC), FALSE, det_MBC)
  )

BINs_total_by_order <- BINs_det_by_method %>% 
  group_by(class, order) %>% 
  summarise(n_BINs_det_in_order = n())

BINs_total_by_family <- BINs_det_by_method %>% 
  select(bin_uri, Kingdom:Species) %>% 
  unique() %>% 
  group_by(Order, Family) %>% 
  summarise(n_BINs_det_in_fam = n()) %>% 
  rename(order = Order, family = Family)

write_tsv(BINs_det_by_method, 'data/derived_data/BINs_det_by_method.tsv')
  
## Missed - BINs ##
missed_by_BC_BINs <- BINs_det_by_method %>% 
  filter(det_BC == F)

missed_by_MBC_BINs <- BINs_det_by_method %>% 
  filter(det_MBC == F)

## Missed - Orders ##
missed_by_BC_orders <- missed_by_BC_BINs %>% 
  select(1:10) %>% 
  group_by(Class, Order) %>% 
  summarise(n_missed = n()) %>% 
  rename(class = Class, order = Order) %>%
  arrange(desc(n_missed)) %>% 
  left_join(., BINs_total_by_order, by = c("class", "order")) %>% 
  na.omit() %>% 
  mutate(prop_missed = n_missed / n_BINs_det_in_order) %>% 
  mutate(missed_by = "barcode")
  
missed_by_MBC_orders <- missed_by_MBC_BINs %>% 
  group_by(class, order) %>% 
  summarise(n_missed = n()) %>% 
  arrange(desc(n_missed)) %>% 
  left_join(., BINs_total_by_order, by = c("class", "order")) %>% 
  na.omit() %>% 
  mutate(prop_missed = n_missed / n_BINs_det_in_order) %>% 
  mutate(missed_by = "metabarcode")

## Missed - Family ##
missed_by_BC_family <- missed_by_BC_BINs %>% 
  select(bin_uri:Species) %>% 
  rename(class = Class, order = Order, family = Family) %>%
  group_by(order, family) %>% 
  summarise(n_missed = n()) %>% 
  arrange(desc(n_missed)) %>% 
  left_join(., BINs_total_by_family, by = c("order", "family")) %>% 
  na.omit() %>% 
  mutate(prop_missed = n_missed / n_BINs_det_in_fam) %>% 
  mutate(missed_by = "barcode")

missed_by_MBC_family <- missed_by_MBC_BINs %>% 
  group_by(order, family) %>% 
  summarise(n_missed = n()) %>% 
  arrange(desc(n_missed)) %>% 
  left_join(., BINs_total_by_family, by = c("order", "family")) %>% 
  na.omit() %>% 
  mutate(prop_missed = n_missed / n_BINs_det_in_fam) %>% 
  mutate(missed_by = "metabarcode")

## Data for Taxonomy figure ####
missing_tax_df <- bind_rows(missed_by_BC_orders, missed_by_MBC_orders) 

# Step 1: Get total BINs detected per order
bins_per_order <- missing_tax_df %>%
  distinct(class, order, n_BINs_det_in_order)

# Step 2: Get number missed by each method
missed_summary <- missing_tax_df %>%
  select(order, missed_by, n_missed) %>%
  pivot_wider(names_from = missed_by, values_from = n_missed, values_fill = 0)

# Step 3: Combine and compute detected_by_both
stacked_df <- bins_per_order %>%
  left_join(missed_summary, by = c("order", "class")) %>%
  mutate(
    detected_by_both = n_BINs_det_in_order - barcode - metabarcode,
    commonness = if_else(n_BINs_det_in_order >= 50, "Major Orders (≥50 BINs)", "Other Orders (<50 BINs)")
  ) %>%
  pivot_longer(cols = c(metabarcode, barcode, detected_by_both),
               names_to = "detection_status", values_to = "n") %>%
  mutate(
    order = fct_reorder(order, n_BINs_det_in_order),
    detection_status = factor(detection_status, levels = c("metabarcode", "barcode", "detected_by_both"))
  )

# Hymenoptera difference
PHAUS_Hymenoptera_Det <- PHAUS_res_BOLDistilled %>% 
  filter(order == "Hymenoptera") %>% 
  select(bin_uri, family) %>% 
  unique() %>% 
  count(family) %>% 
  filter(n >= 10)%>% 
  rename(n_BC = n)

MAPLE_Hymenoptera_Det <- mbc_results_MAPLE %>% 
  filter(Order == "Hymenoptera") %>% 
  select(bin_uri, Family) %>% 
  unique() %>% 
  count(Family) %>% 
  filter(n >= 10) %>% 
  rename(n_MBC = n, 
         family = Family)

Hymenoptera_Compare_DF <- left_join(PHAUS_Hymenoptera_Det, MAPLE_Hymenoptera_Det, by = "family") %>% 
  mutate(factor_off = abs(log10(n_BC / n_MBC)))

# Detection Probability ####
unique_tax_bc_tally <- PHAUS_res_BOLDistilled %>% 
  rename(sample = fieldid) %>% 
  drop_na(bin_uri) %>%
  drop_na(phylum) %>% 
  select(bin_uri, kingdom:species) %>% 
  group_by(bin_uri) %>% 
  summarise(n_barcoded = n())

MBC_Detections <- mbc_results_MAPLE %>% 
  select(bin_uri) %>% 
  unique() %>% 
  mutate(det_MBC = T)

det_prob <- unique_tax_bc_tally %>% 
  left_join(., MBC_Detections, by = "bin_uri") %>% 
  mutate(det_MBC = ifelse(is.na(det_MBC)==T, F, T))

det_prob_binned <- det_prob %>%
  mutate(n_barcoded_bin = case_when(
    n_barcoded == 1 ~ "1",
    n_barcoded == 2 ~ "2",
    n_barcoded >= 3  & n_barcoded <= 6   ~ "3–6",
    n_barcoded >= 7  & n_barcoded <= 10  ~ "7–10",
    n_barcoded >= 11 & n_barcoded <= 20  ~ "11–20",
    n_barcoded >= 21 & n_barcoded <= 50  ~ "21–50",
    n_barcoded >= 51 & n_barcoded <= 100 ~ "51–100",
    n_barcoded > 100 ~ ">100"
  )) %>% 
  group_by(n_barcoded_bin) %>%
  summarise(
    n = n(),
    n_true = sum(det_MBC),
    frac_true = mean(det_MBC),
    .groups = "drop"
  ) %>%
  arrange(factor(n_barcoded_bin, levels = c("1", "2", "3–5", "6–10", "11–20", "21–50", "51–100", ">100")))

## Logistic regression model ####
det_prob_model <- glm(det_MBC ~ n_barcoded, data = det_prob, family = binomial(link = "logit"))
summary(det_prob_model)
visreg::visreg(det_prob_model)



# DARK BIN Filtering ####

# what are the ASVs that got 'resolved' by barcoding the Traps?
# These are the ASVs that did vs. didn't have a BIN_Match before barcoding the traps
MBC_noPark_BIN_MatchOnly <- MBC_noPark %>% 
  select(ASV_Name, `BIN Match?`) %>% 
  rename(BIN_Match_noTraps = `BIN Match?`) %>% 
  unique() 

MBC_ASVs_NewFromMuo <- MBC_full %>% 
  select(ASV_Name, `BIN Match?`) %>% 
  rename(BIN_Match_Full = `BIN Match?`) %>% 
  unique() %>% 
  filter(BIN_Match_Full == "BIN_MATCH") %>% 
  # this is the ones that had a BIN match AFTER barcoding the traps
  left_join(., MBC_noTraps_BIN_MatchOnly, by = "ASV_Name") %>% 
  filter(BIN_Match_Full == "BIN_MATCH" & BIN_Match_noTraps == "NO_MATCH") %>% 
  select(ASV_Name) %>% 
  unique()
# The ASVs above are the 'REAL' Dark BINs,
# which gained a BIN match after barcoding traps.   

# We want to know: OF ALL THE ASVs identified, how many were BIN HITS (real) vs not (likely false)
# Consider changing to take 'Full vs ref' in case we want to swap between park and sample exclusion.
All_ASVs_BIN_Match_noMuo <- MBC_noPark %>% 
  select(ASV_Name, `BIN Match?`) %>% 
  group_by(ASV_Name) %>% 
  summarise(ASV_Matches_BIN = any(`BIN Match?` == "BIN_MATCH"))

MBC_noPark_ASV_Summary <- MBC_noPark %>% 
  group_by(ASV_Name) %>% 
  summarise(max_reps = max(Replicates),
            tot_reads = sum(Reads),
            n_samples = length(unique(Sample))) %>% 
  ungroup() %>% 
  left_join(., All_ASVs_BIN_Match_noMuo, by = "ASV_Name") #%>% 
  # filter(tot_reads >= 25)

n_Dark_BIN_ASVs <- NULL
prop_Dark_BIN_ASVs <- NULL

for (i in 1:8) {

  ASV_DF_2samp <- MBC_noPark_ASV_Summary %>% 
    filter(ASV_Matches_BIN == F) %>% 
    # Drop the ones that would have already matched a BIN.
    # the remaining are Candidate DARK BINs
    # Dark BINs are those ASVs that have no BIN in the noMuo lib
    # But do have a BIN in the full lib.
    mutate(enough_reps = max_reps >= i,
           enough_samples = n_samples >= 2) %>% 
    filter(enough_reps == T | enough_samples == T) 
  
  ASV_DF_1samp <- MBC_noPark_ASV_Summary %>% 
    # Drop the ones that would have already matched a BIN.
    # the remaining are Candidate DARK BINs
    # Dark BINs are those ASVs that have no BIN in the noMuo lib
    # But do have a BIN in the full lib.
    filter(ASV_Matches_BIN == F) %>% 
    mutate(enough_reps = max_reps >= i) %>% 
  filter(enough_reps == T)
  
  ASVs_for_RAMS_1samp <- ASV_DF_1samp$ASV_Name
  
  n_ASVs_for_RAMS_1samp <- length(ASVs_for_RAMS_1samp) # Number of ASVs that would go to RAMS
  n_TRUE_ASVs_1samp <- sum(ASVs_for_RAMS_1samp %in% MBC_ASVs_NewFromMuo$ASV_Name) # Number of RAMS-bound ASVs that are new BINs
  n_FALSE_ASVs_1samp <- n_ASVs_for_RAMS_1samp - n_TRUE_ASVs_1samp
  prop_ASVs_true_1samp <- n_TRUE_ASVs_1samp / n_ASVs_for_RAMS_1samp
  
  # final vectors
  n_Dark_BIN_ASVs <- c(n_Dark_BIN_ASVs, n_ASVs_for_RAMS_1samp)
  prop_Dark_BIN_ASVs <- c(prop_Dark_BIN_ASVs, prop_ASVs_true_1samp)
          
}

# Fig X: Dark BIN filtering
Fig_S_DarkBIN_DF <- tibble(n_reps = sort(rep(1:8)),
                           n_Dark_BIN_ASVs, 
                           prop_Dark_BIN_ASVs)

# Prepare data for stacked bar plot
# Prepare stacked data with correct ordering
# Prepare stacked data with correct ordering
DarkBIN_plot_data <- Fig_S_DarkBIN_DF %>%
  mutate(
    dark_bin = n_Dark_BIN_ASVs * prop_Dark_BIN_ASVs,
    non_dark_bin = n_Dark_BIN_ASVs - dark_bin) %>% 
  select(n_reps, dark_bin, non_dark_bin) %>%
  pivot_longer(cols = c(dark_bin, non_dark_bin), names_to = "type", values_to = "count") %>%
  mutate(
    type = factor(type, levels = c("non_dark_bin", "dark_bin"))
  ) %>%
  arrange(n_reps, type)

# proportions
DarkBIN_PCT <- DarkBIN_plot_data %>% 
  group_by(n_reps) %>% 
  mutate(pct = count / sum(count) * 100)
  

# Output Dark BINs ####
Dark_BIN_DF <- MBC_noPark_ASV_Summary %>% 
  filter(ASV_Matches_BIN == F) %>% 
  mutate(Dark_BIN = ASV_Name %in% MBC_ASVs_NewFromMuo$ASV_Name) %>% 
  mutate(ASV_type = ifelse(Dark_BIN == T, "Dark_BIN", "Artifact")) %>% 
  select(ASV_Name, ASV_type) 


# Alpha and Beta Diversity ####
BC_Alpha_BIN <- PHAUS_res_BOLDistilled %>% 
  drop_na(bin_uri) %>% 
  select(fieldid, bin_uri) %>% 
  unique() %>% 
  count(fieldid) %>% 
  rename(Sample = fieldid) %>% 
  rename(BIN_Count_BC = n)

MBC_Alpha_BIN <- mbc_results_MAPLE %>% 
  filter(`BIN Match Status` == "BIN_MATCH") %>%
  filter(Reads >= 12 | Replicates >= 2) %>%
  # filter(
  #   Reads >= 250 |
  #     Replicates >= 6 |
  #     (`BIN Match Status` == "BIN_MATCH" & (Reads >= 25 | Replicates >= 2))
  # ) %>% 
  drop_na(bin_uri) %>% 
  select(Sample, bin_uri) %>% 
  unique() %>% 
  count(Sample) %>% 
  rename(BIN_Count_MBC = n)

Alpha_Comp_DF <- left_join(BC_Alpha_BIN, MBC_Alpha_BIN) %>% 
  mutate(duration = ifelse(Sample %in% which_samples_shortDur, "short", "long"))

summary(lm(BIN_Count_MBC ~ BIN_Count_BC, Alpha_Comp_DF))
cor.test(Alpha_Comp_DF$BIN_Count_BC, Alpha_Comp_DF$BIN_Count_MBC)
  

## Beta ####
MBC_for_PPM <- mbc_results_MAPLE %>% 
  filter(`BIN Match Status` == "BIN_MATCH") %>%
  filter(Reads >= 25 | Replicates >= 2) %>%
  # filter(
  #   Reads >= 250 |
  #     Replicates >= 6 |
  #     (`BIN Match Status` == "BIN_MATCH" & (Reads >= 25 | Replicates >= 2))
  # ) %>% 
  rename(fieldid = Sample,
         family = Family) %>% 
  mutate(method = "Metabarcode") %>% 
  select(fieldid, family, bin_uri, method)

Data_for_PPM <- PHAUS_res_BOLDistilled %>% 
  select(fieldid, family, bin_uri) %>% 
  mutate(method = "Barcode") %>% 
  bind_rows(., MBC_for_PPM) %>% 
  mutate(unique_sample = paste0(method, "_", fieldid))

# Beta mantel #
Barcode_df <- Data_for_PPM %>% 
  filter(method == "Barcode")

Metabarcode_df <- Data_for_PPM %>% 
  filter(method == "Metabarcode")

# Matrix
barcode_mat <- Barcode_df %>% 
  select(fieldid, bin_uri) %>%
  mutate(value = 1) %>%
  distinct() %>%                        # prevents duplicates
  pivot_wider(names_from = bin_uri,
              values_from = value,
              values_fill = 0) %>%
  arrange(fieldid)

metabarcode_mat <- Metabarcode_df %>% 
  select(fieldid, bin_uri) %>%
  mutate(value = 1) %>%
  distinct() %>%
  pivot_wider(names_from = bin_uri,
              values_from = value,
              values_fill = 0) %>%
  arrange(fieldid)

# Filter samples 
shared_ids <- intersect(barcode_mat$fieldid, metabarcode_mat$fieldid)

barcode_mat <- barcode_mat %>% 
  filter(fieldid %in% shared_ids) %>% 
  arrange(fieldid)

metabarcode_mat <- metabarcode_mat %>% 
  filter(fieldid %in% shared_ids) %>% 
  arrange(fieldid)

# Drop Fieldid
barcode_comm <- barcode_mat %>% select(-fieldid)
metabarcode_comm <- metabarcode_mat %>% select(-fieldid)

# Distance matrix
D_barcode <- vegdist(barcode_comm, method = "jaccard", binary = TRUE)
D_meta    <- vegdist(metabarcode_comm, method = "jaccard", binary = TRUE)

# Mantel
mantel_out <- mantel(D_barcode, D_meta, method = "pearson", permutations = 10000)

mantel_out

# Procrustes
ord_bar  <- cmdscale(D_barcode, k = 2)
ord_meta <- cmdscale(D_meta, k = 2)

proc <- procrustes(ord_bar, ord_meta)

protest_out <- protest(ord_bar, ord_meta, permutations = 999)

# Print results
mantel_out
protest_out


# BIN matching and Tax Res ####
mbc_results_MAPLE_AllTax <- mbc_results_MAPLE %>% 
  select(ASV_ID, bin_uri, `BIN Match Status`, Kingdom:Species, `BIN Kingdom`:`BIN Species`) %>% 
  group_by(ASV_ID) %>% 
  slice(1)

# Tax match
(1 - mbc_results_MAPLE_AllTax$Family %>% is.na(.) %>% mean(.))*100
(1 - mbc_results_MAPLE_AllTax$`BIN Genus` %>% is.na(.) %>% mean(.))*100
(1 - mbc_results_MAPLE_AllTax$Species %>% is.na(.) %>% mean(.))*100


mbc_results_MAPLE_BINsOnlyTax <- mbc_results_MAPLE %>% 
  filter(`BIN Match Status` == "BIN_MATCH") %>% 
  select(ASV_ID, bin_uri, `BIN Match Status`, Kingdom:Species, `BIN Kingdom`:`BIN Species`) %>% 
  group_by(ASV_ID) %>% 
  slice(1)

# Family match
(1 - mbc_results_MAPLE_BINsOnlyTax$Family %>% is.na(.) %>% mean(.))*100
(1 - mbc_results_MAPLE_BINsOnlyTax$`BIN Genus` %>% is.na(.) %>% mean(.))*100
(1 - mbc_results_MAPLE_BINsOnlyTax$Species %>% is.na(.) %>% mean(.))*100


# Were the BINs missed by barcoding known from Australia?

missed_by_BC_BIN_vec <- missed_by_BC_BINs %>% 
  filter(Kingdom == "Animalia") %>% 
  .$bin_uri



write_lines(missed_by_BC_BIN_vec, file = 'data/derived_data/BINs_missed_by_BC.txt')

n_BINs_missed_BC <- length(missed_by_BC_BIN_vec)

missed_BINs <- bold.fetch(get_by = "bin_uris", identifiers = missed_by_BC_BIN_vec, filt_geography = "Australia")





n_missed_known_Aus <-  unique(missed_BINs$bin_uri) %>% length(.)
# pct known from australis

n_missed_known_Aus / n_BINs_missed_BC


# Specimens without seqs 
# Seqs without BINs
PHAUS_res_badSeq <- PHAUS_res %>% 
  filter(is.na(nuc)==T & is.na(nuc_basecount) == F)
  
write_tsv(PHAUS_res_badSeq, 'data/derived_data/PHAUS_res_badSeq.tsv')
write_tsv(PHAUS_res_noSeq, 'data/derived_data/PHAUS_res_noSeq.tsv')
write_tsv(missed_BINs, 'data/derived_data/missed_BINs.tsv')
