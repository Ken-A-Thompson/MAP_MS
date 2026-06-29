# MAP Paper ####
# Read & Clean Data #

# load packages #
library(BOLDconnectR)
library(DescTools)
library(flextable)
library(eulerr)
library(readxl)
library(tidyverse)
library(vegan)

collections_sample_data <- read_excel('data/collections_data/CBGMB01277 sample data.xlsx')

# Download PHAUS BOLD Data ####
# If downloaded already, can skip to LOAD 

# bold.apikey('00FE0FAE-94A2-4E94-9224-4CB34BC96A91')

# PHAUS_res <- bold.fetch(get_by = "project_codes", identifiers = "PHAUS") # complete on internet

# PHAUS_CBGMB01277 <- filter(PHAUS_res, fieldid %in% collections_sample_data$Lot_fieldID)
# write_tsv(PHAUS_CBGMB01277, 'data/BOLD_BCDM/PHAUS_CBGMB01277.tsv') 

# Read in all data ####
PHAUS_BOLD <- read_tsv('data/BOLD_BCDM/PHAUS_CBGMB01277.tsv') 

PHAUS_BOLDistilled <- read_tsv('data/BOLDistilled_ID/BOLDistilled_ID_PHAUS_BOLD_Jan2025/BIN_MATCH_20260414_121531.tsv') %>% 
  select(Query, BIN, kingdom:species) %>% 
  rename(processid = Query,
         bin_uri = BIN)

PHAUS_BOLD_Clean <- PHAUS_BOLD %>% 
  select(processid, fieldid) %>% 
  left_join(., PHAUS_BOLDistilled) %>% 
  select(-subfamily, -tribe) %>% 
  mutate(method = "BOLD") %>% 
  drop_na(bin_uri)

PHAUS_BOLD_Clean_NTS <- read_excel('data/BOLD_BCDM/PHAUS_parent_nts_BCDM.xlsx') %>%
  rename(bin_uri = bin) %>%
  drop_na(bin_uri) %>%
  filter(kingdom == "Animalia") %>%
  filter(is.na(flagged)==T)

# Vector of all BINs on BOLD
PHAUS_BOLD_BINs <- unique(PHAUS_BOLD_Clean_NTS$bin_uri)

# Make Fasta #
# PHAUS_BOLD %>% 
  # drop_na(bin_uri) %>% 
  # mutate(fasta_entry = paste0(">", processid, "\n", nuc)) %>%
  # pull(fasta_entry) %>%
  # writeLines("data/BOLDistilled_ID/PHAUS_BOLD.fa")

# BOLD summary stats ####
# How many specimens
N_Specimens_PHAUS <- nrow(PHAUS_BOLD_Clean)

# How many had BINs ####
N_Specimens_w_BINs <- na.omit(PHAUS_BOLD_Clean$bin_uri) %>% length(.)
N_Specimens_w_BINs /  N_Specimens_PHAUS

# Look at also mean sample gamma

# Load Data ####
## MAP ONT ####
PHAUS_MBC_ONT_MAP <- read_tsv('data/MAP_output/PHAUS_ONT.MAP2026_03_24/2-TSV Versions of Results/Metabarcoding_Results_PHAUS_COI-5P_658_BySample.tsv') %>% 
  rename(fieldid = Sample, 
         bin_uri = BIN_Hit) %>% 
  filter(fieldid %in% collections_sample_data$Lot_fieldID) %>% 
  rename_with(tolower, Kingdom:Species) %>% 
  rename(tot_reads = Reads,
         replicates = Replicates) %>% 
  select(fieldid, tot_reads, ASV_ID, replicates, kingdom:bin_uri) %>%
  drop_na(bin_uri)

## MAP ILL ####
PHAUS_MBC_ILL_MAP <- readxl::read_excel('data/MAP_output/PHAUS_ILL.MAP2026-05-26/Metabarcoding Results - PHAUS_COI-5P_418.xlsx') %>%
  rename(fieldid = Sample,
         bin_uri = `BIN Hit`) %>%
  mutate(fieldid = str_replace_all(fieldid, "-", "#")) %>%
  filter(fieldid %in% collections_sample_data$Lot_fieldID) %>%
  rename_with(tolower, Kingdom:Species) %>%
  rename(tot_reads = Reads) %>%
  rename(replicates = Replicates) %>% 
  select(fieldid, tot_reads, OTU_ID, replicates, kingdom:bin_uri) %>%
  filter(kingdom == "Animalia")

## Spcfy ####
PHAUS_MBC_ILL_spcfy <- readxl::read_xlsx('data/benchmarking/spcfy/PHAUS_Illumina_Benchmarking_SPARK.xlsx', skip = 1) %>% 
  mutate(OTU_NUMBER = sub(";.*", "", OTU_ID))

PHAUS_MBC_ILL_spcfy %>%
  drop_na(BOLD_BIN_uri) %>%
  mutate(fasta_entry = paste0(">", OTU_NUMBER, "\n", OTU_fasta_sequence)) %>%
  pull(fasta_entry) %>%
  writeLines("data/BOLDistilled_ID/spcfy.fa")

PHAUS_SPCFY_BOLDistilled <- read_tsv('data/BOLDistilled_ID/BOLDistilled_ID_Spcfy_Jan2025/BIN_MATCH_20260415_095915.tsv') %>% 
  select(Query, BIN, kingdom:species) %>% 
  rename(OTU_NUMBER = Query,
         bin_uri = BIN)

PHAUS_MBC_ILL_SPCFY_OTU_BOLDTax <- PHAUS_MBC_ILL_spcfy %>% 
  select(OTU_ID, OTU_NUMBER) %>% 
  left_join(., PHAUS_SPCFY_BOLDistilled, by = "OTU_NUMBER")

PHAUS_MBC_SPCFY_long <- PHAUS_MBC_ILL_spcfy %>% 
  select(BOLD_Process_ID:OTU_ID, GMP_58219_Rep1:GMP_58228_Rep8) %>% 
  select(-`BOLD_Grade%ID`) %>% 
  pivot_longer(GMP_58219_Rep1:GMP_58228_Rep8, names_to = "sample", values_to = "reads") %>% 
  filter(!reads == 0) %>% 
  separate(sample, into = c("fieldid", "rep"), sep = "_Rep", remove = T) %>% 
  group_by(fieldid, OTU_ID) %>% 
  summarise(replicates = n(),
            tot_reads = sum(reads)) %>% 
  left_join(., PHAUS_MBC_ILL_SPCFY_OTU_BOLDTax, by = "OTU_ID") %>% 
  mutate(fieldid = gsub("_", "#", fieldid)) 

PHAUS_MBC_ILL_SPCFY <- PHAUS_MBC_SPCFY_long %>%
  select(fieldid, replicates, tot_reads, bin_uri, phylum:species) %>%
  ungroup() %>%
  mutate(species = gsub("_", " ", species)) %>%
  drop_na(bin_uri)

## mBRAVE ####
PHAUS_MBC_ILL_mBRAVE_Raw <- read_tsv('data/benchmarking/mBRAVE/All_Sets_of_Data_Illumina.tsv')

PHAUS_ILL_mBRAVE <- PHAUS_MBC_ILL_mBRAVE_Raw %>% 
  select(sampleId, bin_uri, sequences, phylum:bin_uri) %>% 
  separate(sampleId, into = c('prefix', 'fieldid', 'rep'), sep = "_") %>% 
  mutate(fieldid = gsub("-", "#", fieldid)) %>% 
  filter(fieldid %in% collections_sample_data$Lot_fieldID) %>% 
  select(-prefix) %>% 
  group_by(fieldid, phylum, class, order, family, genus, species, bin_uri) %>% 
  summarise(tot_reads = sum(sequences),
            replicates = length(unique(rep))) %>% 
  mutate(method = "mBRAVE_ILL") %>%
  ungroup()

# Evaluate Performance ####
## MAP ONT ####

# 1. Determine filtering
PHAUS_MAP_ONT_ABG <- compute_a_b_g(MBC_data = PHAUS_MBC_ONT_MAP, GT_data = PHAUS_BOLD_Clean_NTS)

PHAUS_MAP_ONT_BEST <- slice(PHAUS_MAP_ONT_ABG, which.min(rank_sum))
PHAUS_MAP_ONT_BEST <- slice(PHAUS_MAP_ONT_ABG, which(tot_reads == min(tot_reads) & as.integer(replicates) == min(as.integer(replicates))))
PHAUS_MAP_ONT_BEST <- compute_MinFilt_Metrics(
  MBC_data     = PHAUS_MBC_ONT_MAP,
  ground_truth = PHAUS_BOLD_Clean_NTS,
  threshold_factor = 1e-05
)

# 2a Summary Stats
PHAUS_MAP_ONT_Filt <- PHAUS_MBC_ONT_MAP %>% 
  filter(
    tot_reads >= PHAUS_MAP_ONT_BEST$tot_reads[1],  # dynamic threshold
         replicates >= as.integer(PHAUS_MAP_ONT_BEST$replicates),  # dynamic threshold
         # tot_reads >= min(tot_reads), replicates >= min(replicates)  # minimum (keep all)
  ) %>% 
  mutate(method = "MAP_ONT")



## MAP ILL ####
PHAUS_MAP_ILL_ABG <- compute_a_b_g(MBC_data = PHAUS_MBC_ILL_MAP, GT_data = PHAUS_BOLD_Clean_NTS)

PHAUS_MAP_ILL_BEST <- slice(PHAUS_MBC_ILL_ABG, which.min(rank_sum))
PHAUS_MAP_ILL_BEST <- slice(PHAUS_MBC_ILL_ABG, which(tot_reads == min(tot_reads) & as.integer(replicates) == min(as.integer(replicates))))
PHAUS_MAP_ILL_BEST <- compute_MinFilt_Metrics(
  MBC_data     = PHAUS_MBC_ILL_MAP,
  ground_truth = PHAUS_BOLD_Clean_NTS,
  threshold_factor = 1e-05
)

# MAP_ILL: 2a PCT BOLD Species Captured
PHAUS_MAP_ILL_Filt <- PHAUS_MBC_ILL_MAP %>% 
  filter(
    tot_reads >= PHAUS_MAP_ILL_BEST$tot_reads[1],
         replicates >= as.integer(PHAUS_MAP_ILL_BEST$replicates),
    # tot_reads >= min(tot_reads), replicates >= min(replicates)
    ) %>% 
  mutate(method = "MAP_ILL")

## SPCFY ####
PHAUS_ILL_SPCFY_ABG <- compute_a_b_g(MBC_data = PHAUS_MBC_ILL_SPCFY, GT_data = PHAUS_BOLD_Clean_NTS)

PHAUS_ILL_SPCFY_BEST <- slice(PHAUS_ILL_SPCFY_ABG, which.min(rank_sum))
PHAUS_ILL_SPCFY_BEST <- slice(PHAUS_ILL_SPCFY_ABG, which(tot_reads == min(tot_reads) & as.integer(replicates) == min(as.integer(replicates))))
PHAUS_ILL_SPCFY_BEST <- compute_MinFilt_Metrics(MBC_data = PHAUS_MBC_ILL_SPCFY, ground_truth = PHAUS_BOLD_Clean_NTS, threshold_factor = 1e-05)


# SPCY: 2a PCT BOLD Species Captured
PHAUS_SPCFY_ILL_Filt <- PHAUS_MBC_ILL_SPCFY %>% 
  filter(
    tot_reads >= PHAUS_ILL_SPCFY_BEST$tot_reads[1],
         replicates >= as.integer(PHAUS_ILL_SPCFY_BEST$replicates),
    # tot_reads >= min(tot_reads), replicates >= min(replicates)
    ) %>% 
  mutate(method = "SPCFY")

## mBRAVE ####
PHAUS_ILL_mBRAVE_ABG <- compute_a_b_g(MBC_data = PHAUS_ILL_mBRAVE, GT_data = PHAUS_BOLD_Clean_NTS)

PHAUS_mBRAVE_ILL_BEST <- slice(PHAUS_ILL_mBRAVE_ABG, which.min(rank_sum))
PHAUS_mBRAVE_ILL_BEST <- slice(PHAUS_ILL_mBRAVE_ABG, which(tot_reads == min(tot_reads) & as.integer(replicates) == min(as.integer(replicates))))
PHAUS_mBRAVE_ILL_BEST <- compute_MinFilt_Metrics(MBC_data = PHAUS_ILL_mBRAVE, ground_truth = PHAUS_BOLD_Clean_NTS, threshold_factor = 1e-05)


# MAP_ILL: 2a PCT BOLD Species Captured
PHAUS_mBRAVE_ILL_Filt <- PHAUS_ILL_mBRAVE %>% 
  filter(
    tot_reads >= PHAUS_ILL_mBRAVE_BEST$tot_reads[1],
         # replicates >= as.integer(PHAUS_ILL_mBRAVE_BEST$replicates),
    tot_reads >= min(tot_reads), replicates >= min(replicates)) %>% 
  mutate(method = "mBRAVE")

mean(PHAUS_BOLD_BINs %in% unique(PHAUS_ILL_mBRAVE_Filt$bin_uri))

# Beta Dist All Filt ####
# use this instead if you want the smart filter
all_filt_data <- bind_rows(PHAUS_BOLD_Clean_NTS,
                           PHAUS_MAP_ONT_Filt, 
                           PHAUS_MAP_ILL_Filt,
                           PHAUS_mBRAVE_ILL_Filt,
                           PHAUS_SPCFY_Filt) %>% 
  mutate(method_fieldid = paste0(method, ".", fieldid)) %>% 
  select(method_fieldid, bin_uri) %>% 
  unique() %>% 
  drop_na(method_fieldid) %>% 
  drop_na(bin_uri)

all_filt_commat <- make_comm_matrix(all_filt_data, "method_fieldid", sp_col = "bin_uri")
dist_all_filt <- vegan::vegdist(all_filt_commat, method = "bray")
nmds_all_filt <- vegan::metaMDS(dist_all_filt, k = 2, trymax = 100)

# --- scores ---
scores_all <- as.data.frame(vegan::scores(nmds_all_filt, display = "sites")) %>%
  tibble::rownames_to_column("sample") %>% 
  separate(sample, into = c("method", "fieldid"), sep = "\\.")

scores_bold <- scores_all %>%
  filter(method == "BOLD") %>%
  rename(NMDS1_BOLD = NMDS1,
         NMDS2_BOLD = NMDS2) %>% 
  select(-method)

Euclid_dist_all_filt <- scores_all %>%
  filter(method != "BOLD") %>%
  left_join(scores_bold, by = "fieldid") %>% 
  mutate(
    euclid_dist = sqrt(
      (NMDS1 - NMDS1_BOLD)^2 +
        (NMDS2 - NMDS2_BOLD)^2
    )
  )

EUCLID_MEAN_SD <- Euclid_dist_all_filt %>% 
  group_by(method) %>% 
  summarise(mean_Euc = mean(euclid_dist),
            sd_Euc = sd(euclid_dist))

# Figures ####

Fig2A <- ggplot(PHAUS_MBC_ONT_FiltStats_Joint, aes(x = tot_reads, y = a_avg_beta_alpha, colour = replicates)) + 
  geom_line(linewidth = 1.25) + 
  theme_bw() +
  ylab('parameter rank') + 
  xlab('min. OTU read depth') + 
  labs(colour = "Min. OTU # Reps") +
  scale_colour_brewer(palette = "Set2") + 
  scale_y_reverse()
Fig2A


ONT_MAP_MDS_plot <-  plot_mds_comparison(focal = PHAUS_ONT_MAP_Filt)
ILL_MAP_MDS_plot <-  plot_mds_comparison(focal = PHAUS_ILL_MAP_Filt)
ILL_mBRAVE_NMDS_plot <-  plot_mds_comparison(focal = PHAUS_ILL_mBRAVE_Filt)
ILL_SPCFY_NMDS_plot <-  plot_mds_comparison(focal = PHAUS_SPCFY_Filt)

ONT_MAP_Cor_plot <- plot_alpha_slope(PHAUS_ONT_MAP_Filt)
ILL_MAP_Cor_plot <- plot_alpha_slope(PHAUS_ILL_MAP_Filt)
ILL_mBRAVE_Cor_plot <- plot_alpha_slope(PHAUS_ILL_mBRAVE_Filt)
ILL_SPCFY_Cor_plot <- plot_alpha_slope(PHAUS_ILL_mBRAVE_Filt)

PHAUS_MAP_BINs_Filt <- unique(PHAUS_MBC_ONT_filtered$bin_uri)

fit <- euler(list(
  BAR = PHAUS_BOLD_BINs,
  META  = PHAUS_MAP_BINs_Filt
))

# fit <- euler(list(
#   BAR = PHAUS_BOLD_BINs,
#   META  = PHAUS_MAP_BINs_Unfilt
# ))

plot(fit,
     fills = list(fill = c("lightseagreen", "tomato"), alpha = 0.8),
     labels = TRUE,
     quantities = TRUE)

# Total Number of BINs
length(unique(c(PHAUS_BOLD_BINs, PHAUS_MAP_BINs_Filt)))
sum(PHAUS_BOLD_BINs %in% PHAUS_MAP_BINs_Filt)

# write_tsv(PHAUS_MBC_ONT_MAP_Clean, 'data/derived_data/PHAUS_MBC_ONT_MAP_Clean.tsv')
# PHAUS_MBC_ONT_MAP_Clean <- read_tsv('data/derived_data/PHAUS_MBC_ONT_MAP_Clean.tsv')

## Correlation alpha div
PHAUS_BOLD_SpRich <- PHAUS_BOLD_Clean %>% 
  select(fieldid, bin_uri) %>% 
  unique() %>% 
  count(fieldid, name = "BOLD_Rich")

PHAUS_MAP_ONT_SpRich <- PHAUS_MBC_ONT_filtered %>% 
  select(fieldid, bin_uri) %>% 
  unique() %>% 
  count(fieldid, name = "MAP_ONT_Rich")

# Correlation
Sp_Rich_DF_All <- left_join(PHAUS_BOLD_SpRich, PHAUS_MAP_ONT_SpRich)

summary(lm(BOLD_Rich ~ MAP_ONT_Rich, data = Sp_Rich_DF_All))
cor.test(x = Sp_Rich_DF_All$MAP_ONT_Rich, y = Sp_Rich_DF_All$BOLD_Rich)




