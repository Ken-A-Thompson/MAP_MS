# MAPLE & RAMS Validation #
# Author: Ken A. Thompson
# Script 2 of X: Summary of Samples

# Load Data ####
library(tidyverse)

PHAUS_res <- read_csv('data/BOLD_downloads/PHAUS_res.csv')
nrow(PHAUS_res) # n records submitted

# Number of BINs by samples
PHAUS_res_Sample_BINcount <- PHAUS_res %>% 
  drop_na(bin_uri) %>%
  select(fieldid, bin_uri) %>% 
  unique() %>% 
  count(fieldid) %>%
  arrange(n)

# Benchmarking will use GMP#58226; has just over the median number of BINs
  

# Which specimens failed to deliver sequence

PHAUS_res_noSeq <- PHAUS_res %>% 
  filter(is.na(nuc_basecount) == T) 

PHAUS_res_noSeq_Class <- PHAUS_res_noSeq %>% 
  group_by(class) %>% 
  summarise(n = n())
  

# Sequencing runs
table(PHAUS_res$processid_minted_date)[1:13] %>% sum(.) # number in first barcoding run
table(PHAUS_res$processid_minted_date)[14:15] %>% sum(.) # number in second
table(PHAUS_res$processid_minted_date) %>% sum(.) # grand total

# Load PHAUS_res_BIN_Match
PHAUS_res_BIN_Match <- read_tsv('vsearch/BIN_match_results/full_BC/PHAUS_res_BIN_MATCH_RESULTS.tsv') %>% 
  rename(processid = Query) %>% 
  select(processid, `Hit (BIN)`, `BIN Match?`, `%ID`, Kingdom:Species) %>% 
  rename_with(tolower, Kingdom:Species)

PHAUS_res_BOLDistilled <- PHAUS_res %>% 
  drop_na(nuc_basecount) %>% 
  select(processid, fieldid) %>% 
  left_join(., PHAUS_res_BIN_Match, by = "processid") %>% 
  filter(!class %in% c("Alphaproteobacteria", "Gammaproteobacteria")) %>% 
  filter(`BIN Match?` == "BIN_MATCH") %>% 
  rename(bin_uri = `Hit (BIN)`)

# Fraction of samples with sequences ####
sum(is.na(PHAUS_res_BOLDistilled$bin_uri)==F)
sum(is.na(PHAUS_res_BOLDistilled$bin_uri)==F) / nrow(PHAUS_res)

# Breakdown by class ####
PHAUS_res_BOLDistilled %>% 
  group_by(class) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  na.omit() %>% 
  mutate(pct = 100* n / sum(n)) %>% 
  arrange(desc(n))

## Insect order breakdown ####
PHAUS_res %>% 
  drop_na(bin_uri) %>% 
  filter(class == "Insecta") %>% 
  group_by(order) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  na.omit() %>% 
  mutate(pct = 100* n / sum(n)) %>% 
  arrange(desc(pct))

## Collembola order breakdown #####
PHAUS_res_BOLDistilled %>% 
  filter(class == "Collembola") %>% 
  group_by(order) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  na.omit() %>% 
  mutate(pct = 100* n / sum(n)) %>% 
  arrange(desc(pct))

## Arachnid order breakdown #####
PHAUS_res_BOLDistilled %>% 
  filter(class == "Arachnida") %>% 
  group_by(order) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  na.omit() %>% 
  mutate(pct = 100* n / sum(n))%>% 
  arrange(desc(pct))

## How many with sequences LACK BIN ####
PHAUS_res_w_seq <- PHAUS_res %>% 
  drop_na(nuc)
sum(is.na(PHAUS_res_w_seq$bin_uri)==T)
sum(is.na(PHAUS_res_w_seq$bin_uri)==F)

# How many w/ BIN Match? ####
PHAUS_res_BOLDistilled %>% 
  filter(`BIN Match?` == "BIN_MATCH") %>% 
  select(bin_uri) %>% 
  unique() %>% 
  nrow(.)

# Frequencies of singles, etc. ####
PHAUS_res_BOLDistilled %>% 
  filter(`BIN Match?` == "BIN_MATCH") %>% 
  group_by(bin_uri) %>% 
  summarise(n=n()) %>% 
  group_by(n) %>% 
  summarise(n_BINs = n()) %>% 
  mutate(pct = n_BINs / sum(n_BINs) * 100)
  
# Frequencies of BINs ####
PHAUS_res_BOLDistilled %>% 
  filter(`BIN Match?` == "BIN_MATCH") %>% 
  group_by(`Hit (BIN)`) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

# Trap days summary ####
PHAUS_res %>% 
  select(fieldid, collection_date_start, collection_date_end) %>% 
  unique() %>% 
  mutate(collection_date_start = as_date(collection_date_start),
         collection_date_end = as_date(collection_date_end),
         days_between = as.numeric(collection_date_end - collection_date_start))

# Metabarcoding results summary ####
MBC_Animalia <- mbc_results_MAPLE %>% 
  filter(Kingdom == "Animalia") 

# How many ASVs
unique(mbc_results_MAPLE$ASV_ID) %>% length(.)

median(mbc_results_MAPLE$Reads)
median(mbc_results_MAPLE$Replicates)

# How many ASVs matched a BIN?
mbc_results_MAPLE %>%
  select(ASV_ID, `BIN Match Status`) %>% 
  group_by(ASV_ID) %>% 
  mutate(BIN_MATCH = as.integer(any(`BIN Match Status` == "BIN_MATCH"))) %>% 
  filter(BIN_MATCH == max(BIN_MATCH)) %>% 
  ungroup() %>% 
  select(ASV_ID, BIN_MATCH) %>% 
  group_by(BIN_MATCH) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(pct = n / sum(n))
  
# Parameters BIN matches
MAPLE_BIN_Match_only <- filter(mbc_results_MAPLE, `BIN Match Status` == "BIN_MATCH")

median(MAPLE_BIN_Match_only$Reads)
median(MAPLE_BIN_Match_only$Replicates)

 