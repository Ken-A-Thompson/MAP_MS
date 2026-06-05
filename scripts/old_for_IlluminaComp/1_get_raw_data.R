# MAPLE / RAMS Validation Analysis
# Script author: Ken A. Thompson
# Script 1 of 4: Get Raw Data for PHAUS Analysis
# This script does the following:
## 1.1: Downloads data from 

# load packages
library(BOLDconnectR)
library(readxl)
library(tidyverse)

# BARCODE ####

## Load data
# We load metabarcode data because we use it to filter BOLD data
mbc_output <- read_excel('data/metabarcode_data/Metabarcoding Results - PHAUS_CBGMB01277_COI-5P_658.xlsx')
mbc_samples <- unique(mbc_output$Sample)[1:10]

# Install BOLDconnectR ####
# Note, this code used my PERSONAL API key, which has access to the project.
# That API key has been deleted from this script
# To reproduce the analysis replace 'your-api-key-here' with your API key, and keep the quotes.

bold.apikey('00FE0FAE-94A2-4E94-9224-4CB34BC96A91') # Ken

# Fetch 'res'ults from Project PHAUS
PHAUS_res <- bold.fetch(get_by = "project_codes", 
                        identifiers = "PHAUS") %>% 
  filter(fieldid %in% mbc_samples) 

# then write file
readr::write_tsv(PHAUS_res, 'data/BOLD_downloads/PHAUS_res.tsv')

writexl::write_xlsx(PHAUS_res, 'data/BOLD_downloads/PHAUS_res.xlsx')
mean(is.na(PHAUS_res$nuc)==F)

# Make fasta of PHAUS_res ####
PHAUS_fasta_lines <- PHAUS_res %>%
  dplyr::mutate(
    # Clean sequence: replace anything not A/T/C/G/N with N and convert to uppercase
    clean_nuc = toupper(gsub("[^ATCGN]", "N", nuc)),
    # Make FASTA entry
    fasta_entry = paste0(">", processid, "\n", clean_nuc)
  ) %>%
  dplyr::pull(fasta_entry)

writeLines(PHAUS_fasta_lines, "data/BOLD_downloads/PHAUS_res.fasta")
writeLines(PHAUS_fasta_lines, "vsearch/PHAUS_res.fasta")

# METABARCODE ####

# Make fasta of MBC results ####
## MAPLE ####
mbc_results_MAPLE <- readxl::read_excel('data/metabarcode_data/Metabarcoding Results - PHAUS_CBGMB01277_COI-5P_658.xlsx') %>% 
  mutate(row_number = row_number())


# Which BINs known from countries ####
all_BINs_in_study <- read_tsv('data/derived_data/BINs_det_by_method.tsv')

# Need to loop this;
BIN_geog_results <- bold.fetch(get_by = "bin_uris",
                               identifiers = all_BINs_in_study$BIN[1:10])

BINs_by_country <- BIN_geog_results %>% 
  group_by(bin_uri, country.ocean) %>% 
  summarise(n = n()) %>% 
  na.omit()

write_tsv(BINs_by_country, 'data/derived_data/BINs_by_country.tsv')