library(dplyr)
library(tidyr)
library(vegan)
library(purrr)

#------------------------------
# Function to compute Mantel r
# for a given Reads + Replicates cutoff
#------------------------------

compute_mantel <- function(read_cutoff, rep_cutoff) {
  
  # Build Metabarcode subset for this iteration
  MBC_for_PPM <- mbc_results_MAPLE %>% 
    filter(`BIN Match Status` == "BIN_MATCH") %>%
    filter(Reads >= read_cutoff | Replicates >= rep_cutoff) %>% 
    rename(fieldid = Sample,
           family = Family) %>% 
    mutate(method = "Metabarcode") %>% 
    select(fieldid, family, bin_uri, method)
  
  # Build combined Barcode + Metabarcode dataset
  Data_for_PPM <- PHAUS_res_BOLDistilled %>% 
    select(fieldid, family, bin_uri) %>% 
    mutate(method = "Barcode") %>% 
    bind_rows(., MBC_for_PPM)
  
  # Split by method
  Barcode_df <- Data_for_PPM %>% filter(method == "Barcode")
  Metabarcode_df <- Data_for_PPM %>% filter(method == "Metabarcode")
  
  # Make presence–absence matrices
  barcode_mat <- Barcode_df %>% 
    select(fieldid, bin_uri) %>%
    mutate(value = 1) %>%
    distinct() %>% 
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
  
  # Identify shared samples
  shared_ids <- intersect(barcode_mat$fieldid, metabarcode_mat$fieldid)
  
  # Skip iteration if no shared samples
  if (length(shared_ids) < 3) {
    return(tibble(
      Reads = read_cutoff,
      Replicates = rep_cutoff,
      mantel_r = NA_real_,
      mantel_p = NA_real_,
      n_shared_samples = length(shared_ids)
    ))
  }
  
  # Filter matrices
  barcode_comm <- barcode_mat %>% filter(fieldid %in% shared_ids) %>% arrange(fieldid) %>% select(-fieldid)
  metabarcode_comm <- metabarcode_mat %>% filter(fieldid %in% shared_ids) %>% arrange(fieldid) %>% select(-fieldid)
  
  # Distance matrices
  D_barcode <- vegdist(barcode_comm, method = "jaccard", binary = TRUE)
  D_meta    <- vegdist(metabarcode_comm, method = "jaccard", binary = TRUE)
  
  # Mantel test
  mantel_out <- mantel(D_barcode, D_meta, method = "pearson", permutations = 999)
  
  # Return results
  tibble(
    Reads = read_cutoff,
    Replicates = rep_cutoff,
    mantel_r = as.numeric(mantel_out$statistic),
    mantel_p = mantel_out$signif,
    n_shared_samples = length(shared_ids)
  )
}

#------------------------------
# Run the loop over all combinations
#------------------------------

results_all <- expand_grid(
  Reads = 2:50,
  Replicates = 1:8
) %>%
  pmap_dfr(~ compute_mantel(..1, ..2))

results_all
