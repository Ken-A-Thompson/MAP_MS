beta_surface <- function(
    data,
    ground_truth = PHAUS_BOLD_Clean_NTS,
    sample_col          = "fieldid",
    sp_col              = "bin_uri",
    tr_vals             = 1:40,
    rep_vals            = 1:4,
    mantel_permutations = 999,
    seed                = 142
) {
  
  set.seed(seed)
  
  tr_vals <- tr_vals[tr_vals %in% unique(data$tot_reads)]
  
  ground_truth_BINs_Samples <- ground_truth %>% 
    mutate(method_fieldid = paste0("BOLD_", fieldid)) %>% 
    select(method_fieldid, bin_uri) %>% 
    unique() 
  
  
  out <- vector("list", length(tr_vals) * length(rep_vals))
  idx <- 1L
  
  for (tr in tr_vals) {
    for (rep in rep_vals) {
      
      focal_filt_BINs_Samples <- data %>%
        dplyr::filter(tot_reads >= tr, replicates >= rep) %>%
        mutate(method_fieldid = paste0("MBC_", fieldid)) %>% 
        select(method_fieldid, bin_uri) %>% 
        unique() %>% 
        na.omit()
        
      # --- combine data ---
      all_data <- dplyr::bind_rows(ground_truth_BINs_Samples, focal_filt_BINs_Samples)
      
      all_comm <- make_comm_matrix(all_data, "method_fieldid", sp_col)
      
      BOLD_comm <- all_comm[1:10,]
      MBC_comm <- all_comm[11:20,]
      
      BOLD_dist <- vegan::vegdist(BOLD_comm, method = "bray")
      MBC_dist <- vegan::vegdist(MBC_comm, method = "bray")
      
      mantel_stat <- vegan::mantel(
        BOLD_dist,
        MBC_dist,
        method = "pearson",
        permutations = mantel_permutations
      )$statistic
      
      # Per-sample Bray-Curtis distance between matched GT and MBC samples
      gt_samples  <- sub("^BOLD_", "", rownames(BOLD_comm))
      mbc_samples <- sub("^MBC_",  "", rownames(MBC_comm))
      shared_samp <- intersect(gt_samples, mbc_samples)
      bc_per_samp <- sapply(shared_samp, function(s) {
        pair_mat <- rbind(BOLD_comm[paste0("BOLD_", s), ],
                          MBC_comm[paste0("MBC_",  s), ])
        as.numeric(vegan::vegdist(pair_mat, method = "bray"))
      })

      out[[idx]] <- data.frame(
        tot_reads   = tr,
        replicates  = rep,
        mantel_stat = mantel_stat,
        bc_mean     = mean(bc_per_samp),
        bc_sd       = sd(bc_per_samp)
        # n_samples   = nrow(comm_focal)
      )
      
      idx <- idx + 1L
    }
  }
  
  out <- dplyr::bind_rows(out)
  out_tib <- out %>% 
    mutate(replicates = as_factor(replicates)) %>% 
    tibble(.)
  return(out_tib)
}

alpha_surface <- function(
    data,
    ground_truth = PHAUS_BOLD_Clean_NTS,
    sample_col = "fieldid",
    sp_col     = "bin_uri",
    tr_vals = c(1:40),
    rep_vals = 1:4
) {
  
  tr_vals <- tr_vals[tr_vals %in% unique(data$tot_reads)]
  
  # --- ground truth richness (fixed) ---
  gt_richness <- ground_truth %>%
    dplyr::select(fieldid, bin_uri) %>% 
    unique() %>% 
    dplyr::count(fieldid)
    
  
  results <- data.frame()
  
  for (tr in tr_vals) {
    for (rep in rep_vals) {
      
      # --- filter focal ---
      focal_filtered <- data %>%
        filter(tot_reads >= tr, replicates >= rep)
      
      # --- focal richness ---
      focal_richness <- focal_filtered %>% 
        select(fieldid, bin_uri) %>% 
        unique() %>% 
        dplyr::count(fieldid)
  
      
      # --- join on shared samples ---
      joined <- left_join(gt_richness, focal_richness, by = "fieldid")
      
      # --- correlation ---
      cor_val      <- cor(joined$n.x, joined$n.y)
      spearman_val <- cor(joined$n.x, joined$n.y, method = "spearman")

      lm_obj <- lm(n.y ~ n.x, data = joined)

      # Root mean square error
      y = joined$n.y
      x = joined$n.x

      rmse_1to1 <- sqrt(mean((y - x)^2))
      ccc <- DescTools::CCC(x, y)

      # Also calculate BIN_pct
      PCT_BINs <- mean(PHAUS_BOLD_BINs %in% unique(focal_filtered$bin_uri))

      results <- rbind(results, data.frame(
        tot_reads   = tr,
        replicates  = rep,
        # n_samples   = nrow(joined),
        correlation = cor_val,
        spearman    = spearman_val,
        slope = lm_obj$coefficients[2],
        rmse = rmse_1to1,
        ccc = ccc$rho.c[1],
        PCT_BINs
      ))
    }
  }
  results <- results %>% 
    mutate(replicates = as_factor(replicates)) %>% 
    tibble(.)
  return(results)
}

make_comm_matrix <- function(data, sample_col = "fieldid", sp_col = "bin_uri") {
  data %>%
    dplyr::count(.data[[sample_col]], .data[[sp_col]]) %>%
    tidyr::pivot_wider(
      names_from = all_of(sp_col),
      values_from = n,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames(sample_col) %>%
    as.matrix()
}

plot_alpha_slope <- function(MBC_data) {
  
  BOLD_BIN_Rich <- PHAUS_BOLD_Clean_NTS %>% 
    select(fieldid, bin_uri) %>% 
    unique() %>% 
    dplyr::count(fieldid, name = "BAR")
  
  MBC_BIN_Rich <- MBC_data %>% 
    select(fieldid, bin_uri) %>% 
    unique() %>% 
    dplyr::count(fieldid, name = "META")
  
  BOTH_Rich <- left_join(BOLD_BIN_Rich, MBC_BIN_Rich)
  
  Alpha_Mod <- lm(BAR ~ META, data = BOTH_Rich)
  
  vr <- visreg::visreg(Alpha_Mod, plot = FALSE)
  
  p <- ggplot() +
    geom_line(data = vr$fit,
              aes(x = META, y = visregFit)) +
    geom_ribbon(data = vr$fit,
                aes(x = META, ymin = visregLwr, ymax = visregUpr),
                alpha = 0.2) +
    geom_point(data = vr$res,
               aes(x = META, y = visregRes),
               alpha = 1) +
    theme_bw() + 
    coord_fixed() + 
    geom_abline(intercept = 0, slope = 1, colour = "red")
return(p)
}

plot_mds_comparison <- function(
    ground_truth  = PHAUS_BOLD_Clean_NTS,
    focal,
    sample_col    = "fieldid",
    sp_col        = "bin_uri",
    gt_label      = "GT",
    focal_label   = "Focal",
    point_size    = 2.5,
    alpha         = 0.8,
    k             = 2
) {
  
  ground_truth_BINs_Samples <- ground_truth %>% 
    mutate(method_fieldid = paste0("BOLD_", fieldid)) %>% 
    select(method_fieldid, bin_uri) %>% 
    unique() 
  
  focal_BINs_Samples <- focal %>% 
    mutate(method_fieldid = paste0("MAP_", fieldid)) %>% 
    select(method_fieldid, bin_uri) %>% 
    unique() 
  
  # --- combine data ---
  all_data <- dplyr::bind_rows(ground_truth_BINs_Samples, focal_BINs_Samples)
  
  # --- single community matrix ---
  comm_all <- make_comm_matrix(all_data, sample_col = "method_fieldid", sp_col = "bin_uri")
  
  # --- classical MDS (deterministic PCoA) ---
  dist_all  <- vegan::vegdist(comm_all, method = "bray")
  mds_all   <- cmdscale(dist_all, k = k)                  # base R, no seed needed
  
  # --- scores ---
  scores_all <- as.data.frame(mds_all) %>%
    setNames(paste0("MDS", seq_len(k))) %>%
    tibble::rownames_to_column("sample") %>%
    separate(sample, into = c("method", "fieldid"), sep = "_")
  
  # --- segments (match base sample IDs) ---
  df_segments <- scores_all %>%
    select(method, fieldid, MDS1, MDS2) %>%
    pivot_wider(
      names_from  = method,
      values_from = c(MDS1, MDS2)
    ) %>% 
    transmute(
      fieldid,
      x    = MDS1_BOLD,
      y    = MDS2_BOLD,
      xend = MDS1_MAP,
      yend = MDS2_MAP
    )
  
  # --- plot ---
  p <- ggplot() +
    
    geom_segment(
      data = df_segments,
      aes(x = x, y = y, xend = xend, yend = yend, colour = fieldid),
      linewidth  = 0.4,
      alpha      = 0.45,
      show.legend = FALSE
    ) +
    geom_point(
      data = scores_all,
      aes(MDS1, MDS2, fill = fieldid, shape = method),
      colour = "black",
      size   = point_size,
      alpha  = 0.9,
      stroke = 0.5
    ) +
    
    scale_shape_manual(values = c(21, 22)) +
    scale_fill_viridis_d(option = "turbo") +
    
    guides(
      fill  = "none",
      shape = guide_legend(
        title = "Dataset",
        override.aes = list(size = 3.5, fill = "grey60", colour = "black")
      )
    ) +
    
    labs(x = "MDS1", y = "MDS2", shape = NULL) +
    coord_fixed() + 
    theme_bw(base_size = 12)
  
  return(p)
}

compute_opt_filt <- function(data_alpha, data_beta, data_gamma) {

  joint_data <- left_join(data_alpha, data_beta, by = c("tot_reads", "replicates")) %>%
    left_join(., data_gamma, by = c("tot_reads", "replicates")) %>%
    mutate(ccc_rank  = rank(-est),
           mntl_rank = as.integer(rank(-mantel_stat)),
           pb_rank   = as.integer(rank(-recall)),
           f1_rank   = as.integer(rank(f1))) %>%
    mutate(rank_sum = as.integer(ccc_rank + mntl_rank + f1_rank)) %>%
    select(tot_reads, replicates, est, spearman, mantel_stat, bc_mean, bc_sd,
           recall, f1, ccc_rank, mntl_rank, pb_rank, f1_rank, rank_sum) %>%
    rename(mntl = mantel_stat) %>%
    mutate(across(where(is.double), ~ round(.x, 3)))

}


f1_surface <- function(
    data,
    ground_truth        = PHAUS_BOLD_Clean_NTS,
    sample_col          = "fieldid",
    sp_col              = "bin_uri",
    tr_vals             = 1:40,
    rep_vals            = 1:4,
    seed                = 142
) {
  
  set.seed(seed)
  
  tr_vals <- tr_vals[tr_vals %in% unique(data$tot_reads)]
  
  
  # Ground truth: all BINs for F1; parent-record BINs only for recall (% Target)
  gt_bins     <- ground_truth %>% pull(bin_uri) %>% unique() %>% na.omit()
  parent_bins <- ground_truth %>% filter(record_type == "parent") %>% pull(bin_uri) %>% unique() %>% na.omit()

  out <- vector("list", length(tr_vals) * length(rep_vals))
  idx <- 1L

  for (tr in tr_vals) {
    for (rep in rep_vals) {

      # MBC BINs under this threshold combination
      mbc_bins <- data %>%
        dplyr::filter(tot_reads >= tr, replicates >= rep) %>%
        pull(bin_uri) %>%
        unique() %>%
        na.omit()

      TP <- length(intersect(mbc_bins, gt_bins))
      FP <- length(setdiff(mbc_bins,  gt_bins))
      FN <- length(setdiff(gt_bins,   mbc_bins))

      recall <- mean(parent_bins %in% mbc_bins)
      f1     <- if ((2 * TP + FP + FN) == 0) NA_real_ else (2 * TP) / (2 * TP + FP + FN)

      out[[idx]] <- data.frame(
        tot_reads  = tr,
        replicates = rep,
        TP         = TP,
        FP         = FP,
        FN         = FN,
        precision  = TP / (TP + FP),
        recall     = recall,
        f1         = f1
      )
      
      idx <- idx + 1L
    }
  }
  
  dplyr::bind_rows(out) %>%
    mutate(replicates = as_factor(replicates)) %>%
    tibble()
}



compute_a_b_g <- function(MBC_data, GT_data) {
  
  # 1. Determine filtering
  # 1a Alpha: We use 'CCC'
  PHAUS_MBC_ONT_FiltStats_Alpha <- alpha_surface(data = MBC_data, ground_truth = GT_data) %>%
    select(tot_reads, replicates, est, spearman)
  
  # 1b Beta
  PHAUS_MBC_ONT_FiltStats_Beta <- beta_surface(data = MBC_data, ground_truth = GT_data)
  
  # 1c Gamma
  PHAUS_MBC_ONT_FiltStats_Gamma <- f1_surface(data = MBC_data, ground_truth = GT_data) %>% 
    select(tot_reads, replicates, recall, f1)
  
  # 1c Joint
  PHAUS_MBC_ONT_FiltStats_Joint <- compute_opt_filt(data_alpha = PHAUS_MBC_ONT_FiltStats_Alpha,
                                                    data_beta = PHAUS_MBC_ONT_FiltStats_Beta,
                                                    data_gamma = PHAUS_MBC_ONT_FiltStats_Gamma)
  
}

compute_MinFilt_Metrics <- function(
    MBC_data,
    ground_truth        = PHAUS_BOLD_Clean_NTS,
    threshold_factor    = 1e-05,
    sample_col          = "fieldid",
    sp_col              = "bin_uri",
    mantel_permutations = 999,
    seed                = 142
) {
  set.seed(seed)

  # 1. Proportional threshold filter
  filt_data <- MBC_data %>%
    group_by(.data[[sample_col]]) %>%
    mutate(
      sample_reads = sum(tot_reads),
      threshold    = threshold_factor * sample_reads
    ) %>%
    filter(tot_reads >= threshold) %>%
    ungroup()

  # ── Alpha: CCC of per-sample BIN richness ─────────────────────────────────────
  gt_richness <- ground_truth %>%
    select(all_of(c(sample_col, sp_col))) %>%
    distinct() %>%
    count(.data[[sample_col]])

  focal_richness <- filt_data %>%
    select(all_of(c(sample_col, sp_col))) %>%
    distinct() %>%
    count(.data[[sample_col]])

  joined       <- left_join(gt_richness, focal_richness, by = sample_col)
  ccc_val      <- DescTools::CCC(joined$n.x, joined$n.y)$rho.c[[1]]
  spearman_val <- cor(joined$n.x, joined$n.y, method = "spearman")

  # ── Beta: Mantel R on Bray-Curtis dissimilarity ───────────────────────────────
  gt_long <- ground_truth %>%
    mutate(method_sample = paste0("BOLD_", .data[[sample_col]])) %>%
    select(method_sample, all_of(sp_col)) %>%
    distinct()

  focal_long <- filt_data %>%
    mutate(method_sample = paste0("MBC_", .data[[sample_col]])) %>%
    select(method_sample, all_of(sp_col)) %>%
    distinct() %>%
    na.omit()

  all_comm <- make_comm_matrix(
    bind_rows(gt_long, focal_long),
    sample_col = "method_sample",
    sp_col     = sp_col
  )

  n_gt      <- n_distinct(gt_long$method_sample)
  BOLD_comm <- all_comm[seq_len(n_gt), ]
  MBC_comm  <- all_comm[seq(n_gt + 1L, nrow(all_comm)), ]

  # Restrict to samples present in both matrices (Mantel requires equal dimensions)
  gt_samples  <- sub("^BOLD_", "", rownames(BOLD_comm))
  mbc_samples <- sub("^MBC_",  "", rownames(MBC_comm))
  shared_samp <- intersect(gt_samples, mbc_samples)
  BOLD_comm   <- BOLD_comm[paste0("BOLD_", shared_samp), , drop = FALSE]
  MBC_comm    <- MBC_comm[ paste0("MBC_",  shared_samp), , drop = FALSE]

  mantel_stat <- vegan::mantel(
    vegan::vegdist(BOLD_comm, method = "bray"),
    vegan::vegdist(MBC_comm,  method = "bray"),
    method = "pearson", permutations = mantel_permutations
  )$statistic
  bc_per_samp <- sapply(shared_samp, function(s) {
    pair_mat <- rbind(BOLD_comm[paste0("BOLD_", s), ],
                      MBC_comm[paste0("MBC_",  s), ])
    as.numeric(vegan::vegdist(pair_mat, method = "bray"))
  })
  bc_mean_val <- mean(bc_per_samp)
  bc_sd_val   <- sd(bc_per_samp)

  # ── Gamma: recall (% Target) uses parent-record BINs; F1 uses full BIN set ───
  gt_bins     <- ground_truth %>% pull(all_of(sp_col)) %>% unique() %>% na.omit()
  parent_bins <- ground_truth %>% filter(record_type == "parent") %>% pull(all_of(sp_col)) %>% unique() %>% na.omit()
  focal_bins  <- filt_data    %>% pull(all_of(sp_col)) %>% unique() %>% na.omit()

  recall <- mean(parent_bins %in% focal_bins)

  TP <- length(intersect(focal_bins, gt_bins))
  FP <- length(setdiff(focal_bins,  gt_bins))
  FN <- length(setdiff(gt_bins,     focal_bins))
  f1 <- if ((2 * TP + FP + FN) == 0) NA_real_ else (2 * TP) / (2 * TP + FP + FN)

  tibble(
    threshold_factor = threshold_factor,
    tot_reads        = min(filt_data$tot_reads),
    replicates       = factor(min(filt_data$replicates)),
    n_reads          = sum(filt_data$tot_reads),
    est              = round(ccc_val,      3),
    spearman         = round(spearman_val, 3),
    mntl             = round(mantel_stat,  3),
    bc_mean          = round(bc_mean_val,  3),
    bc_sd            = round(bc_sd_val,    3),
    recall           = round(recall,       3),
    f1               = round(f1,           3)
  )
}