# figures

library(eulerr) # Venn
library(grid)
library(gridExtra) # double check this one
library(tidyverse)
library(ggalluvial) # for 'flow'
library(patchwork)
library(scales)  # for comma()

my_colors <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # green
  "#0072B2",  # blue
  "#F0E442",  # yellow
  "#D55E00",  # red-orange
  "#CC79A7",  # pink/magenta
  "#999999",  # gray
  "#F39B7F",  # peach
  "#66CC99"   # teal
)
# Load data

mbc_results_MAPLE <- read_excel('data/metabarcode_data/Metabarcoding Results - PHAUS_CBGMB01277_COI-5P_658.xlsx', sheet = 1) %>% 
  mutate(row_number = row_number()) %>% 
  rename(bin_uri = `BIN Hit`)

# Fig 1: Venn / False Pos Neg ####


## 1A: BINs ####
mbc_results_MAPLE_filt <- filter(mbc_results_MAPLE, Reads >= 7 | Replicates >= 2) %>% 
  filter(`BIN Match Status` == "BIN_MATCH") %>% 
  filter(Kingdom == "Animalia")
# Prepare the BIN lists
BC_BINs <- na.omit(unique(PHAUS_res_BOLDistilled$bin_uri))
MBC_BINs <- na.omit(unique(mbc_results_MAPLE_filt$bin_uri))

# How many unique BINs
length(unique(c(BC_BINs, MBC_BINs)))

# Fit the Euler diagram
fit_BIN <- euler(list(BAR = BC_BINs, META = MBC_BINs))
fit_BIN$original.values[3] / sum(fit_BIN$original.values) * 100 # What PCT of BINs were found in both datasets?

# Format the counts with commas
formatted_labels <- format(fit_BIN$original.values, big.mark = "", scientific = FALSE)

# Plot with formatted quantities
venn_plot_BIN <- plot(
  fit_BIN,
  fills = c("skyblue", "orange"),
  quantities = list(labels = formatted_labels),
  labels = TRUE
)
## 1B: BIN/sample combo ####
BC_full_BINs_samples <- PHAUS_res_BOLDistilled %>% 
  select(bin_uri, fieldid) %>% 
  unique()
  
MBC_full_BINs_samples <- mbc_results_MAPLE_filt %>% 
  select(Sample, bin_uri) %>% 
  unique()


BC_BIN_x_sample <- paste(BC_full_BINs_samples$fieldid, BC_full_BINs_samples$bin_uri, sep = "_")
MBC_BIN_x_sample <- paste(MBC_full_BINs_samples$Sample, MBC_full_BINs_samples$bin_uri, sep = "_")

# Fit the Euler diagram
fit_BIN_x_sample <- euler(list(BAR = BC_BIN_x_sample, META = MBC_BIN_x_sample))

fit_BIN_x_sample$original.values[1] + fit_BIN_x_sample$original.values[3]
fit_BIN_x_sample$original.values[2] + fit_BIN_x_sample$original.values[3]

# Format the region counts with commas
formatted_labels <- format(fit_BIN_x_sample$original.values, big.mark = "", scientific = FALSE)

# Plot with formatted quantities
venn_plot_BIN_x_sample <- plot(
  fit_BIN_x_sample,
  fills = c("skyblue", "orange"),
  quantities = list(labels = formatted_labels),
  labels = TRUE
)

# Assemble 1A and 1B
# Wrap each plot with a title
# Create individual plots with left-justified titles
p1 <- arrangeGrob(
  venn_plot_BIN,
  top = textGrob("(A) BINs", x = 0, hjust = 0)
)

p2 <- arrangeGrob(
  venn_plot_BIN_x_sample,
  top = textGrob("(B) BIN-by-sample\ncombinations", x = 0, hjust = 0)
)

library(grid)
library(gridExtra)

# 1️⃣ Draw the combined plot
grid.arrange(p1, p2, ncol = 2)

# 2️⃣ Force grid to expose nested grobs
grid.force()

# 3️⃣ Move the "META" set label (tag.label.2) mostly up, tiny right
grid.edit(
  gPath(
    "diagram.grob.1",
    "tags",
    "tag.number.2",
    "tag.label.2"
  ),
  x = unit(0.92, "npc"),  # tiny nudge to the right
  y = unit(0.15, "npc")    # move mostly up
)

# 4️⃣ Move the corresponding number (tag.quantity.2) up and tiny right
grid.edit(
  gPath(
    "diagram.grob.1",
    "tags",
    "tag.number.2",
    "tag.quantity.2"
  ),
  x = unit(0.95, "npc"),  # same tiny nudge
  y = unit(0.5, "npc")   # slightly below the label
)



# Fig 1C
# Pivot to long format
long_df_Fig2C <- Data_For_Fig2_Plot %>%
  pivot_longer(3:7,
               names_to = "metric",
               values_to = "value")

long_df_Fig2C_ERR <- filter(long_df_Fig2C, metric == "NUM_ERRORS")
long_df_Fig2C_TPOS <- filter(long_df_Fig2C, metric == "NUM_TRUE_POS")
long_df_Fig2C_ACC <- filter(long_df_Fig2C, metric == "ACCURACY_SCORE")
long_df_Fig2C_FPOS <- filter(long_df_Fig2C, metric == "NUM_FALSE_POS")
long_df_Fig2C_NEG <- filter(long_df_Fig2C, metric == "NUM_FALSE_NEG")

# Create heatmap
# Define axis limits based on your data range
x_limits <- c(0, 101)
y_limits <- c(0.5,8.5)

# Left plot: keep axis titles and ticks
Fig2C <- ggplot(long_df_Fig2C_TPOS, aes(x = reads, y = reps, fill = value)) +
  geom_tile(color = "white", linewidth = 0.01) +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = "Min. # reads", y = "Min. # reps", fill = "Count") +
  xlim(x_limits) + ylim(y_limits) +
  theme_minimal() +
  ggtitle('(C) number of true positives') + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(angle = 45, hjust = 1))

# Middle plot: hide y-axis labels
Fig2D <- ggplot(long_df_Fig2C_ERR, aes(x = reads, y = reps, fill = value)) +
  geom_tile(color = "white", linewidth = 0.01) +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = "Min. # reads", y = NULL, fill = "Count") +
  xlim(x_limits) + ylim(y_limits) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggtitle('(D) number of possible errors\n(false positive + false negative)') + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(angle = 45, hjust = 1))

# Right plot: hide y-axis labels
Fig2E <- ggplot(long_df_Fig2C_ACC, aes(x = reads, y = reps, fill = value)) +
  geom_tile(color = "white", linewidth = 0.01) +
  # scale_fill_viridis_b(option = "plasma") +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = "Min. # reads", y = NULL, fill = "Count") +
  xlim(x_limits) + ylim(y_limits) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggtitle('(E) dataset similarity\n(true positives \u2212 errors)') + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(angle = 45, hjust = 1))

# Combine plots with patchwork
Fig2_CtoE <- Fig2C | Fig2D | Fig2E
ggsave(Fig2_CtoE, filename = 'figures/main/Fig2_C2E.png', width = 220, height = 100, units = "mm")
# copy at 1000x433

## Fig S_E: False Pos/Neg ####
Fig_S_Error_Pos <- ggplot(long_df_Fig2C_FPOS, aes(x = reads, y = reps, fill = value)) +
  geom_tile(color = "white", linewidth = 0.01) +
  scale_fill_viridis_c(option = "plasma") + 
  labs(x = "Min. # reads", y = "Min. # reps", fill = "Count") +
  theme_minimal() + 
  ggtitle('(A) # false positives')


Fig_S_Error_Neg <- ggplot(long_df_Fig2C_NEG, aes(x = reads, y = reps, fill = value)) +
  geom_tile(color = "white", linewidth = 0.01) +
  scale_fill_viridis_c(option = "plasma") + 
  labs(x = "Min. # reads", y = "Min. # reps", fill = "Count") +
  theme_minimal() + 
  ggtitle('(B) # false negatives')

Fig_SX_Error <- Fig_S_Error_Pos + Fig_S_Error_Neg
ggsave(Fig_SX_Error, filename='figures/supp/Fig_SX_FalsePos_FalseNeg.png', width = 160, height = 80, units = "mm")
# Copy at 800x333

# Fig 3: Taxonomy ####
Fig_3_Tax <- ggplot(stacked_df, aes(x = order, y = n, fill = detection_status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      "metabarcode" = "#377EB8",     # blue
      "barcode" = "#E41A1C",         # red
      "detected_by_both" = "#4DAF4A" # green
    ),
    labels = c(
      "metabarcode" = "Missed by Metabarcode",
      "barcode" = "Missed by Barcode",
      "detected_by_both" = "Detected by Both"
    )
  ) +
  labs(
    x = "Order",
    y = "Number of BINs",
    fill = "Detection Status",
    title = "BIN Detection by Method and Order"
  ) +
  facet_wrap(~ commonness, scales = "free") +
  theme_minimal() +
  coord_flip()

ggsave(Fig_3_Tax, filename = 'figures/main/Fig3.png', width = 200, height = 80, units = "mm")

# Fig 4: BOLDistilled Flow ####

sequential_dat_forplot <- sequential_dat %>% 
  pivot_longer(1:2, names_to = "groups", values_to = "value") %>% 
  group_by(groups) %>%
  mutate(id = cur_group_id()) %>%
  ungroup() %>%
  mutate(
    BD_lib = fct_recode(BD_lib, "no site" = "no park"),
    BD_lib = factor(BD_lib, levels = c("no Aus", "no NSW", "no site", "complete"))
  ) %>% 
  mutate(groups_2 = ifelse(groups == "BIN_MATCH", "BIN match", "no BIN match"))

Flow_Fig_All_Libs <- ggplot(sequential_dat_forplot, aes(x = BD_lib, stratum = forcats::fct_rev(groups_2), alluvium = id, y = value, fill = groups_2)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", alpha = 0.5) +
  geom_stratum() +
  theme_minimal() +
  scale_fill_manual(values = c("#2ca02c", "black")) + 
  ylab('# ASVs') + 
  xlab('library resolution') + 
  theme(legend.title = element_blank())

ggsave(plot = Flow_Fig_All_Libs, filename = 'figures/main/Fig5.png', height = 90, width = 150, units = "mm")
# clip 650 x 350

# Fig 3: Dark BIN ASVs ####
# Plot with correct stacking and colors
Fig_3_ggplot <- ggplot(DarkBIN_plot_data, aes(x = factor(n_reps), y = count, fill = type)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      dark_bin = "lightseagreen",       # green = Dark BIN
      non_dark_bin = "tomato"    # black = Non-Dark BIN
    ),
    labels = c(
      non_dark_bin = "Lacking BIN",
      dark_bin = "Dark BIN"
    )
  ) +
  labs(
    x = "Minimum # Replicates",
    y = "# unidentified ASVs",
    fill = "ASV Type",
  ) +
  theme_minimal()

# Fig 4: Alpha/Beta ####
alpha_plot <- 
  ggplot(Alpha_Comp_DF, aes(x = BIN_Count_BC, BIN_Count_MBC)) + 
  geom_point(aes(colour = duration)) + 
  geom_smooth(method = "lm") + 
  coord_fixed() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  scale_colour_manual(values = c("tomato", "lightseagreen")) + 
  xlim(c(0, 600)) +
  ylim(c(0, 575)) + 
  ylab('BIN count (metabarcode)') + 
  xlab('BIN count (barcode)') + 
  ggtitle("(A) α-diversity") + 
  theme(
    legend.position = c(0.95, 0.05),      # lower right, inside plot
    legend.justification = c(1, 0),
    legend.background = element_blank(),  # no legend box background
    legend.key = element_blank(),         # no key background  
    legend.title = element_blank()
    )

beta_plot <- plot_PA_MDS(data = Data_for_PPM, plot_title = NA, group_col = "bin_uri") + 
  ggtitle('(B) β-diversity') + 
  xlab('MDS axis 1') + 
  ylab('MDS axis 2')

Fig4_Alpha_Beta <- alpha_plot + beta_plot

ggsave(Fig4_Alpha_Beta, filename = 'figures/main/Fig4_A_B_Div.png', width = 180, height = 100, units = "mm")

# 750 x 400

# SUPPLEMENTARY FIGS ####
## S1 Specimen Bar Chart of Class/Order ####
class_counts <- PHAUS_res_BOLDistilled %>%
  filter(!is.na(class)) %>%
  count(class)

order_counts <- PHAUS_res_BOLDistilled %>%
  filter(!is.na(order)) %>%
  count(order)

FigS1A_Class <- ggplot(class_counts, aes(x = fct_reorder(class, n, .desc = TRUE), y = n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = comma(n)), vjust = -0.5) +
  theme_bw() +
  labs(x = "Class", y = "Specimen Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  expand_limits(y = max(class_counts$n) * 1.1)

FigS1B_Order <- ggplot(order_counts, aes(x = fct_reorder(order, n), y = n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = comma(n)), hjust = -0.1) +
  theme_bw() +
  labs(x = "Order", y = "Specimen Count") +
  coord_flip(clip = "off") +
  expand_limits(y = max(order_counts$n) * 1.2)

FigS1 <- FigS1A_Class + FigS1B_Order + 
  plot_layout(widths = c(1, 2)) + 
  plot_annotation(tag_levels = 'A')

ggsave(plot = FigS1, filename = 'figures/supp/Fig_S1_Specimens_Class_Order.png', width = 200, height = 90, units = "mm")

## S2 BINs Bar Chart of Class/Order ####
class_counts_BINs <- PHAUS_res_BOLDistilled %>%
  select(bin_uri, class) %>% 
  unique() %>% 
  filter(!is.na(class)) %>%
  count(class)

order_counts_BINs <- PHAUS_res_BOLDistilled %>%
  select(bin_uri, order) %>% 
  unique() %>% 
  filter(!is.na(order)) %>%
  count(order)

FigS2_BINs_Class <- ggplot(class_counts_BINs, aes(x = fct_reorder(class, n, .desc = TRUE), y = n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = comma(n)), vjust = -0.5) +
  theme_bw() +
  labs(x = "Class", y = "BIN Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  expand_limits(y = max(class_counts_BINs$n) * 1.1)

FigS2_BINs_Order <- ggplot(order_counts_BINs, aes(x = fct_reorder(order, n), y = n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = comma(n)), hjust = -0.1) +
  theme_bw() +
  labs(x = "Order", y = "BIN Count") +
  coord_flip(clip = "off") +
  expand_limits(y = max(order_counts_BINs$n) * 1.2)

FigS2 <- FigS2_BINs_Class + FigS2_BINs_Order + 
  plot_layout(widths = c(1, 2)) + 
  plot_annotation(tag_levels = 'A')

ggsave(plot = FigS2, filename = 'figures/supp/Fig_S2_BINs_Class_Order.png', width = 200, height = 90, units = "mm")

## S3: BARCODING BIN Abundance ####
BC_full_BIN_freq <- PHAUS_res_BOLDistilled %>% 
  group_by(bin_uri) %>% 
  summarise(n = n()) %>% 
  na.omit()

Most_Freq_BINs <- filter(BC_full_BIN_freq, n > 200)

Most_Freq_BINs_Tax <- filter(PHAUS_res_BOLDistilled, bin_uri %in% Most_Freq_BINs$bin_uri) %>% 
  select(bin_uri, phylum:species) %>% 
  select(-subfamily) %>% 
  unique() %>% 
  left_join(., Most_Freq_BINs, by = "bin_uri") %>% 
  arrange(desc(n))

Fig_S3 <- ggplot(BC_full_BIN_freq, aes(x = n)) +
  geom_histogram(bins = 50) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), 
                     breaks = c(1, 10, 100, 1000),
                     labels = c("1", "10", "100", "1,000")) +
  labs(y = "Number of BINs (log10-scaled)") +
  theme_minimal() + 
  geom_hline(yintercept = 0) + 
  xlab('number of observations')

ggsave(plot = Fig_S3, filename = 'figures/supp/Fig_S3_BIN_Freq.png', width = 110, height = 90, units = "mm")


## S4: Read Depth & Replicates ####
median(mbc_results_MAPLE$Replicates)
median(mbc_results_MAPLE$Reads)

Fig_S4A <- ggplot(mbc_results_MAPLE, aes(x = Replicates)) + 
  geom_histogram() + 
  geom_histogram(binwidth = 0.5, boundary = 0.3) +
  scale_x_continuous(breaks = 1:8) +
  scale_y_log10() + 
  labs(x = "# technical replicates", y = "ASV count") +
  theme_minimal() 

Fig_S4B <- ggplot(mbc_results_MAPLE, aes(x = Reads)) + 
  geom_histogram(bins = 40) + 
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = comma_format(accuracy = 1)
  ) +
  labs(x = "# reads", y = "ASV count") +
  scale_y_log10() + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 

FigS4 <- Fig_S4A + Fig_S4B + 
  plot_annotation(tag_levels = 'A') 

ggsave(plot = FigS4, filename = 'figures/supp/Fig_S4_MBC_Tech_All_ASVs.png', width = 160, height = 90, units = "mm")

# Fig SX: Det Prob ####
# Make n_barcoded_bin an ordered factor
det_prob_binned_plot <- det_prob_binned %>%
  mutate(n_barcoded_bin = factor(n_barcoded_bin, 
                                 levels = c("1", "2", "3–6", "7–10", "11–20", "21–50", "51–100", ">100"),
                                 ordered = TRUE))

# Plot
Fig_S_Det_Prob <- ggplot(det_prob_binned_plot, aes(x = n_barcoded_bin, y = frac_true)) + 
  geom_col() +
  xlab("Number of barcoded specimens") +
  ylab("Fraction detected via metabarcoding") +
  geom_text(aes(label = n), vjust = -0.5) +  # numbers above bars
  ylim(0, 1.05) + 
  theme_minimal()

Fig_S_Det_Prob

ggsave(Fig_S_Det_Prob, filename = 'figures/supp/Fig_S_P_Det.png', width = 120, height = 90, units = "mm")

# Fig S5 ####
mbc_results_MAPLE_BINsONLY <- filter(mbc_results_MAPLE, `BIN Match Status` == "BIN_MATCH") %>% 
  group_by(Sample, `BIN Hit`) %>% 
  slice_max(order_by = Reads, n = 1, with_ties = FALSE) %>%
  ungroup()

Fig_S5A <- ggplot(mbc_results_MAPLE_BINsONLY, aes(x = Replicates)) + 
  geom_histogram() + 
  geom_histogram(binwidth = 0.5, boundary = 0.3) +
  scale_x_continuous(breaks = 1:8) +
  labs(x = "# technical replicates", y = "ASV count") +
  theme_minimal()


Fig_S5B <- ggplot(mbc_results_MAPLE_BINsONLY, aes(x = Reads)) + 
  geom_histogram(bins = 40) + 
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = comma_format(accuracy = 1)
  ) +
  labs(x = "# reads", y = "ASV count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

FigS5 <- Fig_S5A + Fig_S5B + 
  plot_annotation(tag_levels = 'A') 

ggsave(plot = FigS5, filename = 'figures/supp/Fig_S5_MBC_Tech_BINs.png', width = 160, height = 90, units = "mm")

# copy to clip; 725x275

# S Map ####

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Get Australia map
aus_map <- ne_countries(country = "Australia", scale = "medium", returnclass = "sf")

# Define bounding boxes
aus_bbox <- data.frame(
  lon = c(112, 154, 154, 112, 112),
  lat = c(-44, -44, -10, -10, -44)
)

nsw_bbox <- data.frame(
  lon = c(141, 154, 154, 141, 141),
  lat = c(-37.5, -37.5, -28, -28, -37.5)
)


muog_bbox <- data.frame(lon = c(151.056393, 151.272193, 151.272193, 151.056393, 151.056393),
                        lat = c(-33.688903, -33.688903, -33.509503, -33.509503, -33.688903))

# Convert to sf polygons
aus_poly <- st_polygon(list(as.matrix(aus_bbox))) %>% st_sfc(crs = st_crs(aus_map))
nsw_poly <- st_polygon(list(as.matrix(nsw_bbox))) %>% st_sfc(crs = st_crs(aus_map))
muo_poly <- st_polygon(list(as.matrix(muog_bbox))) %>% st_sfc(crs = st_crs(aus_map))

muo_bbox <- st_bbox(muo_poly)

# Top-left corner of the hotpink square
muo_top_left_x <- as.numeric(muo_bbox["xmin"])
muo_top_left_y <- as.numeric(muo_bbox["ymax"])

arrow_start_x <- as.numeric(
  nsw_bbox["xmin"] + 0.6 * (nsw_bbox["xmax"] - nsw_bbox["xmin"])
)
arrow_start_y <- as.numeric(
  nsw_bbox["ymin"] + 0.6 * (nsw_bbox["ymax"] - nsw_bbox["ymin"])
)

# Add a small offset (in degrees, adjust as needed)
offset_x <- -0.3   # moves the arrow end slightly to the right
offset_y <- 0.3  # moves the arrow end slightly down

arrow_end_x <- muo_top_left_x + offset_x
arrow_end_y <- muo_top_left_y + offset_y


# Plot
ggplot() +
  geom_sf(data = aus_map, fill = "lightyellow", color = "black") +
  geom_sf(data = aus_poly, fill = NA, color = "blue", linewidth = 1.2) +
  geom_sf(data = nsw_poly, fill = NA, color = "red", linewidth = 1.2) +
  geom_sf(data = muo_poly, fill = "hotpink", color = "hotpink") +
  geom_segment(
    aes(
      x = arrow_start_x,
      y = arrow_start_y,
      xend = arrow_end_x,
      yend = arrow_end_y
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 1,
    color = "black"
  ) + 
  coord_sf(xlim = c(110, 160), ylim = c(-45, -5)) +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 14)
  )

# Copy to clip at 800 x 300



# Revision figures
