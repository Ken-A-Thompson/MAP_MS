# 1.5 barcode summary metrics
# Claude will later clean this up

# Number of BINs total
PHAUS_BOLD_Clean_NTS$bin_uri %>% unique(.) %>% length(.)

# Parent only, NTS only, both
PHAUS_BOLD_Clean_NTS %>% 
  select(record_type, bin_uri) %>% 
  unique() %>% 
  group_by(bin_uri) %>% 
  summarise(
    has_parent = any(record_type == "parent"),
    has_nts    = any(record_type == "nts")
  ) %>% 
  mutate(category = case_when(
    has_parent & has_nts  ~ "both",
    has_parent & !has_nts ~ "only parent",
    !has_parent & has_nts ~ "only NTS"
  )) %>% 
  count(category) %>% 
  mutate(tot = sum(n))

PHAUS_BOLD_Clean_NTS %>% 
  select(bin_uri, class) %>% 
  unique() %>% 
  count(class) %>% 
  mutate(tot = sum(n),
         pct = paste0(round(100*n/tot, 2), '%'))

