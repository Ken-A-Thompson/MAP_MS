# Paul Hemiptera request

bin_uris <- c(
  "BOLD:ACQ9652",
  "BOLD:AGQ4665",
  "BOLD:AED0367",
  "BOLD:AGS2455",
  "BOLD:AAJ7097",
  "BOLD:ACQ2061"
)

Paul_data <- mbc_results_MAPLE %>% 
  filter(`BIN Hit` %in% bin_uris) %>% 
  select(Sample, Reads, Order, Family, Replicates, `BIN Hit`) %>% 
  rename(bin_uri = `BIN Hit`) %>% 
  group_by(Sample, Order, Family, bin_uri) %>% 
  summarise(total_reads = sum(Reads))

write_xlsx(Paul_data, 'data/PHAUS_Paul_Requests/Hemiptera_BIN_Query.xlsx')

# Sampel IDs

hempitera_sample_ids <- c(
  "CBG-A63272-E08",
  "CBG-A63275-H04",
  "CBG-A63261-B10",
  "CBG-A63261-E08",
  "CBG-A63285-B01",
  "CBG-A63286-D10",
  "CBG-A63305-G09",
  "CBG-A63318-C12",
  "CBG-A63320-F04",
  "CBG-A63325-B12"
)

PHAUS_res_Paul <- PHAUS_res %>% 
  filter(sampleid %in% hempitera_sample_ids) %>% 
  select(processid, order, family, sampleid, bin_uri, fieldid)

write_xlsx(PHAUS_res_Paul, 'data/PHAUS_Paul_Requests/Hempitera_Sample_Query.xlsx')
