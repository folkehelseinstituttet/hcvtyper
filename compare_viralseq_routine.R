library(tidyverse)

if (Sys.info()['sysname'] == "Linux") {
  prefix <- "/mnt/N/"
} else if (Sys.info()['sysname'] == "Windows") {
  prefix <- "N:/"
}


# Default Bowtie2 ---------------------------------------------------------



# Read viralseq output
viralseq_Run879 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run879/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>%
  # Add run name
  add_column("Run" = "Run879") %>%
  # Remove duplicated entries
  distinct()
viralseq_Run882 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run882_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run882") %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run891 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run891_MiSeq_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run891") %>% 
  # Remove sample that has be re-run in Run897 with more mapped reads
  filter(Sample != "HCV082023A1") %>% 
  filter(Sample != "HCV082023F3") %>% 
  filter(Sample != "HCV082023A3") %>% 
  filter(Sample != "HCV082023B1") %>% 
  filter(Sample != "HCV082023B2") %>% 
  filter(Sample != "HCV082023D2") %>% 
  filter(Sample != "HCV082023D3") %>% 
  filter(Sample != "HCV082023E1") %>% 
  filter(Sample != "HCV082023E3") %>% 
  filter(Sample != "HCV082023G2") %>% 
  filter(Sample != "HCV082023H1") %>% 
  # Remove duplicated samples
  distinct()
viralseq_Run897 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run897_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run897") %>% 
  # Remove sample that has be run in Run891 with more mapped reads
  filter(Sample != "HCV082023C3") %>% 
  filter(Sample != "HCV082023F1") %>% 
  filter(Sample != "HCV082023H2") %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run898 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run898_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run898") %>% 
  mutate(glecaprevir_mut_short = as.character(glecaprevir_mut_short)) %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  mutate(NS34A_short = as.character(NS34A_short)) %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run902 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run902_virus_og_bakt/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run902") %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run911 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run911_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run911") %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run917 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run917_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run917") %>% 
  # Remove duplicated entries
  distinct()
viralseq_20240416 <- read_csv(paste0(prefix, "Virologi/JonBrate/NGS_SEQ-20240416-01/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "20240416-01") %>% 
  mutate(glecaprevir_mut_short = as.character(glecaprevir_mut_short)) %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  mutate(voxilaprevir_mut_short = as.character(voxilaprevir_mut_short)) %>% 
  # Some samples are duplicated with identical entries in all columns
  distinct()
viralseq_20240617 <- read_csv(paste0(prefix, "Virologi/JonBrate/NGS_SEQ-20240617-01/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "20240627-01") %>% 
  # Some samples are duplicated with identical entries in all columns
  distinct()
viralseq_20240814 <- read_csv(paste0(prefix, "Virologi/JonBrate/NGS_SEQ-20240814-01/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "20240814-01") %>% 
  mutate(glecaprevir_mut_short = as.character(glecaprevir_mut_short)) %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  mutate(voxilaprevir_mut_short = as.character(voxilaprevir_mut_short)) %>% 
  # Some samples are duplicated with identical entries in all columns
  distinct()

viralseq <- bind_rows(viralseq_Run879, viralseq_Run882, viralseq_Run891, viralseq_Run897, viralseq_Run898, viralseq_Run902, viralseq_Run911, viralseq_Run917, viralseq_20240416, viralseq_20240617, viralseq_20240814)
vs_cols <- colnames(viralseq)

# Read routine results
routine_Run879 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run879_HCV_summaries/Run879_HCV_summary_with_glue.tsv"),
                             delim = "\t", escape_double = FALSE,
                             col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                              `Most abundant minority genotype:` = col_character(),
                                              `Percent mapped reads of trimmed:` = col_character(),
                                              `Percent covered:` = col_character(),
                                              `Percent covered above depth=5 without duplicates:` = col_character(),
                                              `Percent covered above depth=9 without duplicates:` = col_character(),
                                              `Percent covered minor:` = col_character(),
                                              `Percent covered above depth=5 minor without duplicates:` = col_character(),
                                              `Percent covered above depth=9 minor without duplicates:` = col_character(),
                                              glecaprevir_mut_short = col_character(),
                                              paritaprevir_mut_short = col_character(),
                                              voxilaprevir_mut_short = col_character(),
                                              NS34A_short = col_character()),
                             trim_ws = TRUE) %>%
  # Convert comma to "." and then to double
  mutate(`Average depth without duplicates:` = str_replace(`Average depth without duplicates:`, ",", "\\.")) %>%
  mutate(`Average depth without duplicates:` = as.double(`Average depth without duplicates:`)) %>%
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>%
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>%
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>%
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>%
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>%
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>%
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>%
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>%
  # Change "N/A" to NA and convert to double
  mutate(`Number of mapped reads minor:` = na_if(`Number of mapped reads minor:`, "N/A")) %>%
  mutate(`Number of mapped reads minor:` = as.numeric(`Number of mapped reads minor:`)) %>%
  mutate(`Number of mapped reads minor without duplicates:` = na_if(`Number of mapped reads minor without duplicates:`, "N/A")) %>%
  mutate(`Number of mapped reads minor without duplicates:` = as.numeric(`Number of mapped reads minor without duplicates:`)) %>%
  # Convert some columns to double
  mutate(`Number of mapped reads:` = as.numeric(`Number of mapped reads:`)) %>%
  mutate(`Number of mapped reads without duplicates:` = as.numeric(`Number of mapped reads without duplicates:`)) %>%
  # Add run name
  add_column("Run" = "Run879") %>%
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))
routine_Run882 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run882_HCV_summaries/Run882_Virus_summary_with_glue.tsv"), 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                              `Most abundant minority genotype:` = col_character(),
                                              `Percent mapped reads of trimmed:` = col_character(), 
                                              `Percent covered:` = col_character(), 
                                              `Percent covered above depth=5 without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 without duplicates:` = col_character(), 
                                              `Percent covered minor:` = col_character(), 
                                              `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                              glecaprevir_mut_short = col_character(),
                                              paritaprevir_mut_short = col_character(),
                                              voxilaprevir_mut_short = col_character(),
                                              NS34A_short = col_character()),
                             trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Average depth without duplicates:` = str_replace(`Average depth without duplicates:`, ",", "\\.")) %>% 
  mutate(`Average depth without duplicates:` = as.double(`Average depth without duplicates:`)) %>% 
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Change "N/A" to NA and convert to double
  mutate(`Number of mapped reads minor:` = na_if(`Number of mapped reads minor:`, "N/A")) %>% 
  mutate(`Number of mapped reads minor:` = as.numeric(`Number of mapped reads minor:`)) %>% 
  mutate(`Number of mapped reads minor without duplicates:` = na_if(`Number of mapped reads minor without duplicates:`, "N/A")) %>% 
  mutate(`Number of mapped reads minor without duplicates:` = as.numeric(`Number of mapped reads minor without duplicates:`)) %>% 
  # Convert some columns to double
  mutate(`Number of mapped reads:` = as.numeric(`Number of mapped reads:`)) %>% 
  mutate(`Number of mapped reads without duplicates:` = as.numeric(`Number of mapped reads without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "Run882") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))
routine_Run891 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run891_HCV_summaries/Run891_HCV_summary_with_glue.tsv"), 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                              `Most abundant minority genotype:` = col_character(),
                                              `Percent mapped reads of trimmed:` = col_character(), 
                                              `Percent covered:` = col_character(), 
                                              `Percent covered above depth=5 without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 without duplicates:` = col_character(), 
                                              `Percent covered minor:` = col_character(), 
                                              `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                              glecaprevir_mut_short = col_character(),
                                              paritaprevir_mut_short = col_character(),
                                              voxilaprevir_mut_short = col_character(),
                                              NS34A_short = col_character()),
                             trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Average depth without duplicates:` = str_replace(`Average depth without duplicates:`, ",", "\\.")) %>% 
  mutate(`Average depth without duplicates:` = as.double(`Average depth without duplicates:`)) %>% 
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Change "N/A" to NA and convert to double
  mutate(`Number of mapped reads minor:` = na_if(`Number of mapped reads minor:`, "N/A")) %>% 
  mutate(`Number of mapped reads minor:` = as.numeric(`Number of mapped reads minor:`)) %>% 
  mutate(`Number of mapped reads minor without duplicates:` = na_if(`Number of mapped reads minor without duplicates:`, "N/A")) %>% 
  mutate(`Number of mapped reads minor without duplicates:` = as.numeric(`Number of mapped reads minor without duplicates:`)) %>% 
  # Convert some columns to double
  mutate(`Number of mapped reads:` = as.numeric(`Number of mapped reads:`)) %>% 
  mutate(`Number of mapped reads without duplicates:` = as.numeric(`Number of mapped reads without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "Run891") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  )) %>% 
  # Remove negative sample with weird entries
  filter(Sample != "HCV082023negH4") %>% 
  # Remove some samples that have been run also in Run897 with more mapped reads
  filter(Sample != "HCV082023A1") %>% 
  filter(Sample != "HCV082023F3") %>% 
  filter(Sample != "HCV082023A3") %>% 
  filter(Sample != "HCV082023B1") %>% 
  filter(Sample != "HCV082023B2") %>% 
  filter(Sample != "HCV082023D2") %>% 
  filter(Sample != "HCV082023D3") %>% 
  filter(Sample != "HCV082023E1") %>% 
  filter(Sample != "HCV082023E3") %>% 
  filter(Sample != "HCV082023G2") %>% 
  filter(Sample != "HCV082023H1")
routine_Run897 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run897_HCV_summaries/Run897_HCV_summary_with_glue.tsv"), 
                              delim = "\t", escape_double = FALSE, 
                              col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                               `Most abundant minority genotype:` = col_character(),
                                               `Percent mapped reads of trimmed:` = col_character(), 
                                               `Percent covered:` = col_character(), 
                                               `Percent covered above depth=5 without duplicates:` = col_character(), 
                                               `Percent covered above depth=9 without duplicates:` = col_character(), 
                                               `Percent covered minor:` = col_character(), 
                                               `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                               `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                               glecaprevir_mut_short = col_character(),
                                               paritaprevir_mut_short = col_character(),
                                               voxilaprevir_mut_short = col_character(),
                                               NS34A_short = col_character()),
                              trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Average depth without duplicates:` = str_replace(`Average depth without duplicates:`, ",", "\\.")) %>% 
  mutate(`Average depth without duplicates:` = as.double(`Average depth without duplicates:`)) %>% 
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Change "N/A" to NA and convert to double
  mutate(`Number of mapped reads minor:` = na_if(`Number of mapped reads minor:`, "N/A")) %>% 
  mutate(`Number of mapped reads minor:` = as.numeric(`Number of mapped reads minor:`)) %>% 
  mutate(`Number of mapped reads minor without duplicates:` = na_if(`Number of mapped reads minor without duplicates:`, "N/A")) %>% 
  mutate(`Number of mapped reads minor without duplicates:` = as.numeric(`Number of mapped reads minor without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "Run897") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  )) %>% 
  # Remove sample HCV082023H2. Has been run also in Run891, with more mapped reads and an additional RAS in NS34A. Run891 is similar to Bowtie2
  filter(Sample != "HCV082023H2") %>% 
  # Remove some samples that have been run also in Run891 with more mapped reads
  filter(Sample != "HCV082023C3") %>% 
  filter(Sample != "HCV082023F1")
routine_Run898 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run898_HCV_summaries/Run898_HCV_summary_with_glue.tsv"), 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                              `Most abundant minority genotype:` = col_character(),
                                              `Percent mapped reads of trimmed:` = col_character(), 
                                              `Percent covered:` = col_character(), 
                                              `Percent covered above depth=5 without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 without duplicates:` = col_character(), 
                                              `Percent covered minor:` = col_character(), 
                                              `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                              glecaprevir_mut_short = col_character(),
                                              paritaprevir_mut_short = col_character(),
                                              voxilaprevir_mut_short = col_character(),
                             NS34A_short = col_character()),
                             trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Average depth without duplicates:` = str_replace(`Average depth without duplicates:`, ",", "\\.")) %>% 
  mutate(`Average depth without duplicates:` = as.double(`Average depth without duplicates:`)) %>% 
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "Run898") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))
routine_Run902 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run902_HCV_summaries/Run902_HCV_summary_with_glue.tsv"), 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                              `Percent mapped reads of trimmed:` = col_character(), 
                                              `Percent covered:` = col_character(), 
                                              `Percent covered above depth=5 without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 without duplicates:` = col_character(), 
                                              `Percent covered minor:` = col_character(), 
                                              `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                              glecaprevir_mut_short = col_character(),
                                              paritaprevir_mut_short = col_character(),
                                              voxilaprevir_mut_short = col_character()), 
                             trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Average depth without duplicates:` = str_replace(`Average depth without duplicates:`, ",", "\\.")) %>% 
  mutate(`Average depth without duplicates:` = as.double(`Average depth without duplicates:`)) %>% 
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "Run902") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))
routine_Run911 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run911_HCV_summaries/Run911_HCV_summary_with_glue.tsv"), 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                              `Percent mapped reads of trimmed:` = col_character(), 
                                              `Percent covered:` = col_character(), 
                                              `Percent covered above depth=5 without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 without duplicates:` = col_character(), 
                                              `Percent covered minor:` = col_character(), 
                                              `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                              `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                              glecaprevir_mut_short = col_character(),
                                              paritaprevir_mut_short = col_character(),
                                              voxilaprevir_mut_short = col_character()), 
                             trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "Run911") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))
routine_Run917 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run917_HCV_summaries/Run917_HCV_summary_with_glue.tsv"), 
                               delim = "\t", escape_double = FALSE, 
                               col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                                `Percent mapped reads of trimmed:` = col_character(), 
                                                `Percent covered:` = col_character(), 
                                                `Percent covered above depth=5 without duplicates:` = col_character(), 
                                                `Percent covered above depth=9 without duplicates:` = col_character(), 
                                                `Percent covered minor:` = col_character(), 
                                                `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                                `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                                glecaprevir_mut_short = col_character(),
                                                paritaprevir_mut_short = col_character(),
                                                voxilaprevir_mut_short = col_character()), 
                               trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "Run917") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))
routine_20240416 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2024/NGS_SEQ-20240416-01_summaries/NGS_SEQ-20240416-01_summary_with_glue.tsv"), 
                               delim = "\t", escape_double = FALSE, 
                               col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                                `Percent mapped reads of trimmed:` = col_character(), 
                                                `Percent covered:` = col_character(), 
                                                `Percent covered above depth=5 without duplicates:` = col_character(), 
                                                `Percent covered above depth=9 without duplicates:` = col_character(), 
                                                `Percent covered minor:` = col_character(), 
                                                `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                                `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                                glecaprevir_mut_short = col_character(),
                                                paritaprevir_mut_short = col_character(),
                                                voxilaprevir_mut_short = col_character()), 
                               trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "20240416-01") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))
routine_20240617 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2024/NGS_SEQ-20240617-01_HCV_summaries/NGS_SEQ-20240617-01_HCV_summary_with_glue.tsv"), 
                               delim = "\t", escape_double = FALSE, 
                               col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                                `Percent mapped reads of trimmed:` = col_character(), 
                                                `Percent covered:` = col_character(), 
                                                `Percent covered above depth=5 without duplicates:` = col_character(), 
                                                `Percent covered above depth=9 without duplicates:` = col_character(), 
                                                `Percent covered minor:` = col_character(), 
                                                `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                                `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                                glecaprevir_mut_short = col_character(),
                                                paritaprevir_mut_short = col_character(),
                                                voxilaprevir_mut_short = col_character()), 
                               trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "20240627-01") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))
routine_20240814 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2024/NGS_SEQ-20240814-01_HCV_summaries/NGS_SEQ-20240814-01_HCV_summary_with_glue.tsv"), 
                               delim = "\t", escape_double = FALSE, 
                               col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                                `Percent mapped reads of trimmed:` = col_character(), 
                                                `Percent covered:` = col_character(), 
                                                `Percent covered above depth=5 without duplicates:` = col_character(), 
                                                `Percent covered above depth=9 without duplicates:` = col_character(), 
                                                `Percent covered minor:` = col_character(), 
                                                `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                                `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                                glecaprevir_mut_short = col_character(),
                                                paritaprevir_mut_short = col_character(),
                                                voxilaprevir_mut_short = col_character()), 
                               trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "20240814-01") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))

routine <- bind_rows(routine_Run879, routine_Run882, routine_Run891, routine_Run897, routine_Run898, routine_Run902, routine_Run911, routine_Run917, routine_20240416, routine_20240617)
rt_cols <- colnames(routine)

# Compare columns
setdiff(vs_cols, rt_cols)

# Join data - keep only samples present in both viralseq and routine
joined <- inner_join(viralseq, routine, by = "Sample") 

# Compare
comparison <- joined %>% 
  mutate("same_major_geno" = case_when(
    `Majority genotype:.x` == `Majority genotype:.y` ~ "YES",
    .default = "NO"
  )) %>% 
  mutate("diff_mapped_reads_with_dup" = `Number of mapped reads:.x` - `Number of mapped reads:.y`) %>% 
  mutate("perc_diff_mapped_reads_with_dup" = `Number of mapped reads:.y` / `Number of mapped reads:.x` * 100) %>% 
  mutate("diff_mapped_reads_no_dup" = `Number of mapped reads minor without duplicates:.x` - `Number of mapped reads minor without duplicates:.y`) %>% 
  mutate("perc_diff_mapped_reads_no_dup" = `Number of mapped reads minor without duplicates:.y` / `Number of mapped reads minor without duplicates:.x` * 100) %>% 
  mutate("trim_reads_diff" = `Total number of reads after trim:.x` - `Total number of reads after trim:.y`) %>% 
  mutate("perc_trimmed_reads" = `Total number of reads after trim:.y` / `Total number of reads after trim:.x` * 100 ) %>%   
  select(Sample, 
         same_major_geno,
         `Total number of reads before trim:.x`,
         `Total number of reads before trim:.y`,
         `Total number of reads after trim:.x`,
         `Total number of reads after trim:.y`,
         trim_reads_diff,
         perc_trimmed_reads, 
         #total_classified_reads,
         `Number of mapped reads:.x`,
         `Number of mapped reads:.y`,
         diff_mapped_reads_with_dup,
         perc_diff_mapped_reads_with_dup,
         everything()) %>% 
  mutate("same_NS34A" = case_when(
    NS34A_short.x == NS34A_short.y ~ "YES",
    is.na(NS34A_short.x) & is.na(NS34A_short.y) ~ "YES",
    .default = "NO"
  )) %>% 
  mutate("same_NS5A" = case_when(
    NS5A_short.x == NS5A_short.y ~ "YES",
    is.na(NS5A_short.x) & is.na(NS5A_short.y) ~ "YES",
    .default = "NO"
  )) %>% 
  mutate("same_NS5B" = case_when(
    NS5B_short.x == NS5B_short.y ~ "YES",
    is.na(NS5B_short.x) & is.na(NS5B_short.y) ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(glecaprevir = case_when(
    glecaprevir.x == glecaprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(grazoprevir = case_when(
    grazoprevir.x == grazoprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(paritaprevir = case_when(
    paritaprevir.x == paritaprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(voxilaprevir = case_when(
    voxilaprevir.x == voxilaprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(daclatasvir = case_when(
    daclatasvir.x == daclatasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(elbasvir = case_when(
    elbasvir.x == elbasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(ledipasvir = case_when(
    ledipasvir.x == ledipasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(ombitasvir = case_when(
    ombitasvir.x == ombitasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(pibrentasvir = case_when(
    pibrentasvir.x == pibrentasvir.y ~ "YES",
    .default = "NO"
  )) %>%
  mutate(velpatasvir = case_when(
    velpatasvir.x == velpatasvir.y ~ "YES",
    .default = "NO"
  )) %>%
  mutate(dasabuvir = case_when(
    dasabuvir.x == dasabuvir.y ~ "YES",
    .default = "NO"
  )) %>%
  mutate(sofosbuvir = case_when(
    sofosbuvir.x == sofosbuvir.y ~ "YES",
    .default = "NO"
  ))

#write_csv(comparison, file = "N:/Virologi/JonBrate/2024.06.24_viralseq_routine_comparison.csv")
save(comparison, file = "/home/jonr/ny_hcv_rutine.RData")

# See ny_hcv_rutine.qmd

# Statistics



# Some are not analysed with "routine". Remove those
comparison <- comparison %>% 
  filter(!is.na(`Majority genotype:.y`))

# Alle har lik genotype
# Alle er like på typbar
# 87 av 93 typbare prøver har helt lik GLUE-rapport (dvs. alle NS34A short osv. har de samme mutasjonene. Men kan variere for enkeltmedikamenter, men alle de samme mutasjonene)
# 6 av 93 typbare har ulik GLUE-rapport. 

# Ulik genotype
comparison %>% filter(same_major_geno == "NO") %>% View()
  distinct() %>% 
  select(`Majority quality:.x`, `Majority quality:.y`, `Number of mapped reads without duplicates:.x`, `Number of mapped reads without duplicates:.y`) %>% 
  print(n=100)
  
# Typbare med nytt script, ikke typbare med gammelt
comparison %>% 
    distinct() %>% 
    filter(`Majority quality:.x` == "YES" | `Majority quality:.y` == "Ikke typbar") %>% 
    select(`Majority quality:.x`, `Majority quality:.y`) %>% print(n=300)

# Samme genotype, typbare begge script, like GLUE-rapporter
comparison %>% 
  distinct() %>% 
  filter(same_major_geno == "YES") %>% 
  filter(is.na(`Majority quality:.y`)) %>% 
  filter(`Majority quality:.x` == "YES") %>% 
  filter(same_NS34A == "YES" & same_NS5A == "YES" & same_NS5B == "YES") %>% 
  View()


# Samme genotype, typbar begge script, ulik GLUE-rapport
diff_glue <- comparison %>% 
  distinct() %>% 
  filter(same_major_geno == "YES") %>% 
  filter(is.na(`Majority quality:.y`)) %>% 
  filter(`Majority quality:.x` == "YES") %>% 
  filter(same_NS34A == "NO" | same_NS5A == "NO" | same_NS5B == "NO")
write_csv(diff_glue, file = "N:/Virologi/JonBrate/2024.06.24_diff_glue.csv")



# Investigate GLUE differences --------------------------------------------


# Rsamtools

# HCV122023C3 80K
library(Rsamtools)
library(GenomicRanges)


# 80K i NS34A - posisjon 3503-3511
# Mutasjon til A i posisjon 3506 gir K aa (Lysine)

HCV122023C3_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/HCV122023C3.1a_HQ850279_major.major.nodup.bam" 
HCV122023C3_N_bam <- scanBam(HCV122023C3_N_bam_file, param = )
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV122023C3_N_bam_file)
param <- ScanBamParam(which=GRanges("1a_HQ850279", IRanges(start=3506, end=3508)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("HCV122023C3 - K/Q80K - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV122023C3_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run917_HCV_summaries/bam/HCV122023C3_1a_tanoti_vbest_sorted.marked.bam" 
bf_G <- BamFile(HCV122023C3_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("HCV122023C3 - K/Q80K - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/HCV122023C3_res.RData")
res %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(0, 90, 10)) +
  ggtitle("HCV122023C3 - K/Q80K") +
  xlab("Position") +
  ylab("Reads") +
  facet_wrap(~Mapper)

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
    geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
             position = position_dodge2(padding = 0.01, preserve = "single"),
             width = 0.7) +
    scale_y_continuous(breaks = seq(0, 160, 10)) +
    scale_x_continuous(breaks = c(3506, 3507, 3508)) +
    ggtitle("HCV122023C3 - K/Q80K") +
    labs(x = "Genome position", y = "Reads", title = "HCV122023C3 - K/Q80K", fill = "Nucleotide") +  # Change legend title
    facet_wrap(~Mapper) +
    theme_minimal() +
    theme(
      text = element_text(size = 16),          # General font size
      axis.text = element_text(size = 14),     # Axis labels
      axis.title = element_text(size = 16),    # Axis titles
      strip.text = element_text(size = 18),    # Facet labels
      legend.title = element_text(size = 16),  # Legend title
      legend.text = element_text(size = 14),   # Legend text
      plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
      panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
      plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
    )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/HCV122023C3_KQ80K.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing

# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/HCV122023C3_total_counts.RData")


## HCV122023E4 C/N316N
HCV122023E4_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/HCV122023E4.1b_EU781827_major.major.nodup.bam" 
HCV122023E4_N_bam <- scanBam(HCV122023E4_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV122023E4_N_bam_file)
param <- ScanBamParam(which=GRanges("1b_EU781827", IRanges(start=8518, end=8520)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("HCV122023E4 - C/N316N - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV122023E4_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run917_HCV_summaries/bam/HCV122023E4_1b_tanoti_vbest_sorted.marked.bam" 
bf_G <- BamFile(HCV122023E4_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("HCV122023E4 - C/N316N - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/HCV122023E4_res.RData")
res %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(0, 90, 10)) +
  ggtitle("HCV122023C3 - K/Q80K") +
  xlab("Position") +
  ylab("Reads") +
  facet_wrap(~Mapper)

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 20, 2)) +
  scale_x_continuous(breaks = c(8518, 8519, 8520)) +
  labs(x = "Genome position", y = "Reads", title = "HCV122023E4 - C/N316N", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/HCV122023E4_CN316N.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing


# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/HCV122023E4_total_counts.RData")

## HCV092023kitukG3 - Y93H
# Y93H er i NS5A
HCV092023kitukG3_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/HCV092023kitukG3.3a_D17763_major.major.nodup.bam" 
HCV092023kitukG3_N_bam <- scanBam(HCV092023kitukG3_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV092023kitukG3_N_bam_file)
param <- ScanBamParam(which=GRanges("3a_D17763", IRanges(start=6550, end=6552)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("HCV092023kitukG3 - Y93H - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV092023kitukG3_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run898_HCV_summaries/bam/HCV092023kitukG3_3a_tanoti_vbest_sorted.marked.bam"
bf_G <- BamFile(HCV092023kitukG3_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("HCV092023kitukG3 - Y93H - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/HCV092023kitukG3_res.RData")
res %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(0, 90, 10)) +
  ggtitle("HCV122023C3 - K/Q80K") +
  xlab("Position") +
  ylab("Reads") +
  facet_wrap(~Mapper)

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  scale_x_continuous(breaks = c(6550, 6551, 6552)) +
  labs(x = "Genome position", y = "Reads", title = "HCV092023kitukG3 - Y93H", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/HCV092023kitukG3_Y93H.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing

# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/HCV092023kitukG3_total_counts.RData")

## HCV112023A3 G/S556G
# NS5B
HCV112023A3_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/HCV112023A3.1b_EU781827_major.major.nodup.bam" 
HCV112023A3_N_bam <- scanBam(HCV112023A3_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV112023A3_N_bam_file)
param <- ScanBamParam(which=GRanges("1b_EU781827", IRanges(start=9241, end=9243)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("HCV112023A3 - G/S556G - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV112023A3_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run911_HCV_summaries/bam/HCV112023A3_1b_tanoti_vbest_sorted.marked.bam"
bf_G <- BamFile(HCV112023A3_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("HCV112023A3 - G/S556G - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/HCV112023A3_res.RData")
res %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(0, 90, 10)) +
  ggtitle("HCV112023A3 - G/S556G") +
  xlab("Position") +
  ylab("Reads") +
  facet_wrap(~Mapper)

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 20, 2)) +
  scale_x_continuous(breaks = c(9241, 9242, 9243)) +
  labs(x = "Genome position", y = "Reads", title = "HCV112023A3 - G/S556G", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/HCV112023A3_GS556G.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing

# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/HCV112023A3_total_counts.RData")

## HCV082023B4 Y56Y;Q168Q;I/V170I
# Y56Y i NS34A
HCV082023B4_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/HCV082023B4.3a_D17763_major.major.nodup.bam" 
HCV082023B4_N_bam <- scanBam(HCV082023B4_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV082023B4_N_bam_file)
param <- ScanBamParam(which=GRanges("3a_D17763", IRanges(start=3601, end=3603)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("HCV082023B4 - Y56Y - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV082023B4_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run891_HCV_summaries/bam/HCV082023B4_3a_tanoti_vbest_sorted.marked.bam"
bf_G <- BamFile(HCV082023B4_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("HCV082023B4 - Y56Y - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/HCV082023B4_Y56Y_res.RData")
res %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(0, 90, 10)) +
  ggtitle("HCV112023A3 - G/S556G") +
  xlab("Position") +
  ylab("Reads") +
  facet_wrap(~Mapper)

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 70, 10)) +
  scale_x_continuous(breaks = c(3601, 3602, 3603)) +
  labs(x = "Genome position", y = "Reads", title = "HCV082023B4 - Y56Y", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/HCV082023B4_Y56Y.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing

# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/HCV082023B4_Y56Y_total_counts.RData")

# Q168Q i NS34A
HCV082023B4_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/HCV082023B4.3a_D17763_major.major.nodup.bam" 
HCV082023B4_N_bam <- scanBam(HCV082023B4_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV082023B4_N_bam_file)
param <- ScanBamParam(which=GRanges("3a_D17763", IRanges(start=3937, end=3939)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("HCV082023B4 - Q168Q - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV082023B4_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run891_HCV_summaries/bam/HCV082023B4_3a_tanoti_vbest_sorted.marked.bam"
bf_G <- BamFile(HCV082023B4_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("HCV082023B4 - Q168Q - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/HCV082023B4_Q168Q_res.RData")

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 80, 10)) +
  scale_x_continuous(breaks = c(3937, 3938, 3939)) +
  labs(x = "Genome position", y = "Reads", title = "HCV082023B4 - Q168Q", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/HCV082023B4_Q168Q.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing


# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/HCV082023B4_Q168Q_total_counts.RData")

# I/V170I i NS34A
HCV082023B4_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/HCV082023B4.3a_D17763_major.major.nodup.bam" 
HCV082023B4_N_bam <- scanBam(HCV082023B4_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV082023B4_N_bam_file)
param <- ScanBamParam(which=GRanges("3a_D17763", IRanges(start=3943, end=3945)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("HCV082023B4 - I/V170I - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV082023B4_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run891_HCV_summaries/bam/HCV082023B4_3a_tanoti_vbest_sorted.marked.bam"
bf_G <- BamFile(HCV082023B4_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("HCV082023B4 - I/V170I - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/HCV082023B4_IV170I_res.RData")

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 90, 10)) +
  scale_x_continuous(breaks = c(3943, 3944, 3945)) +
  labs(x = "Genome position", y = "Reads", title = "HCV082023B4 - IV170I", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/HCV082023B4_IV170I.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing


# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/HCV082023B4_IV170I_total_counts.RData")

## 2190609-HCV - A/T/V150V
# Dette er NS5B. Run: 20240416-01. Referansen er: 3a_D17763
HCV2190609_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/2190609-HCV.3a_D17763_major.major.nodup.bam" 
HCV2190609_bam <- scanBam(HCV2190609_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV2190609_N_bam_file)
param <- ScanBamParam(which=GRanges("3a_D17763", IRanges(start=8077, end=8079)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2-default") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("2190609-HCV - A/T/V150V - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV2190609_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2024/NGS_SEQ-20240416-01_summaries/bam/2190609-HCV_3a_tanoti_vbest_sorted.marked.bam"
bf_G <- BamFile(HCV2190609_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("2190609-HCV - A/T/V150V - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/2190609-HCV_V150V_res.RData")

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 50, 10)) +
  scale_x_continuous(breaks = c(8077, 8078, 8079)) +
  labs(x = "Genome position", y = "Reads", title = "2190609-HCV - V150V", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/2190609-HCV_V150V.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing

# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/2190609-HCV_V150V_total_counts.RData")

# Investigate the effect of changing Bowtie2 gap penalties
HCV2190609_new_params_bam_file <- "/home/jonr/ny_hcv_rutine_bams/teste_gap_penalty/new.bam" 
HCV2190609_new_params_bam <- scanBam(HCV2190609_new_params_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N2 <- BamFile(HCV2190609_new_params_bam_file)
param <- ScanBamParam(which=GRanges("3a_D17763", IRanges(start=8077, end=8079)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N2 <- pileup(bf_N2, scanBamParam=param, pileupParam = p_param)
res_N2 <- tibble(res_N2) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2-endret") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N2 %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("2190609-HCV - A/T/V150V - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

res_new <- bind_rows(res_N, res_N2, res_G)

res_new %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 50, 10)) +
  scale_x_continuous(breaks = c(8077, 8078, 8079)) +
  labs(x = "Genome position", y = "Reads", title = "2190609-HCV - V150V", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/2190609-HCV_V150V_endret_parameter.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing


## 2192514-HCV - C/N316N
# Dette er NS5B. Run: 20240416-01. Referansen er: 1b_EU781827
HCV2192514_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/2192514-HCV.1b_EU781827_major.major.nodup.bam" 
HCV2192514_bam <- scanBam(HCV2192514_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV2192514_N_bam_file)
param <- ScanBamParam(which=GRanges("1b_EU781827", IRanges(start=8518, end=8520)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("2192514-HCV - C/N316N - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV2192514_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2024/NGS_SEQ-20240416-01_summaries/bam/2192514-HCV_1b_tanoti_vbest_sorted.marked.bam"
bf_G <- BamFile(HCV2192514_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("2192514-HCV - C/N316N - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/HCV2192514-HCV_res.RData")

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 20, 5)) +
  scale_x_continuous(breaks = c(8518, 8519, 8520)) +
  labs(x = "Genome position", y = "Reads", title = "HCV2192514-HCV - C/N316N", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/HCV2192514-HCV_CN316N.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing

# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/HCV2192514-HCV_total_counts.RData")


# Minor default Bowtie2 ---------------------------------------------------



## Minor genotype
viralseq_HCVkoinf2022 <- read_csv(paste0("/home/jonr/Prosjekter/viralseq/hcv_koinf/summarize/Genotype_mapping_summary_long_LW_import_with_glue_minor.csv")) %>% 
  # Add run name
  add_column("Run" = "HCVkoinf2022") %>% 
  mutate(glecaprevir_mut_short = as.character(glecaprevir_mut_short)) %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  mutate(voxilaprevir_mut_short = as.character(voxilaprevir_mut_short)) %>% 
  # Some samples are duplicated with identical entries in all columns
  distinct()

viralseq <- viralseq_HCVkoinf2022
vs_cols <- colnames(viralseq)

routine_HCVkoinf2022 <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2024/HCVkoinf2022_summaries/HCVkoinf2022_summary_with_glue.tsv"), 
                                   delim = "\t", escape_double = FALSE, 
                                   col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                                    `Percent mapped reads of trimmed:` = col_character(), 
                                                    `Percent covered:` = col_character(), 
                                                    `Percent covered above depth=5 without duplicates:` = col_character(), 
                                                    `Percent covered above depth=9 without duplicates:` = col_character(), 
                                                    `Percent covered minor:` = col_character(), 
                                                    `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                                    `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                                    glecaprevir_mut_short = col_character(),
                                                    paritaprevir_mut_short = col_character(),
                                                    voxilaprevir_mut_short = col_character()), 
                                   trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "HCVkoinf2022") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  )) 

routine_HCVkoinf2022_rest <- read_delim(paste0(prefix, "Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2024/HCVkoinf2022_summaries/HCVkoinf2022.rest_summary_with_glue.tsv"), 
                                        delim = "\t", escape_double = FALSE, 
                                        col_types = cols(`Average depth minor without duplicates:` = col_character(),
                                                         `Percent mapped reads of trimmed:` = col_character(), 
                                                         `Percent covered:` = col_character(), 
                                                         `Percent covered above depth=5 without duplicates:` = col_character(), 
                                                         `Percent covered above depth=9 without duplicates:` = col_character(), 
                                                         `Percent covered minor:` = col_character(), 
                                                         `Percent covered above depth=5 minor without duplicates:` = col_character(), 
                                                         `Percent covered above depth=9 minor without duplicates:` = col_character(), 
                                                         glecaprevir_mut_short = col_character(),
                                                         paritaprevir_mut_short = col_character(),
                                                         voxilaprevir_mut_short = col_character()), 
                                        trim_ws = TRUE) %>% 
  # Convert comma to "." and then to double
  mutate(`Percent mapped reads of trimmed:` = str_replace(`Percent mapped reads of trimmed:`, ",", "\\.")) %>% 
  mutate(`Percent mapped reads of trimmed:` = as.double(`Percent mapped reads of trimmed:`)) %>% 
  mutate(`Percent covered:` = str_replace(`Percent covered:`, ",", "\\.")) %>% 
  mutate(`Percent covered:` = as.double(`Percent covered:`)) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = str_replace(`Percent covered above depth=5 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=5 without duplicates:` = as.double(`Percent covered above depth=5 without duplicates:`)) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = str_replace(`Percent covered above depth=9 without duplicates:`, ",", "\\.")) %>% 
  mutate(`Percent covered above depth=9 without duplicates:` = as.double(`Percent covered above depth=9 without duplicates:`)) %>% 
  # Add run name
  add_column("Run" = "HCVkoinf2022_rest") %>% 
  # Create typbar/ikke typbar because NA is filled in on import
  mutate(`Majority quality:` = case_when(
    `Percent covered:` >= 10 & `Average depth without duplicates:` >=2 ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(`Minor quality:` = case_when(
    `Percent covered minor:` >= 10 & `Average depth minor without duplicates:` >=2 ~ "YES",
    .default = "NO"
  ))

routine <- bind_rows(routine_HCVkoinf2022, routine_HCVkoinf2022_rest)

# The glue reports are for the Major genotype. Remove columns and join with minor Glue reports created outside the routine
routine <- routine %>% 
  select(1:23)

# Read the minor GLUE report
glue_minor <- read_tsv("/home/jonr/ny_hcv_rutine_bams/rutine_json/GLUE_collected_report.tsv") %>% 
  # Trim the sample name
  separate(Sample, into = c("Sample"), sep = "_", remove = TRUE) %>% 
  # Fix a sample name
  mutate(Sample = str_replace(Sample, "Virus220204", "HCVVirus220204"))

# Join with routine
routine <- left_join(routine, glue_minor, by = "Sample")
rt_cols <- colnames(routine)

# Compare columns
setdiff(vs_cols, rt_cols)

# Join data - keep only samples present in both viralseq and routine
joined_minor <- inner_join(viralseq, routine, by = "Sample") 

comparison_minor <- joined_minor %>% 
  mutate("same_minor_geno" = case_when(
    `Most abundant minority genotype:.x` == `Most abundant minority genotype:.y` ~ "YES",
    .default = "NO"
  )) %>% 
  mutate("diff_mapped_reads_with_dup" = `Number of mapped reads:.x` - `Number of mapped reads:.y`) %>% 
  mutate("perc_diff_mapped_reads_with_dup" = `Number of mapped reads:.y` / `Number of mapped reads:.x` * 100) %>% 
  mutate("diff_mapped_reads_no_dup" = `Number of mapped reads minor without duplicates:.x` - `Number of mapped reads minor without duplicates:.y`) %>% 
  mutate("perc_diff_mapped_reads_no_dup" = `Number of mapped reads minor without duplicates:.y` / `Number of mapped reads minor without duplicates:.x` * 100) %>% 
  mutate("trim_reads_diff" = `Total number of reads after trim:.x` - `Total number of reads after trim:.y`) %>% 
  mutate("perc_trimmed_reads" = `Total number of reads after trim:.y` / `Total number of reads after trim:.x` * 100 ) %>%   
  select(Sample, 
         same_minor_geno,
         `Total number of reads before trim:.x`,
         `Total number of reads before trim:.y`,
         `Total number of reads after trim:.x`,
         `Total number of reads after trim:.y`,
         trim_reads_diff,
         perc_trimmed_reads, 
         #total_classified_reads,
         `Number of mapped reads:.x`,
         `Number of mapped reads:.y`,
         diff_mapped_reads_with_dup,
         perc_diff_mapped_reads_with_dup,
         everything()) %>% 
  mutate("same_NS34A" = case_when(
    NS34A_short.x == NS34A_short.y ~ "YES",
    is.na(NS34A_short.x) & is.na(NS34A_short.y) ~ "YES",
    .default = "NO"
  )) %>% 
  mutate("same_NS5A" = case_when(
    NS5A_short.x == NS5A_short.y ~ "YES",
    is.na(NS5A_short.x) & is.na(NS5A_short.y) ~ "YES",
    .default = "NO"
  )) %>% 
  mutate("same_NS5B" = case_when(
    NS5B_short.x == NS5B_short.y ~ "YES",
    is.na(NS5B_short.x) & is.na(NS5B_short.y) ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(glecaprevir = case_when(
    glecaprevir.x == glecaprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(grazoprevir = case_when(
    grazoprevir.x == grazoprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(paritaprevir = case_when(
    paritaprevir.x == paritaprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(voxilaprevir = case_when(
    voxilaprevir.x == voxilaprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(daclatasvir = case_when(
    daclatasvir.x == daclatasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(elbasvir = case_when(
    elbasvir.x == elbasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(ledipasvir = case_when(
    ledipasvir.x == ledipasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(ombitasvir = case_when(
    ombitasvir.x == ombitasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(pibrentasvir = case_when(
    pibrentasvir.x == pibrentasvir.y ~ "YES",
    .default = "NO"
  )) %>%
  mutate(velpatasvir = case_when(
    velpatasvir.x == velpatasvir.y ~ "YES",
    .default = "NO"
  )) %>%
  mutate(dasabuvir = case_when(
    dasabuvir.x == dasabuvir.y ~ "YES",
    .default = "NO"
  )) %>%
  mutate(sofosbuvir = case_when(
    sofosbuvir.x == sofosbuvir.y ~ "YES",
    .default = "NO"
  ))

#write_csv(comparison, file = "N:/Virologi/JonBrate/2024.06.24_viralseq_routine_comparison.csv")
save(comparison_minor, file = "/home/jonr/ny_hcv_rutine_minor.RData")


## HCV082022A3 - S556G - NS5B
HCV122022B5_N_bam_file <- "/home/jonr/Prosjekter/viralseq/hcv_koinf/samtools/HCV122022B5.1a_HQ850279_minor.minor.nodup.bam" 
HCV122022B5_N_bam <- scanBam(HCV122022B5_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV122022B5_N_bam_file)
param <- ScanBamParam(which=GRanges("1a_HQ850279", IRanges(start=6550, end=6552)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("HCV092023kitukG3 - Y93H - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV092023kitukG3_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Run898_HCV_summaries/bam/HCV092023kitukG3_3a_tanoti_vbest_sorted.marked.bam"
bf_G <- BamFile(HCV092023kitukG3_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("HCV092023kitukG3 - Y93H - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/HCV092023kitukG3_res.RData")
res %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(0, 90, 10)) +
  ggtitle("HCV122023C3 - K/Q80K") +
  xlab("Position") +
  ylab("Reads") +
  facet_wrap(~Mapper)

# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/HCV092023kitukG3_total_counts.RData")


# Changed Bowtie2 params --------------------------------------------------
## Undersøke effekten av å endre Bowtie2 gap penalties
## Først undersøke på major prøvene

viralseq_20240416_changed <- read_csv(paste0(prefix, "Virologi/JonBrate/NGS_SEQ-20240416-01_endret/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "20240416-01") %>% 
  mutate(glecaprevir_mut_short = as.character(glecaprevir_mut_short)) %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  mutate(voxilaprevir_mut_short = as.character(voxilaprevir_mut_short)) %>% 
  # Some samples are duplicated with identical entries in all columns
  distinct()

viralseq_changed <- bind_rows(viralseq_20240416_changed)

# Join data - keep only samples present in both viralseq and routine
joined_changed <- inner_join(viralseq_changed, routine, by = "Sample") 

# Compare
comparison_changed <- joined_changed %>% 
  mutate("same_major_geno" = case_when(
    `Majority genotype:.x` == `Majority genotype:.y` ~ "YES",
    .default = "NO"
  )) %>% 
  mutate("diff_mapped_reads_with_dup" = `Number of mapped reads:.x` - `Number of mapped reads:.y`) %>% 
  mutate("perc_diff_mapped_reads_with_dup" = `Number of mapped reads:.y` / `Number of mapped reads:.x` * 100) %>% 
  mutate("diff_mapped_reads_no_dup" = `Number of mapped reads minor without duplicates:.x` - `Number of mapped reads minor without duplicates:.y`) %>% 
  mutate("perc_diff_mapped_reads_no_dup" = `Number of mapped reads minor without duplicates:.y` / `Number of mapped reads minor without duplicates:.x` * 100) %>% 
  mutate("trim_reads_diff" = `Total number of reads after trim:.x` - `Total number of reads after trim:.y`) %>% 
  mutate("perc_trimmed_reads" = `Total number of reads after trim:.y` / `Total number of reads after trim:.x` * 100 ) %>%   
  select(Sample, 
         same_major_geno,
         `Total number of reads before trim:.x`,
         `Total number of reads before trim:.y`,
         `Total number of reads after trim:.x`,
         `Total number of reads after trim:.y`,
         trim_reads_diff,
         perc_trimmed_reads, 
         #total_classified_reads,
         `Number of mapped reads:.x`,
         `Number of mapped reads:.y`,
         diff_mapped_reads_with_dup,
         perc_diff_mapped_reads_with_dup,
         everything()) %>% 
  mutate("same_NS34A" = case_when(
    NS34A_short.x == NS34A_short.y ~ "YES",
    is.na(NS34A_short.x) & is.na(NS34A_short.y) ~ "YES",
    .default = "NO"
  )) %>% 
  mutate("same_NS5A" = case_when(
    NS5A_short.x == NS5A_short.y ~ "YES",
    is.na(NS5A_short.x) & is.na(NS5A_short.y) ~ "YES",
    .default = "NO"
  )) %>% 
  mutate("same_NS5B" = case_when(
    NS5B_short.x == NS5B_short.y ~ "YES",
    is.na(NS5B_short.x) & is.na(NS5B_short.y) ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(glecaprevir = case_when(
    glecaprevir.x == glecaprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(grazoprevir = case_when(
    grazoprevir.x == grazoprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(paritaprevir = case_when(
    paritaprevir.x == paritaprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(voxilaprevir = case_when(
    voxilaprevir.x == voxilaprevir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(daclatasvir = case_when(
    daclatasvir.x == daclatasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(elbasvir = case_when(
    elbasvir.x == elbasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(ledipasvir = case_when(
    ledipasvir.x == ledipasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(ombitasvir = case_when(
    ombitasvir.x == ombitasvir.y ~ "YES",
    .default = "NO"
  )) %>% 
  mutate(pibrentasvir = case_when(
    pibrentasvir.x == pibrentasvir.y ~ "YES",
    .default = "NO"
  )) %>%
  mutate(velpatasvir = case_when(
    velpatasvir.x == velpatasvir.y ~ "YES",
    .default = "NO"
  )) %>%
  mutate(dasabuvir = case_when(
    dasabuvir.x == dasabuvir.y ~ "YES",
    .default = "NO"
  )) %>%
  mutate(sofosbuvir = case_when(
    sofosbuvir.x == sofosbuvir.y ~ "YES",
    .default = "NO"
  ))

save(comparison_changed, file = "/home/jonr/ny_hcv_rutine_changed.RData")


# Changed params - Investigate GLUE differences ---------------------------

## 2190609-HCV - A/T/V150V
# Dette er NS5B. Run: 20240416-01. Referansen er: 3a_D17763
HCV2190609_N_bam_file <- "/home/jonr/ny_hcv_rutine_bams/endret_gap_penalty/2190609-HCV.3a_D17763_major.major.nodup.bam"
HCV2190609_bam <- scanBam(HCV2190609_N_bam_file)
# https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
bf_N <- BamFile(HCV2190609_N_bam_file)
param <- ScanBamParam(which=GRanges("3a_D17763", IRanges(start=8077, end=8079)))
p_param <- PileupParam(max_depth=100000, distinguish_strand=TRUE)
res_N <- pileup(bf_N, scanBamParam=param, pileupParam = p_param)
res_N <- tibble(res_N) %>% 
  # Create colums for + and - strands
  tidyr::unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Bowtie2-default") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_N %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_N$count, na.rm = T), max(res_N$count, na.rm = T), by = 10)) +
  ggtitle("2190609-HCV - A/T/V150V - Bowtie2") +
  xlab("Position") +
  ylab("Reads")

HCV2190609_G_bam_file <- "/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2024/NGS_SEQ-20240416-01_summaries/bam/2190609-HCV_3a_tanoti_vbest_sorted.marked.bam"
bf_G <- BamFile(HCV2190609_G_bam_file)
res_G <- pileup(bf_G, scanBamParam=param, pileupParam = p_param)
res_G <- tibble(res_G) %>% 
  # Create colums for + and - strands
  unite("col_name", c(pos, strand), sep = " ", remove = FALSE) %>% 
  # Fill missing observarions with NA
  complete(col_name, nucleotide) %>% 
  # Add mapper column
  add_column("Mapper" = "Tanoti") %>% 
  # Remove = and - positions
  filter(nucleotide == "A" | nucleotide == "T" | nucleotide == "C" | nucleotide == "G")
res_G %>% ggplot() +
  geom_col(aes(x = col_name, y = count, group = nucleotide, fill = nucleotide), 
           position = "dodge") +
  scale_y_continuous(breaks = seq(min(res_G$count, na.rm = T), max(res_G$count, na.rm = T), by = 10)) +
  ggtitle("2190609-HCV - A/T/V150V - Tanoti") +
  xlab("Position") +
  ylab("Reads")

res <- bind_rows(res_N, res_G)
save(res, file = "/home/jonr/2190609-HCV_V150V_res.RData")

# Plot coverage regardless of strand
res %>% 
  group_by(pos, nucleotide, Mapper) %>% 
  summarize(total_count = sum(count, na.rm = T)) %>% 
  ggplot() +
  geom_col(aes(x = pos, y = total_count, group = nucleotide, fill = nucleotide), 
           position = position_dodge2(padding = 0.01, preserve = "single"),
           width = 0.7) +
  scale_y_continuous(breaks = seq(0, 50, 10)) +
  scale_x_continuous(breaks = c(8077, 8078, 8079)) +
  labs(x = "Genome position", y = "Reads", title = "2190609-HCV - V150V", fill = "Nucleotide") +  # Change legend title
  facet_wrap(~Mapper) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # General font size
    axis.text = element_text(size = 14),     # Axis labels
    axis.title = element_text(size = 16),    # Axis titles
    strip.text = element_text(size = 18),    # Facet labels
    legend.title = element_text(size = 16),  # Legend title
    legend.text = element_text(size = 14),   # Legend text
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Margins: top, right, bottom, left
    panel.background = element_rect(fill = "white", colour = "white"),  # Panel background to white
    plot.background = element_rect(fill = "white", colour = "white")    # Plot background to white
  )

# Set the plot dimensions in centimeters for A4
ggsave("/home/jonr/2190609-HCV_V150V.png", 
       width = 27, height = 14,  # Full width, adjust height as needed
       units = "cm",               # Use centimeters for dimensions
       dpi = 300)                  # High resolution for printing

# Calculate the count difference between Bowtie2 and Tanoti
total_counts <- res %>% 
  filter(!is.na(count)) %>%
  select(pos, count, Mapper) %>% 
  # Count per position regardless of strand and nucleotide
  group_by(pos, Mapper) %>% 
  summarize(total_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = Mapper, values_from = total_count, values_fill = 0) %>% 
  mutate(diff = Bowtie2 - Tanoti)
save(total_counts, file = "/home/jonr/2190609-HCV_V150V_total_counts.RData")



### Undersøke om det er noen prøver som har fått minor med Bowtie2 og ikke Tanoti

### Undersøke prøver med minor gt i rutinen
Først trenger jeg å vite total mapped reads mot alle referanser i den første mappingen. Hvor finner jeg det?
  Deretter coverage mot potensiell minor gt. Dette er i parsefirstmapping csv-fil.
Må nok kjøre pipelinen på nytt... Men jeg kan ihvertfall filtere ut prøver med flere enn 50000 reads mappet til major og minor coverage min 5%
I tillegg kan jeg sjekke at minor genotype må være en annen enn 

# Candidates
major_mapping <- viralseq %>% filter(`Number of mapped reads:` > 49999)

# Lese inn alle "parsefirstmapping"-filene:
viralseq_Run879_csv <- list.files(path = paste0(prefix, "Virologi/JonBrate/Run879/parsefirstmapping/"),
                                  pattern = ".csv$",
                                  full.names = TRUE)
# Samples with min 5% minor coverage
viralseq_Run879_minor <- read_csv(viralseq_Run879_csv) %>% filter(minor_cov >= 5)
# Add run name
  add_column("Run" = "Run879") %>%
  # Remove duplicated entries
  distinct()
viralseq_Run882 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run882_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run882") %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run891 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run891_MiSeq_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run891") %>% 
  # Remove sample that has be re-run in Run897 with more mapped reads
  filter(Sample != "HCV082023A1") %>% 
  filter(Sample != "HCV082023F3") %>% 
  filter(Sample != "HCV082023A3") %>% 
  filter(Sample != "HCV082023B1") %>% 
  filter(Sample != "HCV082023B2") %>% 
  filter(Sample != "HCV082023D2") %>% 
  filter(Sample != "HCV082023D3") %>% 
  filter(Sample != "HCV082023E1") %>% 
  filter(Sample != "HCV082023E3") %>% 
  filter(Sample != "HCV082023G2") %>% 
  filter(Sample != "HCV082023H1") %>% 
  # Remove duplicated samples
  distinct()
viralseq_Run897 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run897_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run897") %>% 
  # Remove sample that has be run in Run891 with more mapped reads
  filter(Sample != "HCV082023C3") %>% 
  filter(Sample != "HCV082023F1") %>% 
  filter(Sample != "HCV082023H2") %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run898 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run898_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run898") %>% 
  mutate(glecaprevir_mut_short = as.character(glecaprevir_mut_short)) %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  mutate(NS34A_short = as.character(NS34A_short)) %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run902 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run902_virus_og_bakt/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run902") %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run911 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run911_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run911") %>% 
  # Remove duplicated entries
  distinct()
viralseq_Run917 <- read_csv(paste0(prefix, "Virologi/JonBrate/Run917_Virus/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "Run917") %>% 
  # Remove duplicated entries
  distinct()
viralseq_20240416 <- read_csv(paste0(prefix, "Virologi/JonBrate/NGS_SEQ-20240416-01/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "20240416-01") %>% 
  mutate(glecaprevir_mut_short = as.character(glecaprevir_mut_short)) %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  mutate(voxilaprevir_mut_short = as.character(voxilaprevir_mut_short)) %>% 
  # Some samples are duplicated with identical entries in all columns
  distinct()
viralseq_20240617 <- read_csv(paste0(prefix, "Virologi/JonBrate/NGS_SEQ-20240617-01/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "20240627-01") %>% 
  # Some samples are duplicated with identical entries in all columns
  distinct()
viralseq_20240814 <- read_csv(paste0(prefix, "Virologi/JonBrate/NGS_SEQ-20240814-01/summarize/Genotype_mapping_summary_long_LW_import.csv")) %>% 
  # Add run name
  add_column("Run" = "20240814-01") %>% 
  mutate(glecaprevir_mut_short = as.character(glecaprevir_mut_short)) %>% 
  mutate(paritaprevir_mut_short = as.character(paritaprevir_mut_short)) %>% 
  mutate(voxilaprevir_mut_short = as.character(voxilaprevir_mut_short)) %>% 
  # Some samples are duplicated with identical entries in all columns
  distinct()

viralseq <- bind_rows(viralseq_Run879, viralseq_Run882, viralseq_Run891, viralseq_Run897, viralseq_Run898, viralseq_Run902, viralseq_Run911, viralseq_Run917, viralseq_20240416, viralseq_20240617, viralseq_20240814)
vs_cols <- colnames(viralseq)
