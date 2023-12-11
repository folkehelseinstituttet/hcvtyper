library(tidyverse)

# TODO

# Hvorfor er Percent_reads_mapped_with_dups_majority alltid 100? Kan ikke bruke Raw total sequences fra samtools stats. Fordi bamfilen ikke inneholder unmapped reads. 
# Må bruke antallet trimmed og Kraken2-classified reads. Hvor kan jeg hente dette tallet...?
# kan lese kraken2 rapporten og ta antallet root classified reads
# Om Minority thresholds. Hvor mange reads uten duplikater er ok for å genotype? Bruke dette som kriterium for antall reads og coverage.

run911 <- read_csv("/home/jonr/Prosjekter/viralseq/Run911_HCV_tanoti/summarize/Genotype_mapping_summary_long.csv")

res2023 <- readxl::read_xlsx("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2023/Resultater_HCV_oversikt_2023.xlsx",
                             skip = 1,
                             sheet = "Resultater")


tmp <- left_join(run911, res2023, by = c("sampleName" = "Sample_ID")) %>% 
  mutate("samme_major" = case_when(
    Majority_genotype_mapping == `Majority genotype` ~ "YES",
    .default = "NO"
  )) 

tmp %>% 
  select(sampleName, samme_major, Majority_genotype_mapping, `Majority genotype`, `Typbar/ \r\nIkke typbar Minor`, Majority_cov_breadth_min_5, `Percent covered above depth=5`, 
         Reads_withdup_mapped_majority, Reads_nodup_mapped_majority, `Number of mapped reads`, `Number of mapped reads 2`) %>% View()


