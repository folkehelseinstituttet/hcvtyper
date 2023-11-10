#!/usr/bin/env Rscript

library(tidyverse)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript GLUE_json_parser.R <GLUE json report>", call.=FALSE)
}

json_file <- args[1]

# Create final data file
df_final <- tibble(
  "Sample" = character(), 
  "GLUE_genotype" = character(),
  "GLUE_subtype" = character(),
  "glecaprevir" = character(), 
  "glecaprevir_mut" = character(), 
  "glecaprevir_mut_short" = character(), 
  "grazoprevir" = character(), 
  "grazoprevir_mut" = character(), 
  "grazoprevir_mut_short" = character(), 
  "paritaprevir" = character(), 
  "paritaprevir_mut" = character(), 
  "paritaprevir_mut_short" = character(), 
  "voxilaprevir" = character(), 
  "voxilaprevir_mut" = character(), 
  "voxilaprevir_mut_short" = character(), 
  "NS34A" = character(),
  "NS34A_short" = character(),
  "daclatasvir" = character(), 
  "daclatasvir_mut" = character(), 
  "daclatasvir_mut_short" = character(), 
  "elbasvir" = character(), 
  "elbasvir_mut" = character(), 
  "elbasvir_mut_short" = character(), 
  "ledipasvir" = character(), 
  "ledipasvir_mut" = character(),
  "ledipasvir_mut_short" = character(),
  "ombitasvir" = character(), 
  "ombitasvir_mut" = character(),
  "ombitasvir_mut_short" = character(),
  "pibrentasvir" = character(),
  "pibrentasvir_mut" = character(),
  "pibrentasvir_mut_short" = character(),
  "velpatasvir" = character(), 
  "velpatasvir_mut" = character(), 
  "velpatasvir_mut_short" = character(), 
  "NS5A" = character(),
  "NS5A_short" = character(),
  "dasabuvir" = character(),
  "dasabuvir_mut" = character(),
  "dasabuvir_mut_short" = character(),
  "sofosbuvir" = character(),
  "sofosbuvir_mut" = character(),
  "sofosbuvir_mut_short" = character(),
  "NS5B" = character(),
  "NS5B_short" = character(),
  "HCV project version" = character(),
  "GLUE engine version" = character(),
  "PHE drug resistance extension version"  = character()
)

json <- read_json(json_file)

  # Only read the proper json GLUE reports (i.e. that there was a good sequence)
  if (names(json) == "phdrReport") {
    # Sample name
    sample <- unlist(strsplit(basename(json_file), "\\."))[[1]]
    
    # Få tak i genotype
    genotype <- json[["phdrReport"]][["samReferenceResult"]][["genotypingResult"]][["genotypeCladeCategoryResult"]][["shortRenderedName"]]
    
    # Få tak i subtype
    subtype <- json[["phdrReport"]][["samReferenceResult"]][["genotypingResult"]][["subtypeCladeCategoryResult"]][["shortRenderedName"]]
    
    # Versjoner:
    projectVersion <- json[["phdrReport"]][["projectVersion"]]
    extensionVersion <- json[["phdrReport"]][["extensionVersion"]]
    engineVersion <- json[["phdrReport"]][["engineVersion"]]
    
    # One row per sample
    # Create a temporary dataframe to populate
    df <- as.data.frame(matrix(nrow = 1, ncol = 48))
    colnames(df) <- c("Sample", 
                      "GLUE_genotype",
                      "GLUE_subtype",
                      "glecaprevir", 
                      "glecaprevir_mut",
                      "glecaprevir_mut_short",
                      "grazoprevir", 
                      "grazoprevir_mut", 
                      "grazoprevir_mut_short", 
                      "paritaprevir", 
                      "paritaprevir_mut", 
                      "paritaprevir_mut_short", 
                      "voxilaprevir", 
                      "voxilaprevir_mut", 
                      "voxilaprevir_mut_short", 
                      "NS34A",
                      "NS34A_short",
                      "daclatasvir", 
                      "daclatasvir_mut", 
                      "daclatasvir_mut_short", 
                      "elbasvir", 
                      "elbasvir_mut", 
                      "elbasvir_mut_short", 
                      "ledipasvir", 
                      "ledipasvir_mut", 
                      "ledipasvir_mut_short", 
                      "ombitasvir", 
                      "ombitasvir_mut", 
                      "ombitasvir_mut_short", 
                      "pibrentasvir",
                      "pibrentasvir_mut",
                      "pibrentasvir_mut_short",
                      "velpatasvir", 
                      "velpatasvir_mut", 
                      "velpatasvir_mut_short", 
                      "NS5A",
                      "NS5A_short",
                      "dasabuvir",
                      "dasabuvir_mut",
                      "dasabuvir_mut_short",
                      "sofosbuvir",
                      "sofosbuvir_mut",
                      "sofosbuvir_mut_short",
                      "NS5B",
                      "NS5B_short",
                      "HCV project version",
                      "GLUE engine version",
                      "PHE drug resistance extension version")
    
    df$Sample <- sample
    df$GLUE_genotype <- genotype
    df$GLUE_subtype <- subtype
    df$`HCV project version` <- projectVersion
    df$`GLUE engine version` <- engineVersion
    df$`PHE drug resistance extension version` <- extensionVersion 
    
    # Dette er underlisten for Drug Scores. Lengden av denne angir hvor mange drugs som er funnet.
    # Under drugScores så er det en ny liste for hver drug category
    if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]]) > 0) {
      for (i in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]])) {
        # Så er det en ny liste innenfor hver drug score igjen med drug for hver kategori. Denne heter drugAssessemnts
        for (k in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]])) {
                  # Hvis det er sufficient coverage (denne evaluerer til TRUE):
                  if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["sufficientCoverage"]]) {
                    if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]] == "No resistance") {
                      df[[json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]]]]  <- "No resistance"
                    } else {  
                      # Skrive inn resistance informasjonen for druget
                      df[[json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]]]]  <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
                      
                      # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i den
                      mut <- vector(mode = "character") # Create empty vector to hold mutations
                      mut_short <- vector(mode = "character") # Create empty vector to hold mutations
                      if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
                        for (n in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]])) {
                          mut <- c(mut, json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[n]][["displayStructure"]])
                          mut_short <- c(mut_short, json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[n]][["structure"]])
                        }
                      } 
                      if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
                        for (n in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]])) {
                          mut <- c(mut, json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[n]][["displayStructure"]])
                          mut_short <- c(mut_short, json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[n]][["structure"]])
                        }
                      } 
                      if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
                        for (n in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]])) {
                          mut <- c(mut, json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[n]][["displayStructure"]])
                          mut_short <- c(mut_short, json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[n]][["structure"]])
                        }
                      }
                      mut <- paste(mut, collapse = ";")
                      mut_short <- paste(mut_short, collapse = ";")
                      df[[paste0(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]], "_mut")]] <- mut
                      df[[paste0(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]], "_mut_short")]] <- mut_short
                    }
                  } else { 
                    df[[json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]]]] <- "Insufficient coverage"
                }
              }
            } 
    } 
  }
  
  # Then join mutations per drug category
  tmp <- as_tibble(df)
  
  tmp <- tmp %>% 
    unite("NS34A", c(glecaprevir_mut, grazoprevir_mut, paritaprevir_mut, voxilaprevir_mut), sep = ";", na.rm = TRUE) %>% 
    unite("NS5A", c(daclatasvir_mut, elbasvir_mut, ledipasvir_mut, ombitasvir_mut, pibrentasvir_mut, velpatasvir_mut), sep = ";", na.rm = TRUE) %>% 
    unite("NS5B", c(dasabuvir_mut, sofosbuvir_mut), sep = ";", na.rm = TRUE) %>% 
    unite("NS34A_short", c(glecaprevir_mut_short, grazoprevir_mut_short, paritaprevir_mut_short, voxilaprevir_mut_short), sep = ";", na.rm = TRUE) %>% 
    unite("NS5A_short", c(daclatasvir_mut_short, elbasvir_mut_short, ledipasvir_mut_short, ombitasvir_mut_short, pibrentasvir_mut_short, velpatasvir_mut_short), sep = ";", na.rm = TRUE) %>% 
    unite("NS5B_short", c(dasabuvir_mut_short, sofosbuvir_mut_short), sep = ";", na.rm = TRUE) 
  
  # Gjøre om innholdet i cellene til en vector, deretter fjerne dupliater i vektoren
  try(df$NS34A <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS34A, ";")), ","))), "\\+"))), collapse = ";"))
  try(df$NS5A <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5A, ";")), ","))), "\\+"))), collapse = ";"))
  try(df$NS5B <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5B, ";")), ","))), "\\+"))), collapse = ";"))
  try(df$NS34A_short <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS34A_short, ";")), ","))), "\\+"))), collapse = ";"))
  try(df$NS5A_short <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5A_short, ";")), ","))), "\\+"))), collapse = ";"))
  try(df$NS5B_short <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5B_short, ";")), ","))), "\\+"))), collapse = ";"))
  
  df <- as_tibble(df)
  
  # Merge with final data structure
  df_final <- bind_rows(df_final, df)

  reference <- unlist(strsplit(basename(json_file), "\\."))[[2]]
  major_minor <- unlist(strsplit(basename(json_file), "\\."))[[3]]

  write_tsv(df_final, file = paste0(sample, ".", reference, ".", major_minor, "_GLUE_report.tsv"))