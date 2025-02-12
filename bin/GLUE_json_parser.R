#!/usr/bin/env Rscript

library(tidyverse)
library(jsonlite)

# Define json parser function
parse_json_files <- function(major_minor) {
  
  json_files <- list.files(pattern = paste0(major_minor, ".nodup.json$"),
                           full.names = TRUE)
  
# Create final data file
df_final <- tibble(
  "Sample" = character(),
  "Reference" = character(),
  "Major_minor" = character(),
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

# Only run if the length of json_files is greater than 0
if (length(json_files) > 0) {
for (x in 1:length(json_files)) {
  # Remove any old objects if reading previous json file failed
  try(rm(json))

  # Need to check if there are any lines starting with "DEBUG" in the json files and remove these

  json <- readLines(json_files[x]) # Read the json file line by line

    # Check if there are any lines starting with "DEBUG"
    if (length(grep("^DEBUG", json)) > 0) {
      # Remove these lines
      json <- json[-grep("^DEBUG", json)]

      # And then check if the first line is empty and remove that
      if (json[1] == "") {
        json <- json[-1]
      }
    }

  # Only parse the json object if the second element contains the string "phdrReport"
  if (grepl("phdrReport", json[2])) {

  # Try to parse json object. Could fail if bam file was not OK for Glue
  try(json <- parse_json(json))

  # Check that the json object exists
  if (exists("json")) {

    # Only read the proper json GLUE reports (i.e. that there was a good sequence)
    if (names(json) == "phdrReport") {
    # Sample name
    sample <- unlist(strsplit(basename(json_files[x]), "\\."))[[1]]
    reference <- unlist(strsplit(basename(json_files[x]), "\\."))[[2]]
    major_minor <- unlist(strsplit(basename(json_files[x]), "\\."))[[3]]

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
    try(rm(df_tmp))
    df_tmp <- as.data.frame(matrix(nrow = 1, ncol = 50))
    colnames(df_tmp) <- c("Sample",
                      "Reference",
                      "Major_minor",
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

    df_tmp$Sample <- sample
    df_tmp$Reference <- reference
    df_tmp$Major_minor <- major_minor
    df_tmp$GLUE_genotype <- genotype
    df_tmp$GLUE_subtype <- subtype
    df_tmp$`HCV project version` <- projectVersion
    df_tmp$`GLUE engine version` <- engineVersion
    df_tmp$`PHE drug resistance extension version` <- extensionVersion

    # Dette er underlisten for Drug Scores. Lengden av denne angir hvor mange drugs som er funnet.
    # Under drugScores så er det en ny liste for hver drug category
    if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]]) > 0) {
      for (i in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]])) {
        # Så er det en ny liste innenfor hver drug score igjen med drug for hver kategori. Denne heter drugAssessemnts
        for (k in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]])) {
                  # Hvis det er sufficient coverage (denne evaluerer til TRUE):
                  if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["sufficientCoverage"]]) {
                    if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]] == "No resistance") {
                      df_tmp[[json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]]]]  <- "No resistance"
                    } else {
                      # Skrive inn resistance informasjonen for druget
                      df_tmp[[json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]]]]  <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]

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
                      df_tmp[[paste0(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]], "_mut")]] <- mut
                      df_tmp[[paste0(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]], "_mut_short")]] <- mut_short
                    }
                  } else {
                    df_tmp[[json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]]]] <- "Insufficient coverage"
                }
              }
            }
    }
  }
  }
}
  # Then join mutations per drug category
  # if not NA in df_tmp$Sample
  if (exists("df_tmp")) {
  tmp <- as_tibble(df_tmp)

  tmp <- tmp %>%
    unite("NS34A", c(glecaprevir_mut, grazoprevir_mut, paritaprevir_mut, voxilaprevir_mut), sep = ";", na.rm = TRUE) %>%
    unite("NS5A", c(daclatasvir_mut, elbasvir_mut, ledipasvir_mut, ombitasvir_mut, pibrentasvir_mut, velpatasvir_mut), sep = ";", na.rm = TRUE) %>%
    unite("NS5B", c(dasabuvir_mut, sofosbuvir_mut), sep = ";", na.rm = TRUE) %>%
    unite("NS34A_short", c(glecaprevir_mut_short, grazoprevir_mut_short, paritaprevir_mut_short, voxilaprevir_mut_short), sep = ";", na.rm = TRUE) %>%
    unite("NS5A_short", c(daclatasvir_mut_short, elbasvir_mut_short, ledipasvir_mut_short, ombitasvir_mut_short, pibrentasvir_mut_short, velpatasvir_mut_short), sep = ";", na.rm = TRUE) %>%
    unite("NS5B_short", c(dasabuvir_mut_short, sofosbuvir_mut_short), sep = ";", na.rm = TRUE)

  # Gjøre om innholdet i cellene til en vector, deretter fjerne dupliater i vektoren
  try(df_tmp$NS34A <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS34A, ";")), ","))), "\\+"))), collapse = ";"))
  try(df_tmp$NS5A <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5A, ";")), ","))), "\\+"))), collapse = ";"))
  try(df_tmp$NS5B <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5B, ";")), ","))), "\\+"))), collapse = ";"))
  try(df_tmp$NS34A_short <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS34A_short, ";")), ","))), "\\+"))), collapse = ";"))
  try(df_tmp$NS5A_short <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5A_short, ";")), ","))), "\\+"))), collapse = ";"))
  try(df_tmp$NS5B_short <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5B_short, ";")), ","))), "\\+"))), collapse = ";"))

  df_tmp <- as_tibble(df_tmp)

  # Merge with final data structure
  df_final <- bind_rows(df_final, df_tmp)
  }
}
}

return(df_final)
}

df_major <- parse_json_files("major")
write_tsv(df_major, file = paste0("GLUE_collected_report_major.tsv"))

df_minor <- parse_json_files("minor")
write_tsv(df_minor, file = paste0("GLUE_collected_report_minor.tsv"))


