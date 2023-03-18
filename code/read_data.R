# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Read the Series Matrix File and create 
#   patients.rds   ... patient characteristics
#   expression.rds ... expression data 
# Sample 1000 probes to create an exploratory subset
#   training.rds   ... AEGIS-1 data
#   validation.rds ... AEGIS-2 data
#
# Date: 19 Feb 2023
#
library(tidyverse)

home <- "C:/Projects/RCourse/Masterclass/CAGE"

# -------------------------------------------------
# dependencies - input files used by this script
#
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210271/matrix/GSE210271_series_matrix.txt.gz"

# -------------------------------------------------
# targets - output files created by this script
#
file_LSM <- file.path(home, "data/rawData/GSE210271_series_matrix.txt.gz")
file_EXN <- file.path(home, "data/cache/expression.rds")
file_PAT <- file.path(home, "data/cache/patients.rds")
file_VAL <- file.path(home, "data/cache/validation.rds")
file_TRN <- file.path(home, "data/cache/training.rds")
file_LGF <- file.path(home, "data/cache/read_data_log.txt")

# --------------------------------------------------
# Divert warning messages to a log file
#
logFile <- file(file_LGF, open = "wt")
sink( logFile, type = "message")

# --------------------------------------------------
# Download the series file from GEO. Save as file_LSM
#
if( !file.exists(file_LSM) ) 
  download.file(url, file_LSM)

# --------------------------------------------------
# Read the file as lines of text for exploration
#

# lines <- readLines( file_LSM )

# substr(lines[1:15], 1, 30)

# --------------------------------------------------
# line 14 contains the patient identifiers
#
read.table(file_LSM,  sep = '\t', 
           header = FALSE, skip = 13, nrows = 1) %>%
  { strsplit(.$V2, " ")[[1]] } -> patientId

# --------------------------------------------------
# lines 41-45 contain the patient characteristics
#
read.table(file_LSM,  sep = '\t', 
                     header = FALSE, skip = 41, nrows = 5) %>%
  as_tibble() %>%
  select( -1) %>%
  mutate( across(everything(), ~ str_replace(., "Sex: ", ""))) %>%
  mutate( across(everything(), ~ str_replace(., "age: ", ""))) %>%
  mutate( across(everything(), ~ str_replace(., "smoking status: ", ""))) %>%
  mutate( across(everything(), ~ str_replace(., "fev1 % predicted: ", ""))) %>%
  mutate( across(everything(), ~ str_replace(., "cancer status: ", ""))) %>%
  mutate( var = c("sex", "age", "smoking", "fev1", "diagnosis")) %>%
  pivot_longer(-var, names_to = "col", values_to = "data") %>%
  pivot_wider( names_from = var, values_from = data) %>%
  mutate( id = patientId ) %>%
  mutate( study = c(rep("AEGIS1", 375), rep("AEGIS2", 130))) %>%
  select(id, study, sex, age, smoking, diagnosis) %>%
  mutate( age = as.numeric(age)) %>%
  saveRDS(file_PAT) 

# --------------------------------------------------
# lines 70-21754 contain gene expressions
# 21685 probes for 505 patients
#
read.table(file_LSM,  sep = '\t', 
           header = TRUE, skip = 70, nrows = 21685) %>%
  as_tibble() %>%
  saveRDS(file_EXN) 

# -----------------------------------------------
# Create a sample of 1000 probes 
#
set.seed(7382)
sample(1:21685, size = 1000, replace = FALSE) %>%
  sort() -> rows

# -----------------------------------------------
# Use AEGIS-1 (columns 2-376) as the training data
# transpose so that rows = patients, cols = probes
#
readRDS(file_EXN) %>%
  .[rows, 1:376] %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  saveRDS(file_TRN) 

# -----------------------------------------------
# USE AEGIS-2 (columns 377-506) as the validation data
# transpose so that rows = patients, cols = probes
#
readRDS(file_EXN) %>%
  .[rows, c(1, 377:506)]  %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  saveRDS(file_VAL) 

# -----------------------------------------------
# Close the log file
#
sink( type = "message" ) 
close(logFile)