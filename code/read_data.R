# -----------------------------------------------------------
# Project CAGE: 
#
# Read the Series Matrix File and create 
#   patients.rds   ... patient characteristics
#   expression.rds ... expression data 
# Sample 1000 probes to create an exploratory dataset
#   training.rds   ... AEGIS-1 data
#   validation.rds ... AEGIS-2 data
#
# Date: 21 March 2023
#
library(tidyverse)
library(fs)

# -------------------------------------------------
# data folders
#
cache    <- "C:/Projects/RCourse/Masterclass/CAGE/data/cache"
rawData  <- "C:/Projects/RCourse/Masterclass/CAGE/data/rawData"

# -------------------------------------------------
# dependencies - input files used by this script
#
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210271/matrix/GSE210271_series_matrix.txt.gz"

# -------------------------------------------------
# targets - output files created by this script
#
serRDS <- path(rawData, "GSE210271_series_matrix.txt.gz")
exnRDS <- path(cache,   "expression.rds")
patRDS <- path(cache,   "patients.rds")
valRDS <- path(cache,   "validation.rds")
trnRDS <- path(cache,   "training.rds")

# --------------------------------------------------
# Divert warning messages to a log file
#
lf <- file(path(cache,   "read_log.txt"), open = "wt")
sink(lf, type = "message")

# --------------------------------------------------
# Download the series file from GEO. Save in rawData
#
if(!file.exists(serRDS) ) 
  download.file(url, serRDS)

# --------------------------------------------------
# Read the file as lines of text for exploration
#
# lines <- readLines(serRDS )
# substr(lines[1:15], 1, 30)

# --------------------------------------------------
# line 14 contains the patient identifiers
#
read.table(serRDS,  
           sep    = '\t', 
           header = FALSE, 
           skip   = 13, 
           nrows  = 1) %>%
  { strsplit(.$V2, " ")[[1]] } -> patientId

# --------------------------------------------------
# lines 41-45 contain the patient characteristics
#
read.table(serRDS,  
           sep    = '\t', 
           header = FALSE, 
           skip   = 41, 
           nrows  = 5) %>%
  as_tibble() %>%
  select(-1) %>%
  mutate(across(everything(), ~ str_replace(., "Sex: ", ""))) %>%
  mutate(across(everything(), ~ str_replace(., "age: ", ""))) %>%
  mutate(across(everything(), ~ str_replace(., "smoking status: ", ""))) %>%
  mutate(across(everything(), ~ str_replace(., "fev1 % predicted: ", ""))) %>%
  mutate(across(everything(), ~ str_replace(., "cancer status: ", ""))) %>%
  mutate(var = c("sex", "age", "smoking", "fev1", "diagnosis")) %>%
  pivot_longer(-var, names_to = "col", values_to = "data") %>%
  pivot_wider(names_from = var, values_from = data) %>%
  mutate(id = patientId ) %>%
  mutate(study = c(rep("AEGIS1", 375), rep("AEGIS2", 130))) %>%
  select(id, study, sex, age, smoking, diagnosis) %>%
  mutate(age = as.numeric(age)) %>%
  saveRDS(patRDS) 

# --------------------------------------------------
# lines 70-21754 contain gene expressions
# 21685 probes for 505 patients
#
read.table(serRDS,  
           sep    = '\t', 
           header = TRUE, 
           skip   = 70, 
           nrows  = 21685) %>%
  as_tibble() %>%
  saveRDS(exnRDS) 

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
readRDS(exnRDS) %>%
  .[rows, 1:376] %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  saveRDS(trnRDS) 

# -----------------------------------------------
# USE AEGIS-2 (columns 377-506) as the validation data
# transpose so that rows = patients, cols = probes
#
readRDS(exnRDS) %>%
  .[rows, c(1, 377:506)]  %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  saveRDS(valRDS) 

# -----------------------------------------------
# Close the log file
#
sink(type = "message" ) 
close(lf)