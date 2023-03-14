# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Read the raw data and create 
#   subjects.rds   ... subject characteristics
#   expression.rds ... expression data (geDF)
# Sample 1000 probes to create an exploratory subset
#   training.rds   ... AEGIS-1 data
#   validation.rds ... AEGIS-2 data
#
# Date: 19 Feb 2023
#
library(tidyverse)
library(magrittr)

URL  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210271/matrix/GSE210271_series_matrix.txt.gz"

home <- "C:/Projects/RCourse/Masterclass/CAGE"
file <- "data/rawData/GSE210271_series_matrix.txt.gz"

# --------------------------------------------------
# Download the series file from GEO. Save as localCopy
#
localCopy <- file.path(home, file)

download.file(URL, localCopy)

# --------------------------------------------------
# Read the file as lines of text for exploration
#
lines <- readLines( localCopy )

substr(lines[1:15], 1, 30)
# --------------------------------------------------
# line 14 contains the subject identifiers
#
read.table(localCopy,  sep = '\t', 
           header = FALSE, skip = 13, nrows = 1) %>%
  { strsplit(.$V2, " ")[[1]] } -> subjectId

# --------------------------------------------------
# lines 41-45 contain the subject characteristics
#
read.table(localCopy,  sep = '\t', 
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
  mutate( id = subjectId ) %>%
  mutate( study = c(rep("AEGIS1", 375), rep("AEGIS2", 130))) %>%
  select(id, study, sex, age, smoking, diagnosis) %>%
  mutate( age = as.numeric(age)) %>%
  print() %>%
  saveRDS( file.path(home, "data/rData/subjects.rds") ) 

# --------------------------------------------------
# lines 70-21754 contain gene expressions
# 21685 probes for 505 subjects
#
read.table(localCopy,  sep = '\t', 
           header = TRUE, skip = 70, nrows = 21685) %>%
  as_tibble() %>%
  print() %T>%
  saveRDS( file.path(home, "data/rData/expression.rds") ) -> geDF

# -----------------------------------------------
# Create a sample of 1000 probes 
#
set.seed(7382)
sample(1:21685, size = 1000, replace = FALSE) %>%
  sort() -> rows

# -----------------------------------------------
# Use AEGIS-1 (columns 2-376) as the training data
# transpose so that rows = subjects, cols = probes
#
geDF[rows, 1:376] %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  print() %>%
  saveRDS( file.path(home, "data/rData/training.rds") ) 

# -----------------------------------------------
# USE AEGIS-2 (columns 377-506) as the validation data
# transpose so that rows = subjects, cols = probes
#
geDF[rows, c(1, 377:506)]  %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  print() %>%
  saveRDS( file.path(home, "data/rData/validation.rds") ) 

