# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Read the raw data and create 
#   subjDF (subjects.rds) ... subject characteristics
#   geDF (expression.rds) ... expression data
# Sample 1000 probes to create an exploratory subset
#   trainDF (training.rds)   ... AEGIS-1 data
#   validDF (validation.rds) ... AEGIS-2 data
#
# Date: 19 Feb 2023
#
library(tidyverse)
library(magrittr)

home <- "C:/Projects/RCourse/Masterclass/CAGE"
file <- "data/rawData/GSE210271_series_matrix.txt.gz"

# --------------------------------------------------
# Read the file as lines of text
#
lines <- readLines( file.path(home, file) )

# --------------------------------------------------
# line 13 contains the subject identifiers
#
read.table(file.path(home, file),  sep = '\t', 
           header = FALSE, skip = 13, nrows = 1) %>%
  { strsplit(.$V2, " ")[[1]] } -> id

# --------------------------------------------------
# lines 41-45 contain the subject characteristics
#
read.table(file.path(home, file),  sep = '\t', 
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
  mutate( id = id) %>%
  mutate( study = c(rep("AEGIS1", 375), rep("AEGIS2", 130))) %>%
  select(id, study, sex, age, smoking, diagnosis) %>%
  mutate( age = as.numeric(age)) %>%
  print() %T>%
  saveRDS( file.path(home, "data/rData/subjects.rds") ) -> subjDF

# --------------------------------------------------
# lines 70-21754 contain gene expressions
# 21685 probes for 505 subjects
#
read.table(file.path(home, file),  sep = '\t', 
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
  print() %T>%
  saveRDS( file.path(home, "data/rData/training.rds") )-> trainDF

# -----------------------------------------------
# USE AEGIS-2 (columns 377-506) as the validation data
# transpose so that rows = subjects, cols = probes
#
geDF[rows, c(1, 377:506)]  %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  print() %T>%
  saveRDS( file.path(home, "data/rData/validation.rds") ) -> validDF

