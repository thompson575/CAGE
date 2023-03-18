# -----------------------------------------------------------
# Project CAGE
# Principal Component Analysis of 1000 Probes
#
# Date: 19 Feb 2023
#
# -----------------------------------------------------------
library(tidyverse)

home <- "C:/Projects/RCourse/Masterclass/CAGE/"

# -------------------------------------------------
# dependencies - input files used by this script
#
file_VAL <- file.path(home, "data/cache/validation.rds")
file_TRN <- file.path(home, "data/cache/training.rds")
file_PAT <- file.path(home, "data/cache/patients.rds")

# -------------------------------------------------
# targets - output files created by this script
# vpca .. PCA of covariance, scale = FALSE
# cpca .. PCA of correlations, scale = TRUE
#
file_VPC <- file.path(home, "data/cache/train_vpca.rds")
file_VST <- file.path(home, "data/cache/train_vpca_scores.rds")
file_VSV <- file.path(home, "data/cache/valid_vpca_scores.rds")
file_CPC <- file.path(home, "data/cache/train_cpca.rds")
file_CST <- file.path(home, "data/cache/train_cpca_scores.rds")
file_CSV <- file.path(home, "data/cache/valid_cpca_scores.rds")
file_LGF <- file.path(home, "data/cache/pca_log.txt")

# --------------------------------------------------
# Divert warning messages to a log file
#
logFile <- file(file_LGF, open = "wt")
sink(logFile, type = "message")

# --------------------------------------------
# Read data on 1000 probes
#
validDF   <- readRDS(file_VAL)
trainDF   <- readRDS(file_TRN)
patientDF <- readRDS(file_PAT)

# --------------------------------------------
# covariance pca of the training data
#
trainDF %>%
  select(-id) %>%
  as.matrix() %>% 
  prcomp(scale = FALSE) %>%
  saveRDS(file_VPC) 

# --------------------------------------------
# vpca scores of the training data
#
readRDS(file_VPC) %>%
  predict() %>%
  as_tibble() %>%
  mutate(id = trainDF$id) %>%
  relocate(id) %>%
  saveRDS(file_VST) 

# --------------------------------------------
# vpca scores of the validation data
#
readRDS(file_VPC) %>%
  predict(newdata = validDF) %>%
  as_tibble() %>%
  mutate(id = validDF$id) %>%
  relocate(id) %>%
  saveRDS(file_VSV) 

# --------------------------------------------
# correlation pca of the training data
#
trainDF %>%
  select(-id) %>%
  as.matrix() %>% 
  prcomp(scale = TRUE) %>%
  saveRDS(file_CPC) 

# --------------------------------------------
# pca scores of the training data
#
readRDS(file_CPC) %>%
  predict() %>%
  as_tibble() %>%
  mutate(id = trainDF$id) %>%
  relocate(id) %>%
  saveRDS(file_CST) 

# --------------------------------------------
# pca scores of the validation data
#
readRDS(file_CPC) %>%
  predict(newdata = validDF) %>%
  as_tibble() %>%
  mutate(id = validDF$id) %>%
  relocate(id) %>%
  saveRDS(file_CSV) 

# -----------------------------------------------
# Close the log file
#
sink(type = "message") 
close(logFile)