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
FILE_VAL <- file.path(home, "data/cache/validation.rds")
FILE_TRN <- file.path(home, "data/cache/training.rds")
FILE_PAT <- file.path(home, "data/cache/patients.rds")

# -------------------------------------------------
# targets - output files created by this script
# vpca .. PCA of covariance, scale = FALSE
# cpca .. PCA of correlations, scale = TRUE
#
FILE_VPC <- file.path(home, "data/cache/train_vpca.rds")
FILE_VST <- file.path(home, "data/cache/train_vpca_scores.rds")
FILE_VSV <- file.path(home, "data/cache/valid_vpca_scores.rds")
FILE_CPC <- file.path(home, "data/cache/train_cpca.rds")
FILE_CST <- file.path(home, "data/cache/train_cpca_scores.rds")
FILE_CSV <- file.path(home, "data/cache/valid_cpca_scores.rds")
FILE_LGF <- file.path(home, "data/cache/pca_log.txt")

# --------------------------------------------------
# Divert warning messages to a log file
#
logFile <- file(FILE_LGF, open = "wt")
sink(logFile, type = "message")

# --------------------------------------------
# Read data on 1000 probes
#
validDF   <- readRDS(FILE_VAL)
trainDF   <- readRDS(FILE_TRN)
patientDF <- readRDS(FILE_PAT)

# --------------------------------------------
# covariance pca of the training data
#
trainDF %>%
  select(-id) %>%
  as.matrix() %>% 
  prcomp(scale = FALSE) %>%
  saveRDS(FILE_VPC) 

# --------------------------------------------
# vpca scores of the training data
#
readRDS(FILE_VPC) %>%
  predict() %>%
  as_tibble() %>%
  mutate(id = trainDF$id) %>%
  relocate(id) %>%
  saveRDS(FILE_VST) 

# --------------------------------------------
# vpca scores of the validation data
#
readRDS(FILE_VPC) %>%
  predict(newdata = validDF) %>%
  as_tibble() %>%
  mutate(id = validDF$id) %>%
  relocate(id) %>%
  saveRDS(FILE_VSV) 

# --------------------------------------------
# correlation pca of the training data
#
trainDF %>%
  select(-id) %>%
  as.matrix() %>% 
  prcomp(scale = TRUE) %>%
  saveRDS(FILE_CPC) 

# --------------------------------------------
# pca scores of the training data
#
readRDS(FILE_CPC) %>%
  predict() %>%
  as_tibble() %>%
  mutate(id = trainDF$id) %>%
  relocate(id) %>%
  saveRDS(FILE_CST) 

# --------------------------------------------
# pca scores of the validation data
#
readRDS(FILE_CPC) %>%
  predict(newdata = validDF) %>%
  as_tibble() %>%
  mutate(id = validDF$id) %>%
  relocate(id) %>%
  saveRDS(FILE_CSV) 

# -----------------------------------------------
# Close the log file
#
sink(type = "message") 
close(logFile)