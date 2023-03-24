# -----------------------------------------------------------
# Project CAGE
# Principal Component Analysis of 1000 Probes
#
# Date: 21 March 2023
#
# -----------------------------------------------------------
library(tidyverse)
library(fs)

# -------------------------------------------------
# data folder
#
cache    <- "C:/Projects/RCourse/Masterclass/CAGE/data/cache"
# -------------------------------------------------
# dependencies - input files used by this script
#
valRDS <- path(cache, "validation.rds")
trnRDS <- path(cache, "training.rds")
patRDS <- path(cache, "patients.rds")

# -------------------------------------------------
# targets - output files created by this script
# vpca .. PCA of covariance, scale = FALSE
# cpca .. PCA of correlations, scale = TRUE
#
vpcRDS <- path(cache, "train_vpca.rds")
vstRDS <- path(cache, "train_vpca_scores.rds")
vsvRDS <- path(cache, "valid_vpca_scores.rds")
cpcRDS <- path(cache, "train_cpca.rds")
cstRDS <- path(cache, "train_cpca_scores.rds")
csvRDS <- path(cache, "valid_cpca_scores.rds")

# --------------------------------------------------
# Divert warning messages to a log file
#
lf <- file(path(cache, "pca_log.txt"), open = "wt")
sink(lf, type = "message")

# --------------------------------------------
# Read data on 1000 probes
#
validDF   <- readRDS(valRDS)
trainDF   <- readRDS(trnRDS)
patientDF <- readRDS(patRDS)

# --------------------------------------------
# vpca = covariance PCA of the training data
#
trainDF %>%
  select(-id) %>%
  as.matrix() %>% 
  prcomp(scale = FALSE) %>%
  saveRDS(vpcRDS) 

# --------------------------------------------
# vpca scores of the training data
#
readRDS(vpcRDS) %>%
  predict() %>%
  as_tibble() %>%
  mutate(id = trainDF$id) %>%
  relocate(id) %>%
  saveRDS(vstRDS) 

# --------------------------------------------
# vpca scores of the validation data
#
readRDS(vpcRDS) %>%
  predict(newdata = validDF) %>%
  as_tibble() %>%
  mutate(id = validDF$id) %>%
  relocate(id) %>%
  saveRDS(vsvRDS) 

# --------------------------------------------
# cpca = correlation PCA of the training data
#
trainDF %>%
  select(-id) %>%
  as.matrix() %>% 
  prcomp(scale = TRUE) %>%
  saveRDS(cpcRDS) 

# --------------------------------------------
# cpca scores of the training data
#
readRDS(cpcRDS) %>%
  predict() %>%
  as_tibble() %>%
  mutate(id = trainDF$id) %>%
  relocate(id) %>%
  saveRDS(cstRDS) 

# --------------------------------------------
# cpca scores of the validation data
#
readRDS(cpcRDS) %>%
  predict(newdata = validDF) %>%
  as_tibble() %>%
  mutate(id = validDF$id) %>%
  relocate(id) %>%
  saveRDS(csvRDS) 

# -----------------------------------------------
# Close the log file
#
sink(type = "message") 
close(lf)