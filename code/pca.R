# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Principal Component Analysis of the exploratory subset
#
# Date: 19 Feb 2023
#
library(tidyverse)
library(glue)
library(magrittr)

home <- "C:/Projects/RCourse/Masterclass/CAGE/"

source( file.path(home, "code/CAGE_functions.R") )

# --------------------------------------------
# Read exploratory subset of probes
#
validDF <- readRDS( file.path(home, "data/rData/validation.rds"))
trainDF <- readRDS( file.path(home, "data/rData/training.rds"))
subjDF  <- readRDS( file.path(home, "data/rData/subjects.rds"))

# =======================================================
# UNSCALED PRCOMP .. PCA of Covariance
#

# --------------------------------------------
# pca of the training data
#
trainDF %>%
  select( -id) %>%
  as.matrix() %>% 
  prcomp() %T>%
  saveRDS( file.path(home, "data/rData/subset_pca.rds")) -> pca

# --------------------------------------------
# pca scores of the training data
#
predict(pca) %>%
  as_tibble() %>%
  mutate( id = trainDF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, 
                    "data/rData/subset_train_pca_scores.rds")) -> scoreTrainDF

# --------------------------------------------
# pca scores of the validation data
#
predict(pca, newdata=validDF) %>%
  as_tibble() %>%
  mutate( id = validDF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, 
                    "data/rData/subset_valid_pca_scores.rds")) -> scoreValidDF


# ---------------------------------------------
# Principal Components Stdev (root eigenvalues)
#
plot_eigenvalues(pca$sdev) +
  ggtitle("Standard deviatons of the Principal Components")
# ---------------------------------------------
# Principal Components Percent Variance
#
plot_eigenvalues(pca$sdev[1:50], pct = TRUE) +
  ggtitle("Percent Explained by the Principal Components")

# ---------------------------------------------
# Illustrative plot of first 5 PCs by diagnosis
# training data
#
glue( "PC{1:5}" ) %>%
  boxplot_features( scoreTrainDF %>%
                      left_join(subjDF, by = "id")) +
  labs(title = "Training data: PCA Scores by diagnosis",
       y     = "Score",
       x     = "Principal Component")

# ---------------------------------------------
# Illustrative plot of first 5 PCs by diagnosis
# validation data
#
glue( "PC{1:5}" ) %>%
  boxplot_features( scoreValidDF %>%
                      left_join(subjDF, by = "id")) +
  labs(title = "Validation data: PCA Scores by diagnosis",
       y     = "Score",
       x     = "Principal Component")


# =======================================================
# SCALED PRCOMP .. PCA of Correlations
#

# --------------------------------------------
# pca of the training data
#
trainDF %>%
  select( -id) %>%
  as.matrix() %>% 
  prcomp(scale=TRUE) %T>%
  saveRDS( file.path(home, 
                     "data/rData/subset_scaled_pca.rds")) -> pcaScaled

# --------------------------------------------
# pca scores of the training data
#
predict(pcaScaled) %>%
  as_tibble() %>%
  mutate( id = trainDF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, 
                    "data/rData/subset_train_scaled_pca_scores.rds")) -> scoreScaledTrainDF

# --------------------------------------------
# pca scores of the validation data
#
predict(pcaScaled, newdata=validDF) %>%
  as_tibble() %>%
  mutate( id = validDF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, 
                    "data/rData/subset_valid_scaled_pca_scores.rds")) -> scoreScaledValidDF

# ---------------------------------------------
# Scaled Principal Components Stdev (root eigenvalues)
#
plot_eigenvalues(pcaScaled$sdev) +
  ggtitle("Standard deviatons of the Scaled Principal Components")
# ---------------------------------------------
# Principal Components Percent Variance
#
plot_eigenvalues(pcaScaled$sdev[1:50], pct = TRUE) +
  ggtitle("Percent Explained by the Scaled Principal Components")

# ---------------------------------------------
# Illustrative plot of first 5 Scaled PCs by diagnosis
# training data
#
glue( "PC{1:5}" ) %>%
  boxplot_features( scoreScaledTrainDF %>%
                      left_join(subjDF, by = "id")) +
  labs(title = "Training data: Scaled PCA Scores by diagnosis",
       y     = "Score",
       x     = "Principal Component")

# ---------------------------------------------
# Illustrative plot of first 5 Scaled PCs by diagnosis
# validation data
#
glue( "PC{1:5}" ) %>%
  boxplot_features( scoreScaledValidDF %>%
                      left_join(subjDF, by = "id")) +
  labs(title = "Validation data: Scaled PCA Scores by diagnosis",
       y     = "Score",
       x     = "Principal Component")

