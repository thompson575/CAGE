# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Principal Component Analysis of the exploratory subset
#
# Date: 19 Feb 2023
#
library(tidyverse)
library(magrittr)

home <- "C:/Projects/RCourse/Masterclass/CAGE/"

# --------------------------------------------
# Read exploratory subset of probes
#
validDF <- readRDS( file.path(home, "data/rData/validation.rds"))
trainDF <- readRDS( file.path(home, "data/rData/training.rds"))
subjDF  <- readRDS( file.path(home, "data/rData/subjects.rds"))

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
  saveRDS(file.path(home, "data/rData/subset_train_pca_scores.rds")) -> scoreTrainDF

# --------------------------------------------
# pca scores of the validation data
#
predict(pca, newdata=validDF) %>%
  as_tibble() %>%
  mutate( id = validDF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, "data/rData/subset_valid_pca_scores.rds")) -> scoreValidDF

# ---------------------------------------------
# Principal Components Stdev (root eigenvalues)
#
tibble( pc = 1:375,
        sd = pca$sdev) %>%
  ggplot( aes(x = pc, y = sd, xend = pc, yend = 0)) +
  geom_segment() +
  labs( title = "Standard deviatons of the Principal Components",
        y     = "Standard deviation",
        x     = "Principal Component Number")

# ---------------------------------------------
# Principal Components Percent Variance
#
tibble( pc = 1:375,
        sd = pca$sdev) %>%
  mutate( pct = 100 * sd * sd / sum(sd * sd)) %>%
  ggplot( aes(x = pc, y = pct, xend = pc, yend = 0)) +
  geom_segment() +
  labs( title = "Standard deviatons of the Principal Components",
        y     = "Percent Variance",
        x     = "Principal Component Number")


# ---------------------------------------------
# Illustrative plot of first 5 PCs by diagnosis
# training data
#
scoreTrainDF %>%
  select( id, PC1:PC5) %>%
  pivot_longer(starts_with("PC"), names_to = "PC", 
               values_to = "score") %>%
  left_join( subjDF %>% 
               select(id, diagnosis), by = "id") %>%
  ggplot( aes(x = PC, y = score, fill = diagnosis)) +
  geom_boxplot() +
  labs(title = "Training data: PCA Scores 1 to 5 by diagnosis",
       y     = "Score",
       x     = "Principal Component") +
  theme( legend.position = c(.85, .9))

# ---------------------------------------------
# Illustrative plot of first 5 PCs by diagnosis
# validation data
#
scoreValidDF %>%
  select( id, PC1:PC5) %>%
  pivot_longer(starts_with("PC"), names_to = "PC", 
               values_to = "score") %>%
  left_join( subjDF %>% 
               select(id, diagnosis), by = "id") %>%
  ggplot( aes(x = PC, y = score, fill = diagnosis)) +
  geom_boxplot() +
  labs(title = "Validation data: PCA Scores 1 to 5 by diagnosis",
       y     = "Score",
       x     = "Principal Component") +
  theme( legend.position = c(.85, .9))
