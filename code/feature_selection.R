# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Feature selection comparison
# Selected probes vs First few PCs vs Selected PCs
#
# Date: 19 Feb 2023
#       22 Feb 2023 .. added the scaled PCA results
#
library(tidyverse)
library(magrittr)

home <- "C:/Projects/RCourse/Masterclass/CAGE"

source("code/CAGE_functions.R")

# --------------------------------------------
# Read data on the subset of 1000 probes
# add diagnosis from subjDF and class: 0=Benign 1=Cancer
#
readRDS( file.path(home, "data/rData/subjects.rds")) %>%
  mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) -> subjDF
  
readRDS( file.path(home, "data/rData/validation.rds")) %>%
  left_join( subjDF %>%
               select(id, diagnosis, class), by = "id")  -> validDF 
readRDS( file.path(home, "data/rData/training.rds")) %>%
  left_join( subjDF %>%
               select(id, diagnosis, class), by = "id") -> trainDF 

readRDS(file.path(home, "data/rData/subset_train_pca_scores.rds")) %>%
  left_join( subjDF %>%
               select(id, diagnosis, class), by = "id") -> scoreTrainDF 

readRDS(file.path(home, "data/rData/subset_valid_pca_scores.rds")) %>%
  left_join( subjDF %>%
               select(id, diagnosis, class), by = "id") -> scoreValidDF 

readRDS(file.path(home, "data/rData/subset_train_scaled_pca_scores.rds")) %>%
  left_join( subjDF %>%
               select(id, diagnosis, class), by = "id") -> scoreTrainScaledDF 

readRDS(file.path(home, "data/rData/subset_valid_scaled_pca_scores.rds")) %>%
  left_join( subjDF %>%
               select(id, diagnosis, class), by = "id") -> scoreValidScaledDF 

# =======================================================
# Selected probes
#
# Probe names
#
names(trainDF)[2:1001] -> pbNames

# ------------------------------------------------
# Loss from univariate logistic regressions
#
pbUniDF <- uni_logistic(trainDF, "class", pbNames) 

# ----------------------------------------------------------
# Plot of the 5 best probes in training and validation data
#
pbUniDF %>% arrange(loss) %>% print() %>% slice(1:5) %>% pull(x) -> topX

trainDF %>%
  pivot_longer( starts_with("ENSG"), names_to="probe", values_to="expression") %>%
  filter( probe %in% topX ) %>%
  ggplot( aes(x=probe, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Training Data: 5 most predictive probes",
        x = "") +
  theme( legend.position = c(0.85, 0.85))

validDF %>%
  pivot_longer( starts_with("ENSG"), names_to="probe", values_to="expression") %>%
  filter( probe %in% topX ) %>%
  ggplot( aes(x=probe, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Validation Data: 5 most predictive probes",
        x = "") +
  theme( legend.position = c(0.85, 0.85))

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 selected predictors
#
pbSelDF <- sel_logistic(trainDF, validDF, "class", pbUniDF, maxSel = 50)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pbSelDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "Loss from Selected Probes",
        y     = "Loss",
        x     = "Number of Probes") +
  geom_text( x=35, y=0.58, label="In-sample", colour="blue") +
  geom_text( x=35, y=0.90, label="Out-of-sample", colour="red") 
  
# =======================================================
# Selected Unscaled Principal Components
#
# PC names
#
pcNames <- paste0("PC", 1:375)

# ------------------------------------------------
# Loss from univariate logistic regressions
#
pcUniDF <- uni_logistic(scoreTrainDF, "class", pcNames)

# ----------------------------------------------------------
# Plot of the 5 best PCs in training and validation data
#
pcUniDF %>% arrange(loss) %>% print() %>% slice(1:5) %>% pull(x) -> topX

scoreTrainDF %>%
  pivot_longer( starts_with("PC"), names_to="pc", values_to="expression") %>%
  filter( pc %in% topX ) %>%
  ggplot( aes(x=pc, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Training Data: 5 most predictive PCs",
        x = "") +
  theme( legend.position = c(0.15, 0.85))

scoreValidDF %>%
  pivot_longer( starts_with("PC"), names_to="pc", values_to="expression") %>%
  filter( pc %in% topX ) %>%
  ggplot( aes(x=pc, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Validation Data: 5 most predictive PCs",
        x = "") +
  theme( legend.position = c(0.15, 0.85))

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 selected PCs
#
pcSelDF <- sel_logistic(scoreTrainDF, scoreValidDF, "class", pcUniDF, maxSel = 50)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pcSelDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "Loss from Selected PCs",
        y     = "Loss",
        x     = "Number of PCs") +
  geom_text( x=35, y=0.58, label="In-sample", colour="blue") +
  geom_text( x=35, y=2.00, label="Out-of-sample", colour="red")


# =======================================================
# Initial Unscaled Principal Components
#

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 ordered PCs
#
pcOrdDF <- sel_logistic(scoreTrainDF, scoreValidDF, "class", pcUniDF, 
                    maxSel = 50, sort = FALSE)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pcOrdDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "Loss from Ordered PCs",
        y     = "Loss",
        x     = "Number of PCs") +
  geom_text( x=35, y=0.58, label="In-sample", colour="blue") +
  geom_text( x=35, y=1.00, label="Out-of-sample", colour="red")

# =======================================================
# Selected Scaled Principal Components
#
# PC names
#
pcNames <- paste0("PC", 1:375)

# ------------------------------------------------
# Loss from univariate logistic regressions
#
pcScdUniDF <- uni_logistic(scoreTrainScaledDF, "class", pcNames)

# ----------------------------------------------------------
# Plot of the 5 best PCs in training and validation data
#
pcScdUniDF %>% arrange(loss) %>% print() %>% slice(1:5) %>% pull(x) -> topX

scoreTrainScaledDF %>%
  pivot_longer( starts_with("PC"), names_to="pc", values_to="expression") %>%
  filter( pc %in% topX ) %>%
  ggplot( aes(x=pc, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Training Data: 5 most predictive Scaled PCs",
        x = "") +
  theme( legend.position = c(0.15, 0.85))

scoreValidScaledDF %>%
  pivot_longer( starts_with("PC"), names_to="pc", values_to="expression") %>%
  filter( pc %in% topX ) %>%
  ggplot( aes(x=pc, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Validation Data: 5 most predictive Scaled PCs",
        x = "") +
  theme( legend.position = c(0.15, 0.85))

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 selected PCs
#
pcScdSelDF <- sel_logistic(scoreTrainScaledDF, scoreValidScaledDF, 
                           "class", pcScdUniDF, maxSel = 50)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pcScdSelDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "Loss from Selected Scaled PCs",
        y     = "Loss",
        x     = "Number of PCs") +
  geom_text( x=35, y=0.58, label="In-sample", colour="blue") +
  geom_text( x=35, y=2.00, label="Out-of-sample", colour="red")


# =======================================================
# Initial Scaled Principal Components
#

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 ordered PCs
#
pcScdOrdDF <- sel_logistic(scoreTrainScaledDF, scoreValidScaledDF, 
                           "class", pcScdUniDF, maxSel = 50, sort = FALSE)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pcScdOrdDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "Loss from Ordered Scaled PCs",
        y     = "Loss",
        x     = "Number of PCs") +
  geom_text( x=35, y=0.58, label="In-sample", colour="blue") +
  geom_text( x=35, y=1.00, label="Out-of-sample", colour="red")

# ---------------------------------------------
# Plot in-sample loss for the five methods
#
pbSelDF %>%
  rename( probe = inloss) %>%
  left_join( pcSelDF %>% rename( pcSel  = inloss), by = "n") %>%
  left_join( pcOrdDF %>% rename( pcOrd = inloss), by = "n") %>%
  left_join( pcScdSelDF %>% rename( pcScdSel  = inloss), by = "n") %>%
  left_join( pcScdOrdDF %>% rename( pcScdOrd = inloss), by = "n") %>%
  pivot_longer(c("probe", "pcSel", "pcOrd", "pcScdSel", "pcScdOrd"), 
               names_to = "method", values_to = "loss") %>%
  mutate( source = str_replace(method, "loss", ""))  %>%
  ggplot( aes(x = n, y = loss, colour=method)) +
  geom_line( size=1.2) +
  labs( title = "In-sample loss for the 5 feature selection methods",
        y     = "Loss",
        x     = "Number of predictors") +
  theme( legend.position = c(0.15, 0.25))

pbSelDF %>%
  rename( probe = inloss) %>%
  left_join( pcSelDF %>% rename( pcSel  = inloss), by = "n") %>%
  left_join( pcOrdDF %>% rename( pcOrd = inloss), by = "n") %>%
  left_join( pcScdSelDF %>% rename( pcScdSel  = inloss), by = "n") %>%
  left_join( pcScdOrdDF %>% rename( pcScdOrd = inloss), by = "n") %>%
  pivot_longer(c("probe", "pcSel", "pcOrd", "pcScdSel", "pcScdOrd"), 
               names_to = "method", values_to = "loss") %>%
  mutate( source = str_replace(method, "loss", ""))  %>%
  filter( n <= 8 ) %>%
  ggplot( aes(x = n, y = loss, colour=method)) +
  geom_line( size=1.2) +
  labs( title = "In-sample loss for the 5 feature selection methods",
        y     = "Loss",
        x     = "Number of predictors") +
  theme( legend.position = c(0.15, 0.25))

pbSelDF %>%
  rename( probe = outloss) %>%
  left_join( pcSelDF %>% rename( pcSel = outloss), by = "n") %>%
  left_join( pcOrdDF %>% rename( pcOrd = outloss), by = "n") %>%
  left_join( pcScdSelDF %>% rename( pcScdSel = outloss), by = "n") %>%
  left_join( pcScdOrdDF %>% rename( pcScdOrd = outloss), by = "n") %>%
  pivot_longer(c("probe", "pcSel", "pcOrd", "pcScdSel", "pcScdOrd"), 
               names_to = "method", values_to = "loss") %>%
  mutate( source = str_replace(method, "loss", ""))  %>%
  ggplot( aes(x = n, y = loss, colour=source)) +
  geom_line( size=1.2) +
  labs( title = "Out-of-sample loss for the 5 feature selection methods",
        y     = "Loss",
        x     = "Number of predictors") +
  theme( legend.position = c(0.15, 0.85))

pbSelDF %>%
  rename( probe = outloss) %>%
  left_join( pcSelDF %>% rename( pcSel = outloss), by = "n") %>%
  left_join( pcOrdDF %>% rename( pcOrd = outloss), by = "n") %>%
  left_join( pcScdSelDF %>% rename( pcScdSel = outloss), by = "n") %>%
  left_join( pcScdOrdDF %>% rename( pcScdOrd = outloss), by = "n") %>%
  pivot_longer(c("probe", "pcSel", "pcOrd", "pcScdSel", "pcScdOrd"), 
               names_to = "method", values_to = "loss") %>%
  mutate( source = str_replace(method, "loss", ""))  %>%
  filter( n <= 8 ) %>%
  ggplot( aes(x = n, y = loss, colour=source)) +
  geom_line( size=1.2) +
  labs( title = "Out-of-sample loss for the 5 feature selection methods",
        y     = "Loss",
        x     = "Number of predictors") +
  theme( legend.position = c(0.15, 0.85))

# ---------------------------------------------
# Save results to dataStore
#
saveRDS( list( pcUniDF = pcUniDF, pcSelDF = pcSelDF, 
               pcOrdDF = pcOrdDF, 
               pbUniDF = pbUniDF, pbSelDF = pbSelDF,
               pcScdUniDF = pcScdUniDF, pcScdSelDF = pcScdSelDF, 
               pcScdOrdDF = pcScdOrdDF),
         file.path(home, "data/dataStore/classification_loss.rds"))
