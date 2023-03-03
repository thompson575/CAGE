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
library(glue)
library(magrittr)

home <- "C:/Projects/RCourse/Masterclass/CAGE"

source( file.path(home, "code/CAGE_functions.R") )

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
names(trainDF)[2:1001] -> probeNames

# ------------------------------------------------
# Loss from univariate logistic regressions 
#
probeUniDF <- univariate_logistic(trainDF, "class", probeNames) 

# ----------------------------------------------------------
# Plot of the 5 best probes in training and validation data
#
probeUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(trainDF) +
    ggtitle("Training Data: 5 most predictive probes")

probeUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(validDF) +
  ggtitle("Validation Data: 5 most predictive probes")

# ------------------------------------------------------------------
# Logistic Regression using 1 to 50 selected predictors
#
probeUniDF %>%
  arrange( loss ) %>%
  multiple_logistic(trainDF, validDF, "class") -> probeSelectedDF

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
plot_in_out_loss(probeSelectedDF) +
  labs( title = "Loss from Selected Probes",
        y     = "Loss",
        x     = "Number of Probes")
  
# =======================================================
# Selected Unscaled Principal Components
#
# PC names
#
pcaNames <- glue("PC{1:375}")

# ------------------------------------------------
# Loss from univariate logistic regressions
#
pcaUniDF <- univariate_logistic(scoreTrainDF, "class", pcaNames)

# ----------------------------------------------------------
# Plot of the 5 best PCs in training and validation data
#
pcaUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(scoreTrainDF) +
  ggtitle("Training Data: 5 most predictive PCs")

pcaUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(scoreValidDF) +
  ggtitle("Validation Data: 5 most predictive PCs")

# ------------------------------------------------------------------
# Multivariate Logistic using 1 to 50 PCs selected by univariate loss
#
pcaUniDF %>%
  arrange( loss ) %>%
  multiple_logistic(scoreTrainDF, scoreValidDF, "class") -> pcaSelectedDF

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
plot_in_out_loss(pcaSelectedDF) +
  labs( title = "Loss from Selected PCs",
        y     = "Loss",
        x     = "Number of PCs")


# =======================================================
# First M Unscaled Principal Components
#

# ------------------------------------------------------------------
# Multivariate Logistic using first 1 to 50 PCs
#
pcaUniDF %>%
   multiple_logistic(scoreTrainDF, scoreValidDF, "class") -> pcaOrderedDF

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
plot_in_out_loss(pcaOrderedDF) +
  labs( title = "Loss from First PCs",
        y     = "Loss",
        x     = "Number of PCs")

# =======================================================
# Selected Scaled Principal Components
#

# ------------------------------------------------
# Loss from univariate logistic regressions
#
pcaScaledUniDF <- univariate_logistic(scoreTrainScaledDF, "class", pcaNames)

# ----------------------------------------------------------
# Plot of the 5 best PCs in training and validation data
#
pcaScaledUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(scoreTrainScaledDF) +
  ggtitle("Training Data: 5 most predictive scaled PCs")

pcaScaledUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(scoreValidScaledDF) +
  ggtitle("Validation Data: 5 most predictive scaled PCs")

# ------------------------------------------------------------------
# Multivariate Logistic using 1 to 50 selected PCs
#
pcaScaledUniDF %>%
  arrange( loss ) %>%
  multiple_logistic(scoreTrainScaledDF, 
                    scoreValidScaledDF, 
                    "class")              -> pcaScaledSelectedDF

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
plot_in_out_loss(pcaScaledSelectedDF) +
  labs( title = "Loss from Selected Scaled PCs",
        y     = "Loss",
        x     = "Number of PCs")

# =======================================================
# First M Scaled Principal Components
#

# ------------------------------------------------------------------
# Multivariate Logistic using first 1 to 50 PCs
#
pcaScaledUniDF %>%
   multiple_logistic(scoreTrainScaledDF, 
                     scoreValidScaledDF, 
                     "class")             -> pcaScaledOrderedDF

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
plot_in_out_loss(pcaScaledOrderedDF) +
  labs( title = "Loss from First Scaled PCs",
        y     = "Loss",
        x     = "Number of PCs")

# ---------------------------------------------
# Plot in-sample loss for the five methods M=1-50
#
list( probe  = probeSelectedDF,
      PCsel  = pcaSelectedDF,
      PCord  = pcaOrderedDF,
      sPCsel = pcaScaledSelectedDF,
      sPCord = pcaScaledOrderedDF)  %>%
  plot_method_comparison(inloss) +
  ggtitle("In-sample loss for the 5 feature selection methods")

# ---------------------------------------------
# Plot in-sample loss for the five methods M=1-8
#
list( probe  = probeSelectedDF     %>% slice(1:8),
      PCsel  = pcaSelectedDF       %>% slice(1:8),
      PCord  = pcaOrderedDF        %>% slice(1:8),
      sPCsel = pcaScaledSelectedDF %>% slice(1:8),
      sPCord = pcaScaledOrderedDF  %>% slice(1:8))  %>%
  plot_method_comparison(inloss) +
  ggtitle("In-sample loss for the 5 feature selection methods")

# ---------------------------------------------
# Plot out-of-sample loss for the five methods M=1-50
#
list( probe  = probeSelectedDF,
      PCsel  = pcaSelectedDF,
      PCord  = pcaOrderedDF,
      sPCsel = pcaScaledSelectedDF,
      sPCord = pcaScaledOrderedDF)  %>%
  plot_method_comparison(outloss) +
  ggtitle("Out-of-sample loss for the 5 feature selection methods")


# ---------------------------------------------
# Plot out-of-sample loss for the five methods M=1-8
#
list( probe  = probeSelectedDF     %>% slice(1:8),
      PCsel  = pcaSelectedDF       %>% slice(1:8),
      PCord  = pcaOrderedDF        %>% slice(1:8),
      sPCsel = pcaScaledSelectedDF %>% slice(1:8),
      sPCord = pcaScaledOrderedDF  %>% slice(1:8))  %>%
  plot_method_comparison(outloss) +
  ggtitle("Out-of-sample loss for the 5 feature selection methods")
# ---------------------------------------------
# Save results to dataStore
#
saveRDS( list( probeUniDF          = probeUniDF, 
               probeSelectedDF     = probeSelectedDF,
               pcaUniDF            = pcaUniDF, 
               pcaSelectedDF       = pcaSelectedDF, 
               pcaOrderedDF        = pcaOrderedDF, 
               pcaScaledUniDF      = pcaScaledUniDF, 
               pcaScaledSelectedDF = pcaScaledSelectedDF, 
               pcaScaledOrderedDF  = pcaScaledOrderedDF),
         file.path(home, "data/dataStore/classification_loss.rds"))
