# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Analysis of full data for AEGIS-1
#
# Date: 19 Feb 2023
#
library(tidyverse)
library(glue)
library(magrittr)

home <- "C:/Projects/RCourse/Masterclass/CAGE/"

source("code/CAGE_functions.R")

# --------------------------------------------
# Read AEGIS-1 & AEGIS2 for all 21685 probes
#
subjDF  <- readRDS( file.path(home, "data/rData/subjects.rds"))

readRDS( file.path(home, "data/rData/expression.rds") ) %>%
  { .[, 1:376] } %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  print() %T>%
  saveRDS( file.path(home, "data/rData/expression_aegis1.rds") ) %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") %>%
  relocate( id, diagnosis, class ) -> aegis1DF

readRDS( file.path(home, "data/rData/expression.rds") ) %>%
  { .[, c(1, 377:506)] } %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  print() %T>%
  saveRDS( file.path(home, "data/rData/expression_aegis1.rds") ) %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") %>%
  relocate( id, diagnosis, class ) -> aegis2DF

# --------------------------------------------
# pca of the aegis1 data
#
aegis1DF %>%
  select( -id, -diagnosis, -class) %>%
  as.matrix() %>% 
  prcomp() %T>%
  saveRDS( file.path(home, "data/rData/aegis1_pca.rds")) -> pca

# --------------------------------------------
# pca scores of the AEGIS-1 data
#
predict(pca) %>%
  as_tibble() %>%
  mutate( id = aegis1DF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, "data/rData/aegis1_pca_scores.rds")) %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> scoreAegis1DF

# --------------------------------------------
# pca scores of the AEGIS-2 data
#
predict(pca, newdata=aegis2DF) %>%
  as_tibble() %>%
  mutate( id = aegis2DF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, "data/rData/aegis2_pca_scores.rds")) %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> scoreAegis2DF

# --------------------------------------------
# Scaled pca of the aegis1 data
#
aegis1DF %>%
  select( -id, -diagnosis, -class) %>%
  as.matrix() %>% 
  prcomp( scale=TRUE) %T>%
  saveRDS( file.path(home, "data/rData/aegis1_scaled_pca.rds")) -> pcaScaled

# --------------------------------------------
# pca scores of the AEGIS-1 data
#
predict(pcaScaled) %>%
  as_tibble() %>%
  mutate( id = aegis1DF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, "data/rData/aegis1_scaled_pca_scores.rds")) %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> scoreScaledAegis1DF

# --------------------------------------------
# pca scores of the AEGIS-2 data
#
predict(pcaScaled, newdata=aegis2DF) %>%
  as_tibble() %>%
  mutate( id = aegis2DF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, "data/rData/aegis2_scaled_pca_scores.rds")) %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> scoreScaledAegis2DF

# =======================================================
# Selected probes
#
# Probe names
#
names(aegis1DF)[c(-1, -2, -3)] -> probeNames

# ------------------------------------------------
# Loss from univariate logistic regressions
#
probeUniDF <- univariate_logistic(aegis1DF, "class", probeNames) 

# ----------------------------------------------------------
# Plot of the 5 best probes in AEGIS-1 and AEGIS-2 data
#
probeUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(aegis1DF) +
  ggtitle("AEGIS-1: 5 most predictive probes")

probeUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(aegis2DF) +
  ggtitle("AEGIS-2: 5 most predictive probes")

# ------------------------------------------------------------------
# Multiple Logistic using 1 to 50 selected predictors
#
probeUniDF %>%
  arrange( loss ) %>%
  multiple_logistic(aegis1DF, aegis2DF, "class") -> probeSelectedDF

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
plot_in_out_loss(probeSelectedDF) +
  labs( title = "AEGIS-1: Loss from Selected Probes",
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
pcaUniDF <- univariate_logistic(scoreAegis1DF, "class", pcaNames)

# ----------------------------------------------------------
# Plot of the 5 best PCs in AEGIS-1 and AEGIS-2 data
#
pcaUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(scoreAegis1DF) +
  ggtitle("AEGIS-1 Data: 5 most predictive PCs")

pcaUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(scoreAegis2DF) +
  ggtitle("AEGIS-2 Data: 5 most predictive PCs")

# ------------------------------------------------------------------
# Multiple Logistic using 1 to 50 selected PCs
#
pcaUniDF %>%
  arrange( loss ) %>%
  multiple_logistic(scoreAegis1DF, scoreAegis2DF, "class") -> pcaSelectedDF

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
plot_in_out_loss(pcaSelectedDF) +
  labs( title = "Loss from Selected PCs",
        y     = "Loss",
        x     = "Number of PCs")

# =======================================================
# Initial Unscaled Principal Components
#

# ------------------------------------------------------------------
# Multiple Logistic using 1 to 50 ordered PCs
#
pcaUniDF %>%
  multiple_logistic(scoreAegis1DF, scoreAegis2DF, "class") -> pcaOrderedDF

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
pcaScaledUniDF <- univariate_logistic(scoreScaledAegis1DF, "class", pcaNames)


# ----------------------------------------------------------
# Plot of the 5 best PCs in AEGIS-1 and AEGIS-2 data
#
pcaScaledUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(scoreScaledAegis1DF) +
  ggtitle("AEGIS-1 Data: 5 most predictive scaled PCs")

pcaScaledUniDF %>%
  arrange(loss) %>%
  slice( 1:5 ) %>%
  pull(x) %>%
  boxplot_features(scoreScaledAegis2DF) +
  ggtitle("AEGIS-2 Data: 5 most predictive scaled PCs")

# ------------------------------------------------------------------
# Multiple Logistic using 1 to 50 selected PCs
#
pcaScaledUniDF %>%
  arrange( loss ) %>%
  multiple_logistic(scoreScaledAegis1DF, 
                    scoreScaledAegis2DF, 
                    "class")              -> pcaScaledSelectedDF


# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
plot_in_out_loss(pcaScaledSelectedDF) +
  labs( title = "Loss from Selected Scaled PCs",
        y     = "Loss",
        x     = "Number of PCs")

# =======================================================
# Initial Scaled Principal Components
#

# ------------------------------------------------------------------
# Multiple Logistic using 1 to 50 ordered PCs
#
pcaScaledUniDF %>%
  multiple_logistic(scoreScaledAegis1DF, 
                    scoreScaledAegis2DF, 
                    "class")             -> pcaScaledOrderedDF

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
plot_in_out_loss(pcaScaledOrderedDF) +
  labs( title = "Loss from First Scaled PCs",
        y     = "Loss",
        x     = "Number of PCs")

# ---------------------------------------------
# Plot in-sample loss for the five methods
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
         file.path(home, "data/dataStore/aegis1_classification_loss.rds"))
