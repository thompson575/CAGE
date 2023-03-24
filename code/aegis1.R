# -----------------------------------------------------------
# Project CAGE:
#
# Analysis of full data from AEGIS-1 
# using AEGIS-2 for validation
#
# Date: 21 March 2023
#
library(tidyverse)
library(glue)
library(fs)
library(magrittr)

# -------------------------------------------------
# data & code folders
#
cache    <- "C:/Projects/RCourse/Masterclass/CAGE/data/cache"
code     <- "C:/Projects/RCourse/Masterclass/CAGE/code"
# -------------------------------------------------
# dependencies - input files used by this script
#
exnRDS <- path(cache, "expression.rds")
patRDS <- path(cache, "patients.rds")

# -------------------------------------------------
# targets - output files created by this script
# vpca .. PCA of covariance,   scale = FALSE
# cpca .. PCA of correlations, scale = TRUE
#
ae1RDS <- path(cache, "expression_aegis1.rds")
ae2RDS <- path(cache, "expression_aegis2.rds")
vpvRDS <- path(cache, "aegis1_vpca.rds")
vs1RDS <- path(cache, "aegis1_vpca_scores.rds")
vs2RDS <- path(cache, "aegis2_vpca_scores.rds")
cpcRDS <- path(cache, "aegis1_cpca.rds")
cs1RDS <- path(cache, "aegis1_cpca_scores.rds")
cs2RDS <- path(cache, "aegis2_cpca_scores.rds")
losRDS <- path(cache, "aegis1_loss.rds")

# --------------------------------------------------
# Divert warning messages to a log file
#
lf <- file(path(cache, "aegis1_log.txt"), open = "wt")
sink(lf, type = "message")

# --------------------------------------------------
# Source the computation functions
#
source(path(code, "calculation_functions.R"))

# --------------------------------------------------
# Extract class for each patient
#
readRDS(patRDS) %>%
  mutate(class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
  select(id, class) -> classDF

# --------------------------------------------
# Read AEGIS-1 & AEGIS2 for all 21685 probes
#
readRDS(exnRDS) %>%
  { .[, 1:376] } %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %T>%
  # --- save expressions for aegis1
  saveRDS(ae1RDS) %>%
  # --- add class info
  left_join(classDF, by = "id") %>%
  relocate(id, class) -> aegis1DF

readRDS(exnRDS) %>%
  { .[, c(1, 377:506)] } %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %T>%
  # --- save expressions for aegis2
  saveRDS(ae2RDS) %>%
  # --- add class info
  left_join(classDF, by = "id") %>%
  relocate(id, class) -> aegis2DF

# --------------------------------------------
# Covariance PCA (scale=FALSE)
# vpca (unscaled) of the aegis1 data
#
aegis1DF %>%
  select(-id, -class) %>%
  as.matrix() %>% 
  prcomp() %T>%
  saveRDS(vpvRDS) -> vpca

# --------------------------------------------
# vpca (unscaled) scores of the AEGIS-1 data
#
predict(vpca) %>%
  as_tibble() %>%
  mutate(id = aegis1DF$id) %>%
  relocate(id) %T>%
  # --- save vpca scores
  saveRDS(vs1RDS) %>%
  # --- add class info
  left_join(classDF, by = "id") -> vpcaScore1DF

# --------------------------------------------
# vpca (unscaled) scores of the AEGIS-2 data
#
predict(vpca, newdata=aegis2DF) %>%
  as_tibble() %>%
  mutate(id = aegis2DF$id) %>%
  relocate(id) %T>%
  # --- save vpca scores
  saveRDS(vs2RDS) %>%
  # --- add class info
  left_join(classDF, by = "id") -> vpcaScore2DF

# --------------------------------------------
# Correlation PCA (scale=TRUE)
# cpca (scaled) of the aegis1 data
#
aegis1DF %>%
  select(-id, -class) %>%
  as.matrix() %>% 
  prcomp(scale=TRUE) %T>%
  saveRDS(cpcRDS) -> cpca
# --------------------------------------------
# cpca (scaled) scores of the AEGIS-1 data
#
predict(cpca) %>%
  as_tibble() %>%
  mutate(id = aegis1DF$id) %>%
  relocate(id) %T>%
  # --- save cpcs scores
  saveRDS(cs1RDS) %>%
  # --- add class info
  left_join(classDF, by = "id") -> cpcaScore1DF

# --------------------------------------------
# cpca (scaled) scores of the AEGIS-2 data
#
predict(cpca, newdata=aegis2DF) %>%
  as_tibble() %>%
  mutate(id = aegis2DF$id) %>%
  relocate(id) %T>%
  # --- save cpcs scores
  saveRDS(cs2RDS) %>%
  # --- add class info
  left_join(classDF, by = "id") -> cpcaScore2DF

# ------------------------------------------------
# Probe names
#
names(aegis1DF)[c(-1, -2)] -> probeNames

# ------------------------------------------------
# Loss from univariate logistic regressions
#
univariate_logistic(aegis1DF, 
                    "class", 
                    probeNames) -> probeUniDF

# ------------------------------------------------
# Multiple Logistic using 1 to 50 selected predictors
#
probeUniDF %>%
  arrange(loss) %>%
  multiple_logistic(aegis1DF, 
                    aegis2DF, 
                    "class")    -> probeSelectedDF

# ------------------------------------------------
# PC names
#
pcaNames <- glue("PC{1:375}")

# ------------------------------------------------
# Covariance PCA (scale=FALSE)
# Loss from univariate logistic regressions
#
univariate_logistic(vpcaScore1DF, 
                    "class", 
                    pcaNames)   -> vpcaUniDF

# ------------------------------------------------
# Multiple Logistic using 1 to 50 selected PCs
#
vpcaUniDF %>%
  arrange(loss) %>%
  multiple_logistic(vpcaScore1DF, 
                    vpcaScore2DF, 
                    "class")     -> vpcaSelectedDF

# ------------------------------------------------------------------
# Multiple Logistic using 1 to 50 ordered PCs
#
vpcaUniDF %>%
  multiple_logistic(vpcaScore1DF, 
                    vpcaScore2DF, 
                    "class")      -> vpcaOrderedDF

# ------------------------------------------------
# Correlation PCA (scale=TRUE)
# Loss from univariate logistic regressions
#
univariate_logistic(cpcaScore1DF, 
                    "class", 
                    pcaNames)     -> cpcaUniDF

# ------------------------------------------------------------------
# Multiple Logistic using 1 to 50 selected PCs
#
cpcaUniDF %>%
  arrange(loss) %>%
  multiple_logistic(cpcaScore1DF, 
                    cpcaScore2DF, 
                    "class")      -> cpcaSelectedDF

# ------------------------------------------------------------------
# Multiple Logistic using 1 to 50 ordered PCs
#
cpcaUniDF %>%
  multiple_logistic(cpcaScore1DF, 
                    cpcaScore2DF, 
                    "class")      -> cpcaOrderedDF

# ---------------------------------------------
# Save results to cache
#
list(probeUniDF      = probeUniDF, 
     probeSelectedDF = probeSelectedDF,
     vpcaUniDF       = vpcaUniDF, 
     vpcaSelectedDF  = vpcaSelectedDF, 
     vpcaOrderedDF   = vpcaOrderedDF, 
     cpcaUniDF       = cpcaUniDF, 
     cpcaSelectedDF  = cpcaSelectedDF, 
     cpcaOrderedDF   = cpcaOrderedDF) %>%
   saveRDS(losRDS)

# -----------------------------------------------
# Close the log file
#
sink(type = "message") 
close(lf)