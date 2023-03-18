# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Analysis of full data from AEGIS-1
#
# Date: 19 Feb 2023
#
library(tidyverse)
library(glue)
library(magrittr)

home <- "C:/Projects/RCourse/Masterclass/CAGE"

# -------------------------------------------------
# dependencies - input files used by this script
#
file_EXN <- file.path(home, "data/cache/expression.rds")
file_PAT <- file.path(home, "data/cache/patients.rds")
file_CFN <- file.path(home, "code/calculation_functions.R")

# -------------------------------------------------
# targets - output files created by this script
# vpca .. PCA of covariance,   scale = Ffile_LS1E
# cpca .. PCA of correlations, scale = TRUE
#
file_EA1 <- file.path(home, "data/cache/expression_aegis1.rds")
file_EA2 <- file.path(home, "data/cache/expression_aegis2.rds")
file_VPC <- file.path(home, "data/cache/aegis1_vpca.rds")
file_VS1 <- file.path(home, "data/cache/aegis1_vpca_scores.rds")
file_VS2 <- file.path(home, "data/cache/aegis2_vpca_scores.rds")
file_CPC <- file.path(home, "data/cache/aegis1_cpca.rds")
file_CS1 <- file.path(home, "data/cache/aegis1_cpca_scores.rds")
file_CS2 <- file.path(home, "data/cache/aegis2_cpca_scores.rds")
file_LS1 <- file.path(home, "data/cache/aegis1_loss.rds")
file_LGF <- file.path(home, "data/cache/aegis1_log.txt")

# --------------------------------------------------
# Divert warning messages to a log file
#
logFile <- file(file_LGF, open = "wt")
sink(logFile, type = "message")

# --------------------------------------------------
# Source the calculation functions
#
source(file_CFN)

# --------------------------------------------------
# Extract class for each patient
#
readRDS(file_PAT) %>%
  mutate(class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
  select(id, class) -> classDF

# --------------------------------------------
# Read AEGIS-1 & AEGIS2 for all 21685 probes
#
readRDS(file_EXN) %>%
  { .[, 1:376] } %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %T>%
  # --- save expressions for aegis1
  saveRDS(file_EA1) %>%
  # --- add class info
  left_join(classDF, by = "id") %>%
  relocate(id, class) -> aegis1DF

readRDS(file_EXN) %>%
  { .[, c(1, 377:506)] } %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %T>%
  # --- save expressions for aegis2
  saveRDS(file_EA2) %>%
  # --- add class info
  left_join(classDF, by = "id") %>%
  relocate(id, class) -> aegis2DF

# --------------------------------------------
# vpca (unscaled) of the aegis1 data
#
aegis1DF %>%
  select(-id, -class) %>%
  as.matrix() %>% 
  prcomp() %T>%
  saveRDS(file_VPC) -> vpca

# --------------------------------------------
# vpca (unscaled) scores of the AEGIS-1 data
#
predict(vpca) %>%
  as_tibble() %>%
  mutate(id = aegis1DF$id) %>%
  relocate(id) %T>%
  # --- save vpca scores
  saveRDS(file_VS1) %>%
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
  saveRDS(file_VS2) %>%
  # --- add class info
  left_join(classDF, by = "id") -> vpcaScore2DF

# --------------------------------------------
# cpca (scaled) of the aegis1 data
#
aegis1DF %>%
  select(-id, -class) %>%
  as.matrix() %>% 
  prcomp(scale=TRUE) %T>%
  saveRDS(file_CPC) -> cpca
# --------------------------------------------
# cpca (scaled) scores of the AEGIS-1 data
#
predict(cpca) %>%
  as_tibble() %>%
  mutate(id = aegis1DF$id) %>%
  relocate(id) %T>%
  # --- save cpcs scores
  saveRDS(file_CS1) %>%
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
  saveRDS(file_CS2) %>%
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
   saveRDS(file_LS1)

# -----------------------------------------------
# Close the log file
#
sink(type = "message") 
close(logFile)