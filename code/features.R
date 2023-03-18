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

home <- "C:/Projects/RCourse/Masterclass/CAGE"

# -------------------------------------------------
# dependencies - input files used by this script
# vpca .. PCA of covariance, scale = FALSE
# cpca .. PCA of correlations, scale = TRUE
#
FILE_VAL <- file.path(home, "data/cache/validation.rds")
FILE_TRN <- file.path(home, "data/cache/training.rds")
FILE_PAT <- file.path(home, "data/cache/patients.rds")
FILE_VST <- file.path(home, "data/cache/train_vpca_scores.rds")
FILE_VSV <- file.path(home, "data/cache/valid_vpca_scores.rds")
FILE_CST <- file.path(home, "data/cache/train_cpca_scores.rds")
FILE_CSV <- file.path(home, "data/cache/valid_cpca_scores.rds")
FILE_CFN <- file.path(home, "code/calculation_functions.R")
# -------------------------------------------------
# targets - output files created by this script
#
FILE_FLS <- file.path(home, "data/cache/feature_loss.rds")
FILE_LGF <- file.path(home, "data/cache/features_log.txt")
# --------------------------------------------------
# Divert warning messages to a log file
#
logFile <- file(FILE_LGF, open = "wt")
sink(logFile, type = "message")

source(FILE_CFN)

# --------------------------------------------
# Read data on the 1000 probes
# add diagnosis from patientDF and class: 0=Benign 1=Cancer
#
readRDS(FILE_PAT) %>%
  mutate(class = ifelse(diagnosis == "Cancer", 1, 0)) %>% 
  select(id, diagnosis, class) -> classDF
  
readRDS(FILE_VAL) %>%
  left_join(classDF, by = "id")  -> validDF 

readRDS(FILE_TRN) %>%
  left_join(classDF, by = "id") -> trainDF 

readRDS(FILE_VST) %>%
  left_join(classDF, by = "id") -> vScoreTrainDF 

readRDS(FILE_VSV) %>%
  left_join(classDF, by = "id") -> vScoreValidDF 

readRDS(FILE_CST) %>%
  left_join(classDF, by = "id") -> cScoreTrainDF 

readRDS(FILE_CSV) %>%
  left_join(classDF, by = "id") -> cScoreValidDF 

# ------------------------------------------------
# Probe names
#
names(trainDF)[2:1001] -> probeNames

# ------------------------------------------------
# Loss from univariate logistic regressions 
#
probeUniDF <- univariate_logistic(trainDF, 
                                  "class", 
                                  probeNames) 


# ------------------------------------------------------------------
# Logistic Regression using 1 to 50 selected predictors
#
probeUniDF %>%
  arrange(loss) %>%
  multiple_logistic(trainDF, 
                    validDF, 
                    "class") -> probeSelectedDF

  
# ------------------------------------------------------------
# PC names
#
pcaNames <- glue("PC{1:375}")

# ------------------------------------------------
# Loss from univariate logistic regressions
#
vpcaUniDF <- univariate_logistic(vScoreTrainDF, 
                                 "class", 
                                 pcaNames)


# ------------------------------------------------------------------
# Multivariate Logistic using 1 to 50 PCs selected by univariate loss
#
vpcaUniDF %>%
  arrange(loss) %>%
  multiple_logistic(vScoreTrainDF, 
                    vScoreValidDF, 
                    "class") -> vpcaSelectedDF

# ------------------------------------------------------------------
# Multivariate Logistic using first 1 to 50 PCs
#
vpcaUniDF %>%
   multiple_logistic(vScoreTrainDF, 
                     vScoreValidDF, 
                     "class") -> vpcaOrderedDF

# ------------------------------------------------
# Loss from univariate logistic regressions
#
cpcaUniDF <- univariate_logistic(cScoreTrainDF, 
                                 "class", 
                                 pcaNames)

# ------------------------------------------------------------------
# Multivariate Logistic using 1 to 50 selected PCs
#
cpcaUniDF %>%
  arrange(loss) %>%
  multiple_logistic(cScoreTrainDF, 
                    cScoreValidDF, 
                    "class") -> cpcaSelectedDF

# ------------------------------------------------------------------
# Multivariate Logistic using first 1 to 50 PCs
#
cpcaUniDF %>%
  multiple_logistic(cScoreTrainDF, 
                    cScoreValidDF, 
                    "class") -> cpcaOrderedDF

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
  saveRDS(FILE_FLS)

# -----------------------------------------------
# Close the log file
#
sink(type = "message") 
close(logFile)