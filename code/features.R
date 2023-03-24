# -----------------------------------------------------------
# Project CAGE:
#
# Feature selection comparison of 5 methods
# Selected probes vs First few vpca PCs vs Selected vpca PCs
#                 vs First few cpca PCs vs Selected cpca PCs
# Date: 21 March 2023
#
library(tidyverse)
library(glue)
library(fs)

# -------------------------------------------------
# data & code folders
#
cache    <- "C:/Projects/RCourse/Masterclass/CAGE/data/cache"
code     <- "C:/Projects/RCourse/Masterclass/CAGE/code"
# -------------------------------------------------
# dependencies - input files used by this script
# vpca .. PCA of covariance,   scale = FALSE
# cpca .. PCA of correlations, scale = TRUE
#
valRDS <- path(cache, "validation.rds")
trnRDS <- path(cache, "training.rds")
patRDS <- path(cache, "patients.rds")
vstRDS <- path(cache, "train_vpca_scores.rds")
vsvRDS <- path(cache, "valid_vpca_scores.rds")
cstRDS <- path(cache, "train_cpca_scores.rds")
csvRDS <- path(cache, "valid_cpca_scores.rds")
# -------------------------------------------------
# targets - output files created by this script
#
losRDS <- path(cache, "feature_loss.rds")
# --------------------------------------------------
# Divert warning messages to a log file
#
lf <- file(path(cache, "features_log.txt"), open = "wt")
sink(lf, type = "message")

# --------------------------------------------------
# Read Computation Functions
#
source(path(code, "calculation_functions.R"))

# --------------------------------------------
# Read data on the 1000 probes
# add class: 0=Benign 1=Cancer
#
readRDS(patRDS) %>%
  mutate(class = ifelse(diagnosis == "Cancer", 1, 0)) %>% 
  select(id, diagnosis, class) -> classDF
  
readRDS(valRDS) %>%
  left_join(classDF, by = "id")  -> validDF 

readRDS(trnRDS) %>%
  left_join(classDF, by = "id") -> trainDF 

readRDS(vstRDS) %>%
  left_join(classDF, by = "id") -> vScoreTrainDF 

readRDS(vsvRDS) %>%
  left_join(classDF, by = "id") -> vScoreValidDF 

readRDS(cstRDS) %>%
  left_join(classDF, by = "id") -> cScoreTrainDF 

readRDS(csvRDS) %>%
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


# ------------------------------------------------
# Logistic Regression using 1 to 50 selected probes
#
probeUniDF %>%
  arrange(loss) %>%
  multiple_logistic(trainDF, 
                    validDF, 
                    "class") -> probeSelectedDF

  
# ------------------------------------------------
# PC names
#
pcaNames <- glue("PC{1:375}")

# ------------------------------------------------
# Covariance PCA (scale=FALSE)
# Loss from univariate logistic regressions
#
vpcaUniDF <- univariate_logistic(vScoreTrainDF, 
                                 "class", 
                                 pcaNames)


# ------------------------------------------------
# Multivariate Logistic using 1 to 50 selected PCs
# by univariate loss
#
vpcaUniDF %>%
  arrange(loss) %>%
  multiple_logistic(vScoreTrainDF, 
                    vScoreValidDF, 
                    "class") -> vpcaSelectedDF

# ------------------------------------------------
# Multivariate Logistic using first 1 to 50 PCs
#
vpcaUniDF %>%
   multiple_logistic(vScoreTrainDF, 
                     vScoreValidDF, 
                     "class") -> vpcaOrderedDF

# ------------------------------------------------
# Correlation PCA (scale=TRUE)
# Loss from univariate logistic regressions
#
cpcaUniDF <- univariate_logistic(cScoreTrainDF, 
                                 "class", 
                                 pcaNames)

# ------------------------------------------------
# Multivariate Logistic using 1 to 50 selected PCs
#
cpcaUniDF %>%
  arrange(loss) %>%
  multiple_logistic(cScoreTrainDF, 
                    cScoreValidDF, 
                    "class") -> cpcaSelectedDF

# ------------------------------------------------
# Multivariate Logistic using first 1 to 50 PCs
#
cpcaUniDF %>%
  multiple_logistic(cScoreTrainDF, 
                    cScoreValidDF, 
                    "class") -> cpcaOrderedDF

# ------------------------------------------------
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
