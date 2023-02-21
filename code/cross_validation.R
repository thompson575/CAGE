# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Feature selection comparison
# Selected Probes vs First few PCs vs Selected PCs
#
# Date: 19 Feb 2023
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
  
readRDS( file.path(home, "data/rData/training.rds")) %>%
  left_join( subjDF %>%
               select(id, diagnosis, class), by = "id") -> trainDF 


# =======================================================
# Selected probes
#
# Probe names
#
names(trainDF)[2:1001] -> pbNames

# =====================================================
# 10-fold CV of probe selection
#
cv_loss(trainDF, "class", pbNames, K = 10, maxSel = 10,
                  pca = FALSE, select = TRUE) %T>%
  saveRDS(file.path(home, "data/dataStore/cv_10_probe.rds")) -> cv10Probe
# ----------------------------------------------------
# Plot 10-fold CV of probe selection
#
cv_plot(cv10Probe, "Selected Probes: 10-fold CV")
# ----------------------------------------------------
#  list selected probes 
#
cv10Probe[[3]] %>%
  filter(order <= 5) %>%
  count(item) %>%
  arrange( desc(n))

# =====================================================
# 10-fold CV of PC selection
#
cv_loss(trainDF, "class", pbNames, K = 10, maxSel = 10,
                   pca = TRUE, select = TRUE) %T>%
  saveRDS(file.path(home, "data/dataStore/cv_10_pcSel.rds")) -> cv10PcSel
# ----------------------------------------------------
# Plot 10-fold CV of PC selection
#
cv_plot(cv10PcSel, "Selected PCs: 10-fold CV") 
# ----------------------------------------------------
#  list selected PCs 
#
cv10PcSel[[3]] %>%
  filter(order <= 5) %>%
  count(item) %>%
  arrange( desc(n))

# =====================================================
# 10-fold CV of leading PCs
#
cv_loss(trainDF, "class", pbNames, K = 10, maxSel = 10,
                   pca = TRUE, select = FALSE) %T>%
  saveRDS(file.path(home, "data/dataStore/cv_10_pcOrd.rds")) -> cv10PcOrd
# ----------------------------------------------------
# Plot 10-fold CV of Initial PCs
#
cv_plot(cv10PcOrd, "Initial PCs: 10-fold CV")

# =====================================================
# Leave one out CV of probe selection
#
cv_loss(trainDF, "class", pbNames, K = 375, maxSel = 10,
        pca = FALSE, select = TRUE) %T>%
  saveRDS(file.path(home, "data/dataStore/cv_loo_probe.rds")) -> cvLooProbe
# ----------------------------------------------------
# Plot Leave one out CV of probe selection
#
cv_plot(cvLooProbe, "Selected Probes: Leave one out CV")
# ----------------------------------------------------
#  list selected probes 
#
cvLooProbe[[3]] %>%
  filter(order <= 5) %>%
  count(item) %>%
  arrange( desc(n))

# =====================================================
# Leave one out CV of PC selection
#
cv_loss(trainDF, "class", pbNames, K = 375, maxSel = 10,
        pca = TRUE, select = TRUE) %T>%
  saveRDS(file.path(home, "data/dataStore/cv_loo_pcSel.rds")) -> cvLooPcSel
# ----------------------------------------------------
# Plot Leave one out CV of PC selection
#
cv_plot(cvLooPcSel, "Selected PCs: Leave one out CV") 
# ----------------------------------------------------
#  list selected PCs 
#
cvLooPcSel[[3]] %>%
  filter(order <= 5) %>%
  count(item) %>%
  arrange( desc(n))

# =====================================================
# Leave one out CV of leading PCs
#
cv_loss(trainDF, "class", pbNames, K = 375, maxSel = 10,
                   pca = TRUE, select = FALSE) %T>%
  saveRDS(file.path(home, "data/dataStore/cv_loo_pcOrd.rds")) -> cvLooPcOrd
# --- plot the results --------------------------------
cv_plot(cvLooPcOrd, "Initial PCs: Leave one out CV")
