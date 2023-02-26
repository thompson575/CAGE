# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Analysis of full data for AEGIS-1
#
# Date: 19 Feb 2023
#
library(tidyverse)
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
  saveRDS( file.path(home, "data/rData/expression_aegis1.rds") ) -> aegis1DF

readRDS( file.path(home, "data/rData/expression.rds") ) %>%
  { .[, c(1, 377:506)] } %>%
  pivot_longer(-ID_REF, names_to = "id", values_to = "expression") %>%
  pivot_wider(names_from = ID_REF, values_from = expression) %>%
  print() %T>%
  saveRDS( file.path(home, "data/rData/expression_aegis1.rds") ) -> aegis2DF

# --------------------------------------------
# pca of the aegis1 data
#
aegis1DF %>%
  select( -id) %>%
  as.matrix() %>% 
  prcomp() %T>%
  saveRDS( file.path(home, "data/rData/aegis1_pca.rds")) -> pca

# --------------------------------------------
# pca scores of the training data
#
predict(pca) %>%
  as_tibble() %>%
  mutate( id = aegis1DF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, "data/rData/aegis1_pca_scores.rds")) -> scoreAegis1DF

# --------------------------------------------
# pca scores of the validation data
#
predict(pca, newdata=aegis2DF) %>%
  as_tibble() %>%
  mutate( id = aegis2DF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, "data/rData/aegis2_pca_scores.rds")) -> scoreAegis2DF

# --------------------------------------------
# Scaled pca of the aegis1 data
#
aegis1DF %>%
  select( -id) %>%
  as.matrix() %>% 
  prcomp( scale=TRUE) %T>%
  saveRDS( file.path(home, "data/rData/aegis1_scaled_pca.rds")) -> pcaScaled

# --------------------------------------------
# pca scores of the training data
#
predict(pcaScaled) %>%
  as_tibble() %>%
  mutate( id = aegis1DF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, "data/rData/aegis1_scaled_pca_scores.rds")) -> scoreAegis1ScaledDF

# --------------------------------------------
# pca scores of the validation data
#
predict(pcaScaled, newdata=aegis2DF) %>%
  as_tibble() %>%
  mutate( id = aegis2DF$id ) %>%
  relocate(id) %>%
  print() %T>%   
  saveRDS(file.path(home, "data/rData/aegis2_scaled_pca_scores.rds")) -> scoreAegis2ScaledDF

# =======================================================
# Selected probes
#
# Probe names
#
names(aegis1DF)[-1] -> pbNames

aegis1DF %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> aegis1DF


aegis2DF %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> aegis2DF

# ------------------------------------------------
# Loss from univariate logistic regressions
#
pbUniDF <- uni_logistic(aegis1DF, "class", pbNames) 

# ----------------------------------------------------------
# Plot of the 5 best probes in training and validation data
#
pbUniDF %>% arrange(loss) %>% print() %>% slice(1:5) %>% pull(x) -> topX

aegis1DF %>%
  pivot_longer( starts_with("ENSG"), names_to="probe", values_to="expression") %>%
  filter( probe %in% topX ) %>%
  ggplot( aes(x=probe, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "AEGIS-1: 5 most predictive probes",
        x = "") +
  theme( legend.position = c(0.85, 0.85))

aegis2DF %>%
  pivot_longer( starts_with("ENSG"), names_to="probe", values_to="expression") %>%
  filter( probe %in% topX ) %>%
  ggplot( aes(x=probe, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "AEGIS-2: 5 most predictive probes",
        x = "") +
  theme( legend.position = c(0.85, 0.85))

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 selected predictors
#
pbSelDF <- sel_logistic(aegis1DF, aegis2DF, "class", pbUniDF, maxSel = 50)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pbSelDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "AEGIS-1: Loss from Selected Probes",
        y     = "Loss",
        x     = "Number of Probes") +
  geom_text( x=35, y=0.58, label="In-sample", colour="blue") +
  geom_text( x=35, y=1.20, label="Out-of-sample", colour="red") 

# =======================================================
# Selected Unscaled Principal Components
#
# PC names
#
pcNames <- paste0("PC", 1:375)

scoreAegis1DF %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> scoreAegis1DF


scoreAegis2DF %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> scoreAegis2DF

# ------------------------------------------------
# Loss from univariate logistic regressions
#
pcUniDF <- uni_logistic(scoreAegis1DF, "class", pcNames)

# ----------------------------------------------------------
# Plot of the 5 best PCs in training and validation data
#
pcUniDF %>% arrange(loss) %>% print() %>% slice(1:5) %>% pull(x) -> topX

scoreAegis1DF %>%
  pivot_longer( starts_with("PC"), names_to="pc", values_to="expression") %>%
  filter( pc %in% topX ) %>%
  ggplot( aes(x=pc, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Aegis-1: 5 most predictive PCs",
        x = "") +
  theme( legend.position = c(0.15, 0.85))

scoreAegis2DF %>%
  pivot_longer( starts_with("PC"), names_to="pc", values_to="expression") %>%
  filter( pc %in% topX ) %>%
  ggplot( aes(x=pc, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Aegis-2: 5 most predictive PCs",
        x = "") +
  theme( legend.position = c(0.15, 0.85))

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 selected PCs
#
pcSelDF <- sel_logistic(scoreAegis1DF, scoreAegis2DF, "class", pcUniDF, maxSel = 50)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pcSelDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "Aegis-1: Loss from Selected PCs",
        y     = "Loss",
        x     = "Number of PCs") +
  geom_text( x=35, y=0.45, label="In-sample", colour="blue") +
  geom_text( x=35, y=0.8, label="Out-of-sample", colour="red")

# =======================================================
# Initial Unscaled Principal Components
#

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 ordered PCs
#
pcOrdDF <- sel_logistic(scoreAegis1DF, scoreAegis2DF, "class", pcUniDF, 
                        maxSel = 50, sort = FALSE)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pcOrdDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "Aegis-1: Loss from Ordered PCs",
        y     = "Loss",
        x     = "Number of PCs") +
  geom_text( x=35, y=0.60, label="In-sample", colour="blue") +
  geom_text( x=35, y=0.95, label="Out-of-sample", colour="red")

# =======================================================
# Selected Unscaled Principal Components
#
# PC names
#
pcNames <- paste0("PC", 1:375)

scoreAegis1ScaledDF %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> scoreAegis1ScaledDF


scoreAegis2ScaledDF %>%
  left_join( subjDF %>%
               mutate( class = ifelse(diagnosis == "Cancer", 1, 0)) %>%
               select(id, diagnosis, class), by = "id") -> scoreAegis2ScaledDF

# ------------------------------------------------
# Loss from univariate logistic regressions
#
pcScdUniDF <- uni_logistic(scoreAegis1ScaledDF, "class", pcNames)

# ----------------------------------------------------------
# Plot of the 5 best PCs in training and validation data
#
pcScdUniDF %>% arrange(loss) %>% print() %>% slice(1:5) %>% pull(x) -> topX

scoreAegis1ScaledDF %>%
  pivot_longer( starts_with("PC"), names_to="pc", values_to="expression") %>%
  filter( pc %in% topX ) %>%
  ggplot( aes(x=pc, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Aegis-1: 5 most predictive scaled PCs",
        x = "") +
  theme( legend.position = c(0.85, 0.25))

scoreAegis2ScaledDF %>%
  pivot_longer( starts_with("PC"), names_to="pc", values_to="expression") %>%
  filter( pc %in% topX ) %>%
  ggplot( aes(x=pc, y=expression, fill=diagnosis)) +
  geom_boxplot() +
  coord_flip() +
  labs( title = "Aegis-2: 5 most predictive scaled PCs",
        x = "") +
  theme( legend.position = c(0.85, 0.25))

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 selected PCs
#
pcScdSelDF <- sel_logistic(scoreAegis1ScaledDF, scoreAegis2ScaledDF, 
                           "class", pcScdUniDF, maxSel = 50)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pcScdSelDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "Aegis-1: Loss from Selected scaled PCs",
        y     = "Loss",
        x     = "Number of PCs") +
  geom_text( x=35, y=0.45, label="In-sample", colour="blue") +
  geom_text( x=35, y=0.9, label="Out-of-sample", colour="red")

# =======================================================
# Initial Unscaled Principal Components
#

# ------------------------------------------------------------------
# Multivariable Logistic using 1 to 50 ordered PCs
#
pcScdOrdDF <- sel_logistic(scoreAegis1ScaledDF, scoreAegis2ScaledDF, 
                           "class", pcScdUniDF, maxSel = 50, sort = FALSE)

# ------------------------------------------------------------------
# Plot of in-sample vs out-of-sample loss
#
pcScdOrdDF %>%
  ggplot( aes(x = n, y = inloss)) +
  geom_line( size=1.2, colour="blue") +
  geom_line( aes(y = outloss), size=1.2, colour="red") +
  labs( title = "Aegis-1: Loss from Ordered scaled PCs",
        y     = "Loss",
        x     = "Number of PCs") +
  geom_text( x=35, y=0.60, label="In-sample", colour="blue") +
  geom_text( x=35, y=0.95, label="Out-of-sample", colour="red")

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
  labs( title = "AEGIS-1: In-sample loss for the 5 feature selection methods",
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
  labs( title = "AEGIS-1: In-sample loss for the 5 feature selection methods",
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
  labs( title = "AEGIS-1: Out-of-sample loss for the 5 feature selection methods",
        y     = "Loss",
        x     = "Number of predictors") +
  theme( legend.position = c(0.15, 0.75))

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
  labs( title = "AEGIS-1: Out-of-sample loss for the 5 feature selection methods",
        y     = "Loss",
        x     = "Number of predictors") +
  theme( legend.position = c(0.15, 0.75))

# ---------------------------------------------
# Save results to dataStore
#
saveRDS( list( pcUniDF = pcUniDF, pcSelDF = pcSelDF, 
               pcOrdDF = pcOrdDF, 
               pbUniDF = pbUniDF, pbSelDF = pbSelDF,
               pcScdUniDF = pcScdUniDF, pcScdSelDF = pcScdSelDF, 
               pcScdOrdDF = pcScdOrdDF),
         file.path(home, "data/dataStore/aegis1_classification_loss.rds"))
