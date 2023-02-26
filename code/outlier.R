# -----------------------------------------------------------
# CAGE: Classification After Gene Expression
#
# Outlier Detection
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
# Univariate Search for Outliers
#
pbNames <- names(trainDF)[-1]

outProbe <- NULL
outExpr  <- NULL
outFactor <- NULL
for( probe in pbNames ) {
  y <- sort(trainDF[[probe]])
  q <- quantile(y, probs=c(0.25, 0.75))
  ub <- q[2] + 10 * (q[2] - q[1])
  lb <- q[1] - 10 * (q[2] - q[1])
  out <- y[ y < lb | y > ub]
  for( i in seq_along(out)) {
    outProbe <- c(outProbe, probe)
    outExpr  <- c(outExpr, out[i])
    f <- ifelse( out[i] > ub, (out[i] - q[2]) / (q[2] - q[1]),
                 (out[i] - q[1]) / (q[2] - q[1]))
    outFactor <- c(outFactor, f)
  }
}

tibble( probe = outProbe,
        expression = outExpr,
        factor = outFactor) %>%
  print() -> outlierDF

trainDF %>%
  pivot_longer(-id, names_to="probe", values_to="expression") %>%
  filter( probe %in% unique(outlierDF$probe)) %>%
  ggplot( aes(x=probe, y=expression)) +
  geom_boxplot()

trainDF %>%
  ggplot( aes(x=ENSG00000016490_at , y=ENSG00000185966_at)) +
  geom_point()

for( probe in pbNames ) {
  y <- trainDF[[probe]]
  q <- quantile(y, probs=c(0.25, 0.75))
  ub <- q[2] + 3 * (q[2] - q[1])
  lb <- q[1] - 3 * (q[2] - q[1])
  trainDF[[probe]] <- pmax( pmin(ub, y), lb )
}