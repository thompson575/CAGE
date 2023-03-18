# --------------------------------------------------
# Pipeline for the CAGE project
#
library(makepipe)

# --- Home directory
home <- "C:/Projects/RCourse/Masterclass/CAGE"
setwd(home)

# --------------------------------------------------
# Read Data
#
make_with_source (
  source       = "code/read_data.R",
  targets      = c("data/cache/patients.rds",
                   "data/cache/expression.rds",
                   "data/cache/training.rds",
                   "data/cache/validation.rds"),
  dependencies = NULL
)

# --------------------------------------------------
# Check the data
#
make_with_recipe(
  recipe       = rmarkdown::render( "reports/data_check.rmd"),
  targets      = "reports/data_check.html",
  dependencies = c("data/cache/patients.rds",
                   "data/cache/expression.rds")
)

# --------------------------------------------------
# Read Data Report
#
make_with_recipe(
  recipe       = rmarkdown::render( "reports/data_report.rmd"),
  targets      = "reports/data_report.html",
  dependencies = c("data/cache/patients.rds",
                   "data/cache/expression.rds")
)

# --------------------------------------------------
# PCA of 1000 probes
#
make_with_source (
  source       = "code/pca.R",
  targets      = c("data/cache/train_vpca.rds",
                   "data/cache/train_vpca_scores.rds",
                   "data/cache/valid_vpca_scores.rds",
                   "data/cache/train_cpca.rds",
                   "data/cache/train_cpca_scores.rds",
                   "data/cache/valid_cpca_scores.rds"),
  dependencies = c("data/cache/patients.rds",
                   "data/cache/training.rds",
                   "data/cache/validation.rds")
)

# --------------------------------------------------
# Informal check of PCA of 1000 probes
#
make_with_recipe(
  recipe       = rmarkdown::render( "reports/pca_check.rmd"),
  targets      = "reports/pca_check.html",
  dependencies = c("data/cache/patients.rds",
                   "data/cache/training.rds",
                   "data/cache/validation.rds",
                   "data/cache/train_vpca.rds",
                   "data/cache/train_vpca_scores.rds",
                   "data/cache/valid_vpca_scores.rds",
                   "data/cache/train_cpca.rds",
                   "data/cache/train_cpca_scores.rds",
                   "data/cache/valid_cpca_scores.rds",
                   "code/display_functions.R")
)
# --------------------------------------------------
# Formal report on PCA of 1000 probes
#
make_with_recipe(
  recipe       = rmarkdown::render( "reports/pca_report.rmd"),
  targets      = "reports/pca_report.html",
  dependencies = c("data/cache/patients.rds",
                   "data/cache/training.rds",
                   "data/cache/validation.rds",
                   "data/cache/train_vpca.rds",
                   "data/cache/train_vpca_scores.rds",
                   "data/cache/valid_vpca_scores.rds",
                   "data/cache/train_cpca.rds",
                   "data/cache/train_cpca_scores.rds",
                   "data/cache/valid_cpca_scores.rds",
                   "code/display_functions.R")
)

# --------------------------------------------------
# Feature Selection of 1000 probes
#
make_with_source (
  source       = "code/features.R",
  targets      = "data/cache/feature_loss.rds",
  dependencies = c("data/cache/patients.rds",
                   "data/cache/training.rds",
                   "data/cache/validation.rds",
                   "data/cache/train_vpca_scores.rds",
                   "data/cache/valid_vpca_scores.rds",
                   "data/cache/train_cpca_scores.rds",
                   "data/cache/valid_cpca_scores.rds",
                   "code/calculation_functions.R"))

# --------------------------------------------------
# Informal check of feature selection
#
make_with_recipe(
  recipe       = rmarkdown::render( "reports/features_check.rmd"),
  targets      = "reports/features_check.html",
  dependencies = c("data/cache/feature_loss.rds",
                   "data/cache/patients.rds",
                   "data/cache/training.rds",
                   "data/cache/validation.rds",
                   "data/cache/train_vpca_scores.rds",
                   "data/cache/valid_vpca_scores.rds",
                   "data/cache/train_cpca_scores.rds",
                   "data/cache/valid_cpca_scores.rds",
                   "code/display_functions.R")
)
# --------------------------------------------------
# Formal report on feature selection
#
make_with_recipe(
  recipe       = rmarkdown::render( "reports/features_report.rmd"),
  targets      = "reports/features_report.html",
  dependencies = c("data/cache/feature_loss.rds",
                   "data/cache/patients.rds",
                   "data/cache/training.rds",
                   "data/cache/validation.rds",
                   "data/cache/train_vpca_scores.rds",
                   "data/cache/valid_vpca_scores.rds",
                   "data/cache/train_cpca_scores.rds",
                   "data/cache/valid_cpca_scores.rds",
                   "code/display_functions.R")
)
# --------------------------------------------------
# Feature Selection for ALL probes
#
make_with_source (
  source       = "code/aegis1.R",
  targets      = c("data/cache/expression_aegis1.rds",
                   "data/cache/expression_aegis2.rds",
                   "data/cache/aegis1_vpca.rds",
                   "data/cache/aegis1_vpca_scores.rds",
                   "data/cache/aegis2_vpca_scores.rds",
                   "data/cache/aegis1_cpca.rds",
                   "data/cache/aegis1_cpca_scores.rds",
                   "data/cache/aegis2_cpca_scores.rds",
                   "data/cache/aegis1_loss.rds"),
  dependencies = c("data/cache/patients.rds",
                   "data/cache/expression.rds",
                   "code/calculation_functions.R")
)

make_with_recipe(
  recipe       = rmarkdown::render( "reports/aegis_presentation.rmd"),
  targets      = "reports/aegis_presentation.html",
  dependencies = c("data/cache/patients.rds",
                   "data/cache/aegis1_loss.rds",
                   "code/display_functions.R")
)

make_with_recipe(
  recipe       = rmarkdown::render( "reports/plos_article.rmd"),
  targets      = "reports/plos_article.pdf",
  dependencies = c("data/cache/patients.rds",
                   "data/cache/aegis1_loss.rds",
                   "code/display_functions.R")
)

get_pipeline()$segments

show_pipeline()
show_pipeline(as = "visnetwork")
