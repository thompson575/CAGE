library(tidyverse)
library(glue)

# =====================================================
# Calculate the mean cross-entropy loss
# 
cross_entropy <- function(thisModel,   # the model
                          thisDF       # the data
                          ) {
  response <- names(thisModel$model)[1]
  y        <- thisDF[[ response]]
  yhat     <- predict(thisModel, newdata=thisDF, type="response")
  return( - mean( y * log(yhat) + (1 - y) * log(1 - yhat) ) )
}

# =====================================================
# Calculate the mean sum of squares loss
# 
sum_squares <- function(thisModel,   # the model
                        thisDF       # the data
                        ) {
  response <- names(thisModel$model)[1]
  y        <- thisDF[[ response]]
  yhat     <- predict(thisModel, newdata=thisDF)
  return( mean( (y - yhat)^2 ) )
}


# =====================================================
# calculate in-sample losses for a set of predictors
# repeats univariate logistic regression  
#            response ~ single predictor
#
univariate_logistic <- 
  function(thisDF,     # the data
           response,   # name of response
           predictors  # names of potential predictors
) {
  # -- number of potential predictors --------------
  nPredictors <- length(predictors)
  # -- DF to contain the results -------------------
  tibble( id   = 1:nPredictors,
          x    = predictors,
          loss = rep(0, nPredictors)) -> lossDF
  # -- loop over potential predictors --------------
  for( i in 1:nPredictors ) {
    fml            <- as.formula( glue("{response} ~ {predictors[i]}"))
    model          <- glm(fml, data=thisDF, family="binomial")
    lossDF$loss[i] <- cross_entropy(model, thisDF)
  }
  # -- return DF of results ------------------------
  return( lossDF )
}

# =====================================================
# Loss (in-sample & out-of-sample) for
# Multivariate Logistic Regression
#        response ~ set of predictors
# predictors used are first M from lossDF  M=1..maxSel
#
multiple_logistic <- 
  function(lossDF,    # results from uni_logistic
           thisDF,    # the data
           validDF,   # validation data
           response,  # name of the response
           maxSel=50  # maximum number of predictors
) {
  # -- DF to contain the results ------------
  tibble( M       = 1:maxSel,
          inloss  = rep(0, maxSel),
          outloss = rep(0, maxSel)) -> selLossDF
  # -- loop over number of predictors -------
  for( M in 1:maxSel ) {
    thisDF %>% 
      select( all_of( c(lossDF$x[1:M], response)) ) %>%
      glm( as.formula( glue("{response} ~ .") ), data=.,
           family="binomial") -> model
    
    selLossDF$inloss[M]  <- cross_entropy(model, thisDF)
    selLossDF$outloss[M] <- cross_entropy(model, validDF)
  }
  # -- return the results -------------------
  return( selLossDF )
}

# =====================================================
# calculate in-sample losses for a set of predictors
# repeats univariate logistic regression  
#            response ~ single predictor
#
univariate_linear <- 
  function(thisDF,     # the data
           response,   # name of response
           predictors  # names of potential predictors
  ) {
    # -- number of potential predictors --------------
    nPredictors <- length(predictors)
    # -- DF to contain the results -------------------
    tibble( id   = 1:nPredictors,
            x    = predictors,
            loss = rep(0, nPredictors)) -> lossDF
    # -- loop over potential predictors --------------
    for( i in 1:nPredictors ) {
      fml            <- as.formula( glue("{response} ~ {predictors[i]}"))
      model          <- lm(fml, data=thisDF)
      lossDF$loss[i] <- sum_squares(model, thisDF)
    }
    # -- return DF of results ------------------------
    return( lossDF )
  }

# =====================================================
# Loss (in-sample & out-of-sample) for
# Multivariate Logistic Regression
#        response ~ set of predictors
# predictors used are first M from lossDF  M=1..maxSel
#
multiple_linear <- 
  function(lossDF,    # results from uni_logistic
           thisDF,    # the data
           validDF,   # validation data
           response,  # name of the response
           maxSel=50  # maximum number of predictors
  ) {
    # -- DF to contain the results ------------
    tibble( M       = 1:maxSel,
            inloss  = rep(0, maxSel),
            outloss = rep(0, maxSel)) -> selLossDF
    # -- loop over number of predictors -------
    for( M in 1:maxSel ) {
      thisDF %>% 
        select( all_of( c(lossDF$x[1:M], response)) ) %>%
        lm( as.formula( glue("{response} ~ .") ), data=.) -> model
      
      selLossDF$inloss[M]  <- sum_squares(model, thisDF)
      selLossDF$outloss[M] <- sum_squares(model, validDF)
    }
    # -- return the results -------------------
    return( selLossDF )
  }

