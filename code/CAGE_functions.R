library(tidyverse)
library(glue)

# =====================================================
# Calculate the mean cross-entropy loss
# 
calc_loss <- function(thisModel,   # the model
                      thisDF       # the data
                      ) {
  response <- names(thisModel$model)[1]
  y        <- thisDF[[ response]]
  yhat     <- predict(thisModel, newdata=thisDF, type="response")
  return( - mean( y * log(yhat) + (1 - y) * log(1 - yhat) ) )
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
    lossDF$loss[i] <- calc_loss(model, thisDF)
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
    
    selLossDF$inloss[M]  <- calc_loss(model, thisDF)
    selLossDF$outloss[M] <- calc_loss(model, validDF)
  }
  # -- return the results -------------------
  return( selLossDF )
}

# =====================================================
# K-fold cross-validation of loss for gene or PCs
#
cv_loss <- function(thisDF,      # the original data
                    response,    # response variable
                    predictors,  # the potential predictors
                    K=10,        # number of folds
                    maxSel=10,   # largest number to select
                    select=TRUE, # whether to select
                    pca=FALSE    # whether to take PCs
) {
  # -- number of subjects -------------------------
  N <- nrow(thisDF)
  # -- check that K <= N --------------------------
  if( K > N ) {
    stop("Cannot run CV with K > N")
  }
  # -- matrices to store results ------------------
  testLoss   <- matrix(0, nrow=maxSel, ncol=K)
  trainLoss  <- matrix(0, nrow=maxSel, ncol=K)
  # -- variable to contain the selected items -----
  selection <- NULL
  # -- fold indicator -----------------------------
  thisDF$fold <- rep(1:K, times=ceiling(N/K))[1:N]
  
  # -- loop over the folds ------------------------
  for( k in 1:K ) {
    # -- extract training and test data -----------
    thisDF %>%
      filter( fold != k ) %>%
      select( -fold ) -> trainDF
    thisDF %>%
      filter( fold == k ) %>%
      select( -fold ) -> testDF
    
    # -- PCs if required --------------------------
    if( pca ) {
      trainDF %>%
        select( all_of(predictors)) %>%
        as.matrix() %>% 
        prcomp() -> pcs    
      
      trainDF %>% 
        select(-all_of(predictors)) %>%
        bind_cols(predict(pcs) %>%
                    as_tibble()) -> trainDF
      testDF %>% 
        select(-all_of(predictors)) %>%
        bind_cols(predict(pcs, newdata=testDF) %>%
                    as_tibble()) -> testDF
      # -- names of the PCs ----------------------
      nPC <- length(pcs$sdev[pcs$sdev > 0])
      x <- glue("PC{1:nPC}" )
    } else {
      # -- names of the genes --------------------
      x <- predictors
    }
    # -- univariate Loss all candidates ------------
    uniDF <- uni_logistic(trainDF, response, x)
    # -- multivariate Loss for selected 
    selDF <- sel_logistic(trainDF, testDF, response, uniDF, 
                      maxSel=maxSel, sort=select)
    trainLoss[, k] <- selDF$inloss
    testLoss[, k]  <- selDF$outloss
    # # -- select by Loss if required ----------------
    if( select ) {
      uniDF %>%
        arrange( loss) %>%
        pull( x ) -> topXs
    } else {
      topXs <- uniDF$x
    }
    # -- add topXs to the record -----------------
    selection <- c(selection, topXs[1:maxSel])
  }
  # -- selected items into a DF ------------------
  selDF <- tibble( k = rep(1:K, each=maxSel),
                   order = rep(1:maxSel, times=K),
                   item = selection)
  # -- return the results ------------------------
  return( list(trainLoss = trainLoss, testLoss = testLoss,
               selDF = selDF) )
}

# =====================================================
# Plot Eigenvalues of a PCA
#
plot_eigenvalues <- function(pcStd, pctVar = FALSE) {
  if( pctVar ) {
    y      <- 100 * pcStd * pcStd / sum(pcStd * pcStd)
    yLabel <- "Percent Variance"
  } else {
    y      <- pcStd
    yLabel <- "Standard Deviation"
  }
  tibble( x  = 1:length(y),
          y  = y) %>%
    ggplot( aes(x = x, y = y, xend = x, yend = 0)) +
    geom_segment() +
    labs( y     = yLabel,
          x     = "Principal Component Number")
}
# =====================================================
# Boxplot of features
#
boxplot_features <- 
  function( predictors,  # names of predictors to be plotted
            thisDF       # data for plotting 
  ) {
    thisDF %>%
      pivot_longer( all_of(predictors), names_to="var", values_to="y") %>%
      ggplot( aes(x=var, y=y, fill=diagnosis)) +
      geom_boxplot() +
      coord_flip() +
      labs( x = "", fill="") +
      theme( legend.position = "bottom")
  }

# =====================================================
# Plot in-sample and out-of-sample loss
#
plot_in_out_loss <- 
  function(selDF  # multivariate losses 
  ) {
    i        <- floor(0.7 * nrow(selDF))
    xText    <- selDF$M[i]
    yInText  <- 1.05 * selDF$inloss[i]
    yOutText <- 0.95 * selDF$outloss[i]
    selDF %>%
      ggplot( aes(x = M, y = inloss)) +
      geom_line( size=1.2, colour="blue") +
      geom_line( aes(y = outloss), size=1.2, colour="red") +
      geom_text( x=xText, y=yInText, label="In-sample", colour="blue") +
      geom_text( x=xText, y=yOutText, label="Out-of-sample", colour="red") 
  }

# =====================================================
# Plot loss comparing 5 methods
#
plot_method_comparison <- 
  function( methodResults, loss 
  ) {
    method_labels <- names(methodResults)
    methodResults[[1]] %>%
      rename( m1 = {{loss}}) %>%
      select( M, m1) %>%
      left_join( methodResults[[2]] %>% 
                   rename( m2 = {{loss}}) %>%
                   select( M, m2), by = "M") %>%
      left_join( methodResults[[3]] %>% 
                   rename( m3 = {{loss}}) %>%
                   select( M, m3), by = "M") %>%
      left_join( methodResults[[4]] %>% 
                   rename( m4 = {{loss}}) %>%
                   select( M, m4), by = "M") %>%
      left_join( methodResults[[5]] %>% 
                   rename( m5 = {{loss}}) %>%
                   select(M, m5), by = "M") %>%
      pivot_longer(m1:m5, names_to = "method", values_to = "loss") %>%
      mutate( method = factor(method, levels=c("m1", "m2", "m3", "m4", "m5"),
                              labels = method_labels)) %>%
      ggplot( aes(x = M, y = loss, colour=method)) +
      geom_line( size=1.2) +
      labs( y      = "Loss",
            x      = "Number of predictors",
            colour = "") +
      theme( legend.position = "bottom")
  }

# =====================================================
# Residual plots for a fitted model
#
resid_plots <- function(model) {
  
  # -- Residuals vs Fitted values --------------
  augment(model) %>%
    ggplot( aes(x=.fitted, y=.resid)) +
    geom_point() +
    geom_hline(yintercept=c(-1.96*glance(model)$sigma, 
                            0, 
                            1.96*glance(model)$sigma ),
               lty=c(2,1,2)) +
    labs(title="Residuals vs Fitted Values",
         y = "Residuals",
         x = "Fitted values") -> p1
  
  # -- QQ plot ---------------------------------
  augment(model) %>%
    arrange( .resid ) %>%
    mutate( p = (row_number() - 0.5)/length(.resid),
            z = qnorm(p)) %>%
    ggplot( aes(x=z, y=.resid)) +
    geom_point() +
    geom_smooth( method="lm") +
    labs(title="Normal QQ plot of the residuals",
         y = "Residuals",
         x = "Expected Normal Values")  -> p2
  # -- return list of plots ---------------------
  return( list(p1, p2) )
}

# =====================================================
# plot results from cv_loss
#
cv_plot <- function( cvLoss,  # result from cv_loss  
                     title  # title for the plot
) {
  # -- number of folds used --------------------- 
  K     <- ncol(cvLoss[[2]])
  xText <- floor( nrow(cvLoss[[2]]) * 0.75)
  # -- Loss values to plot --------------------- 
  tibble( n  = 1:nrow(cvLoss[[2]]),
          y  = apply(cvLoss[[2]], 1, mean),
          se = apply(cvLoss[[2]], 1, sd) / sqrt(K),
          iy = apply(cvLoss[[1]], 1, mean)) -> lossDF
  # -- make the plot -------------------------
  lossDF %>%
    ggplot( aes(x=n, y=y) ) + 
    geom_ribbon( aes(ymin=y-se, ymax=y+se), fill="lightblue") +
    geom_point(size=2, colour="blue") +
    geom_line(colour="blue") +
    geom_line( aes(y=iy), colour="red") +
    geom_point( aes(y=iy), colour="red", size=2) +
    geom_text( aes(x=xText, y=iy[xText] + se[xText]/2, 
                   label="In-sample"), colour="red",
               hjust = 0, nudge_x = 0.1) +
    geom_text( aes(x=xText, y=y[xText] + se[xText]/2, 
                   label="Cross-validation"),
               colour="blue", hjust = 0, nudge_x = 0.1) +
    labs( title = title,
          subtitle = "ribbon shows \u00B1 1 se",
          x = "Number of predictors",
          y = "Loss")
}