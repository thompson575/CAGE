library(tidyverse)

# =====================================================
# Plot Eigenvalues of a PCA
#
plot_eigenvalues <- 
  function(pcStd,          # vector of pc standard deviations
           pctVar = FALSE  # whether to plot percent variance
          ) {
  # --- plot percent variance or standard deviation
  if( pctVar ) {
    y      <- 100 * pcStd * pcStd / sum(pcStd * pcStd)
    yLabel <- "Percent Variance"
  } else {
    y      <- pcStd
    yLabel <- "Standard Deviation"
  }
  # --- make the requested plot
  tibble( x  = 1:length(y),
          y  = y) %>%
    ggplot( aes(x = x, y = y, xend = x, yend = 0)) +
    geom_segment() +
    labs( y     = yLabel,
          x     = "Principal Component Number")
}
# =====================================================
# Boxplot of selected features showing the two diagnoses
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
# Plot in-sample and out-of-sample loss for an
# increasing number of features
#
plot_in_out_loss <- 
  function(selDF  # multivariate in and out losses 
  ) {
    # --- position for the labels
    i        <- floor(0.7 * nrow(selDF))
    xText    <- selDF$M[i]
    yInText  <- 1.05 * selDF$inloss[i]
    yOutText <- 0.95 * selDF$outloss[i]
    # --- make the plot
    selDF %>%
      ggplot( aes(x = M, y = inloss)) +
      geom_line( linewidth=1.2, colour="blue") +
      geom_line( aes(y = outloss), linewidth=1.2, colour="red") +
      geom_text( x=xText, y=yInText,  label="In-sample", colour="blue") +
      geom_text( x=xText, y=yOutText, label="Out-of-sample", colour="red") 
  }

# =====================================================
# Plot loss comparing 5 methods for an
# increasing number of features
#
plot_method_comparison <- 
  function( methodResults, # losses for all methods
            loss           # inloss or outloss
  ) {
    # --- labels for the 5 methods
    method_labels <- names(methodResults)
    # --- extract loss for each method and join
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
      # --- make the plot
      ggplot( aes(x = M, y = loss, colour=method)) +
      geom_line( linewidth=1.2) +
      labs( y      = "Loss",
            x      = "Number of predictors",
            colour = "") +
      theme( legend.position = "bottom")
  }
