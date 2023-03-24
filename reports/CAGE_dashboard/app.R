# -----------------------------------------------------------
# Project CAGE: 
# dashboard for exploration of the 1000 probe training data
#
library(shiny)
library(tidyverse)
library(fs)

# Read the training data
cache <- "C:/Projects/RCourse/Masterclass/CAGE/data/cache"
patientDF <- readRDS(path(cache, "patients.rds") )
trainDF   <- readRDS(path(cache, "training.rds"))

# names of the 1000 probes
probes <- names(trainDF)[-1]

# combine training and patient data
trainDF %>%
  left_join(patientDF, by = "id") -> trainDF

# proportion of Cancer cases in training data
pC <- sum(trainDF$diagnosis == "Cancer") / nrow(trainDF)

# --------------------------------------------------
# Define UI:
#
ui <- fluidPage(

    # Application title
    titlePanel("CAGE: 1000 Probes from AEGIS-1"),

    sidebarLayout(
      # Sidebar with probe selector 
      sidebarPanel(
        selectInput("select", 
                    h3("Probe"), 
                    choices  = probes, 
                    selected = probes[1])
      ),

      # Main panel containing the plots
      mainPanel(
           plotOutput("distPlot"),
           plotOutput("propPlot")
        )
    )
)

# --------------------------------------------------
# Define server 
#
server <- function(input, output) {
    # Histogram of expression faceted by diagnosis
    output$distPlot <- renderPlot({
      tibble( diagnosis = trainDF$diagnosis,
              probe     = trainDF[[input$select]] ) %>%
        ggplot( aes(x=probe)) +
        geom_histogram( bins=50 , fill = "steelblue") +
        facet_grid( diagnosis ~ .) +
        labs( x = "Expression")
    })
    # Plot of the proportion of cancer cases
    output$propPlot <- renderPlot({
       tibble( diagnosis = trainDF$diagnosis,
               x         = trainDF[[input$select]],
               cat       = cut(x,
                               breaks = quantile(x,
                                                 probs = seq(0, 1, 0.1)), 
                               include.lowest = TRUE, 
                               labels = 1:10) ) %>%
        group_by( cat ) %>%
        summarise( n = n(),
                c = sum( diagnosis == "Cancer"),
                x = mean(x)) %>%
        mutate( p = c / n,
                se = sqrt( p * (1-p) / n)) %>%
        ggplot( aes(x=x, y=p)) +
        geom_point(size=2, colour="steelblue") +
        geom_errorbar( aes(ymin=p-se, ymax=p+se), width=0.025,
                       colour="steelblue") +
        geom_hline( yintercept = pC, lty="dashed", colour="red") +
        labs( x = "Mean expression",
              y = "Proportion Cancer")
    })
}

# --------------------------------------------------
# Run the application 
#
shinyApp(ui = ui, server = server)
