#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

# Read the training data

home <- "C:/Projects/RCourse/Masterclass/CAGE"

patientDF  <- readRDS( file.path(home, "data/cache/patients.rds") )

trainDF <- readRDS( file.path(home, "data/cache/training.rds"))

# names of the 1000 probes
probes <- names(trainDF)[-1]

trainDF %>%
  left_join(patientDF, by = "id") -> trainDF

# proportion of Cancer cases in AEGIS-1
pC <- sum(trainDF$diagnosis == "Cancer") / nrow(trainDF)


# Define UI: histograms & plots proportions
ui <- fluidPage(

    # Application title
    titlePanel("CAGE: 1000 Probes from AEGIS-1"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("select", h3("Probe"), 
                        choices = probes, 
                        selected = probes[1])
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           plotOutput("propPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
      tibble( diagnosis = trainDF$diagnosis,
              probe     = trainDF[[input$select]] ) %>%
        ggplot( aes(x=probe)) +
        geom_histogram( bins=50 , fill = "steelblue") +
        facet_grid( diagnosis ~ .) +
        labs( x = "Expression")
    })
    
    output$propPlot <- renderPlot({
       tibble( diagnosis = trainDF$diagnosis,
               x         = trainDF[[input$select]],
               cat       = cut(x,
                               breaks = quantile(x,
                                                 probs = seq(0, 1, 0.1)), 
                                include.lowest = TRUE, labels = 1:10) ) %>%
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

# Run the application 
shinyApp(ui = ui, server = server)
