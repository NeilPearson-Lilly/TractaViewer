
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
# library(shinyRGL)
library(plotly)

dimensions = c("SM Druggability", "Safety", "Feasibility", "AB-ability", "Modality")


shinyUI(fluidPage(

  # Application title
  titlePanel("TractaViewer Data"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      # sliderInput("bins",
      #             "Number of bins:",
      #             min = 1,
      #             max = 50,
      #             value = 30),
      
      fileInput("file1", "Choose TractaViewer file",
                accept = c(
                  "Excel file",
                  ".xlsx")
      ),
      selectInput(
        "dim1", "Metric 1",
        dimensions, dimensions[1]
      ),
      selectInput(
        "dim2", "Metric 2",
        dimensions, dimensions[2]
      ),
      selectInput(
        "dim3", "Metric 3",
        dimensions, dimensions[3]
      ),
      submitButton("Apply Changes", icon("refresh"))
      # tags$hr(),
      # checkboxInput("header", "Header", TRUE)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      # plotOutput("distPlot"),
      # tableOutput("testTable")
      # webGLOutput("scatter3d")
      plotlyOutput("plot", height = 600),
      verbatimTextOutput("click"),
      tableOutput("targetsHoverTable")
    )
  )
))
