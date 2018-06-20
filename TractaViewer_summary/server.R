
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(plotly)
library(plyr)
library(dplyr)
library(readxl)

options(shiny.maxRequestSize=30*1024^2)

testfile = "../Example_output/sample_output_druggability.xlsx"
# fl = "TEST STRING DO NOT USE OR MODIFY"
fl = "INTENTIONALLY A STRING"
freqs = NULL
buckets = NULL

xPts <- runif(1000, 0, 1)
yPts <- runif(1000, 0, 1)
zPts <- runif(1000, 0, 1)

shinyServer(function(input, output, session) {

  dim1_names <- reactive({
    dim1_names <- colnames(buckets)[4:length(colnames(buckets))]
  })

  dim2_names <- reactive({
    dim2_names <- colnames(buckets)[4:length(colnames(buckets))]
  })

  dim3_names <- reactive({
    dim3_names <- colnames(buckets)[4:length(colnames(buckets))]
  })
  
  observe({
    # I'll try to do something cleverer with this - thereshould be no need to re-read the file every time a 
    # change is requested.
    # It can be done, and this code does it. Note that we have to assign to varibles outside the scope of this
    # function using the <<- operator. If we don't do that, values are lost, even if assigned to variables
    # declared in the widest scope at the top!
    
    if (is.null(input$file1)) {
      infile = testfile
    } else {
      infile = input$file1[['datapath']]
    }
    
    if (fl != infile) {
      fl <<- infile
      print(fl)
      buckets <<- read_excel(fl, sheet = "Buckets")
      colnames(buckets) <<- colnames(buckets) %>% stringr::str_replace_all("Buckets.","") # A bit of cleanup needed
      colnames(buckets) <<- colnames(buckets) %>% stringr::str_replace_all(".bucket","")
      # We want to count the number of genes falling into each dimension - because that's going to be our point
      # size. (Plotly can probably work this out if asked, but doing it myself makes a bit more sense to me).
      updateSelectInput(session = session, inputId = "dim1", choices = dim1_names(), selected = colnames(buckets)[4])
      updateSelectInput(session = session, inputId = "dim2", choices = dim2_names(), selected = colnames(buckets)[5])
      updateSelectInput(session = session, inputId = "dim3", choices = dim3_names(), selected = colnames(buckets)[6])
    }
    
  })
  
  output$plot <- renderPlotly({
    # Need to re-do freqs calculation here.
    # If we don't, we'll possibly get multiple points appearing at some grid positions - because records are being split by
    # buckets we aren't currently seeing.
    cnames = make.names(c(input$dim1, input$dim2, input$dim3, "n"))
    print(cnames)
    freqs <- buckets %>% group_by(get(input$dim1), get(input$dim2), get(input$dim3)) %>% tally()
    colnames(freqs) = cnames
    freqs$Overall.score <- apply(freqs[cnames[1:3]], 1, mean)

    # This turns out to be by far the easiest way of dealing with a dynamic range of possible point sizes.
    szs = c(15, 15)
    if (max(freqs$n) < 3) {
      szs = c(10, 15)
    } else if (max(freqs$n) < 10) {
      szs = c(10, 30)
    } else {
      szs = c(5, 50)
    }

    m <- list(
      l = 50,
      r = 50,
      b = 100,
      t = 100,
      pad = 4
    )
    
    print(as.data.frame(freqs))
    
    plot_ly(as.data.frame(freqs),  # This won't play nice with tibbles, it has to have a data frame.
            type = "scatter3d",
            x = as.formula(paste0("~", cnames[1])),
            y = as.formula(paste0("~", cnames[2])),
            z = as.formula(paste0("~", cnames[3])),
            size = ~n,
            color = ~Overall.score,
            marker = list(symbol = 'circle', sizemode = 'diameter'),
            colors=c("green", "red"),
            sizes = szs,
            hoverinfo = 'text',
            text = ~paste(input$dim1, ':', get(cnames[1]),
                          '<br>', input$dim2, ':', get(cnames[2]),
                          '<br>', input$dim3, ':', get(cnames[3]),
                          '<br>Num. targets:', n, sep = "")
    ) %>% add_markers(showlegend=FALSE) %>% colorbar(title = "Overall score")
  })

  output$targetsHoverTable <- renderTable({
    # Read in hover data
    eventdata <- event_data("plotly_click")
    validate(need(!is.null(eventdata), "Click on points in the plot to populate this table"))
    targets = buckets %>% filter(get(input$dim1) == eventdata[1, "x"] &
                                   get(input$dim2) == eventdata[1, "y"] &
                                   get(input$dim3) == eventdata[1, "z"])
    targets[c("series", "HGNC Name", "GeneID")]
  })

  output$click <- renderPrint({
    d <- event_data("plotly_click")
    if (is.null(d)) s = "None" else s = paste(paste0(input$dim1, ": ", d[1, "x"]),
                                              paste0(input$dim2, ": ", d[1, "y"]),
                                              paste0(input$dim3, ": ", d[1, "z"]), sep = ", ")
    paste("Axis values:", s, sep = " ")
  })
  
  
})
