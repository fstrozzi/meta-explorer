
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)
library(SPARQL)
library(GetoptLong)

protein_to_pathway <- function(id){
  # Step 1 - Set up preliminaries and define query
  # Define the data.gov endpoint
  endpoint <- "http://www.ebi.ac.uk/rdf/services/reactome/sparql"
  
  # create query statement
  query <-qq("SELECT DISTINCT ?p ?o
             WHERE 
             {
             <http://www.reactome.org/biopax/49/170905#@{id}> ?p ?o
             }")

  # Step 2 - Use SPARQL package to submit query and save results to a data frame
  qd <- SPARQL(endpoint,query)
  df <- qd$results
  df
}

shinyServer(function(input, output) {

  output$table = DT::renderDataTable({
    print(input$protein)
    datatable(protein_to_pathway(input$protein))
  })
#   output$distPlot <- renderPlot({
# 
#     # generate bins based on input$bins from ui.R
#     x    <- faithful[, 2]
#     bins <- seq(min(x), max(x), length.out = input$bins + 1)
# 
#     # draw the histogram with the specified number of bins
#     hist(x, breaks = bins, col = 'darkgray', border = 'white')
# 
#   })

})
