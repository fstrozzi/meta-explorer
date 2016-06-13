
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)
library(SPARQL)

protein_to_pathway <- function(id){
  # Step 1 - Set up preliminaries and define query
  # Define the data.gov endpoint
  endpoint <- "http://www.ebi.ac.uk/rdf/services/reactome/sparql"
  
  # create query statement
  query <-
    "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  SELECT ?Protein ?label0 ?Provenance ?label1 ?BiochemicalReaction ?label2 ?Pathway 
  ?label3 
  WHERE { 
  ?Protein rdf:type <http://www.biopax.org/release/biopax-level3.owl#Protein>.
  OPTIONAL{
  ?Protein rdfs:label ?label0.}
  ?Protein <http://www.biopax.org/release/biopax-level3.owl#dataSource> ?Provenance.
  ?Provenance rdf:type <http://www.biopax.org/release/biopax-level3.owl#Provenance>.
  OPTIONAL{
  ?Provenance rdfs:label ?label1.}
  ?BiochemicalReaction <http://www.biopax.org/release/biopax-level3.owl#dataSource> ?Provenance.
  ?BiochemicalReaction rdf:type <http://www.biopax.org/release/biopax-level3.owl#BiochemicalReaction>.
  OPTIONAL{
  ?BiochemicalReaction rdfs:label ?label2.}
  ?Pathway <http://www.biopax.org/release/biopax-level3.owl#pathwayComponent> ?BiochemicalReaction.
  ?Pathway rdf:type <http://www.biopax.org/release/biopax-level3.owl#Pathway>.
  OPTIONAL{
  ?Pathway rdfs:label ?label3.}
  }LIMIT 100
  
  "
 
  # Step 2 - Use SPARQL package to submit query and save results to a data frame
  qd <- SPARQL(endpoint,query)
  df <- qd$results
  df
}

shinyServer(function(input, output) {

  output$table = DT::renderDataTable({
    datatable(protein_to_pathway(input$proteins))
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
