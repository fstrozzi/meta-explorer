
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)
library(SPARQL)
library(GetoptLong)
library(data.table)

protein_to_pathway <- function(id){
  # Step 1 - Set up preliminaries and define query
  # Define the data.gov endpoint
  endpoint <- "http://sparql.uniprot.org/sparql"
  
  # create query statement
  query <-qq("
             PREFIX up:<http://purl.uniprot.org/core/> 
             PREFIX keywords:<http://purl.uniprot.org/keywords/> 
             PREFIX uniprotkb:<http://purl.uniprot.org/uniprot/> 
             PREFIX taxon:<http://purl.uniprot.org/taxonomy/> 
             PREFIX ec:<http://purl.uniprot.org/enzyme/> 
             PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
             PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
             PREFIX skos:<http://www.w3.org/2004/02/skos/core#> 
             PREFIX owl:<http://www.w3.org/2002/07/owl#> 
             PREFIX bibo:<http://purl.org/ontology/bibo/> 
             PREFIX dc:<http://purl.org/dc/terms/> 
             PREFIX xsd:<http://www.w3.org/2001/XMLSchema#> 
             PREFIX faldo:<http://biohackathon.org/resource/faldo#> 
             
             
             SELECT ?ko ?unipathway ?pathwayName
             WHERE
             {
             uniprotkb:@{id} rdfs:seeAlso ?ko .
             ?ko up:database <http://purl.uniprot.org/database/KO> .
             uniprotkb:@{id} up:annotation ?node .
             ?node rdf:type up:Pathway_Annotation .
             ?node rdfs:seeAlso ?unipathway .
             ?unipathway rdfs:label ?pathwayName .
             }")

  # Step 2 - Use SPARQL package to submit query and save results to a data frame
  qd <- SPARQL(endpoint,query)
  df <- qd$results
  df
}

shinyServer(function(input, output) {

  output$table = DT::renderDataTable({
    print(input$protein)
    dt = as.data.table(protein_to_pathway("P0A6P9"))
    dt[,ko:=toString(tags$a(href=ko,ko)), by = 1:nrow(dt)]
    datatable(dt)
  })


})
