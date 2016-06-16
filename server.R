
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
library(ggplot2)

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
             
             SELECT ?protein ?ko ?unipathway ?pathwayName ?biocyc ?goLabel
             WHERE
             {
             ?protein rdf:type up:Protein .
             OPTIONAL {
             ?protein rdfs:seeAlso ?ko .
             ?ko up:database <http://purl.uniprot.org/database/KO> .
             }
             OPTIONAL {
             ?protein up:annotation ?node .
             ?node rdf:type up:Pathway_Annotation .
             ?node rdfs:seeAlso ?unipathway .
             ?unipathway rdfs:label ?pathwayName . }
             OPTIONAL {
             ?protein rdfs:seeAlso ?biocyc . 
             ?biocyc up:database <http://purl.uniprot.org/database/BioCyc>
             }
             OPTIONAL {
             ?protein up:classifiedWith ?go .
	           ?go up:database <http://purl.uniprot.org/database/go> .
             ?go rdfs:label ?goLabel .
             }
             VALUES ?protein {@{id}}
             }
             ")

  # Step 2 - Use SPARQL package to submit query and save results to a data frame
  qd <- SPARQL(endpoint,query)
  df <- qd$results
  return(df)
}

uri2url = function(x) {
  if (is.na(x)) {
    return(NA)
  }
  else {
    x = gsub('<','',x)
    x = gsub('>','',x)
    name = tail(unlist(strsplit(x,"/")),n=1)
    return(toString(tags$a(href=x,name)))
  }
}

create_urls = function(dt) {
  dt[,protein:=uri2url(protein), by = 1:nrow(dt)]
  dt[,ko:=uri2url(ko), by = 1:nrow(dt)]
  dt[,unipathway:=uri2url(unipathway), by = 1:nrow(dt)]
  dt[,biocyc:=uri2url(biocyc), by = 1:nrow(dt)]
  return(dt)
}

unique_values = function(values) {
  paste0(unique(unlist(strsplit(values,";"))),collapse=",")
}

collapse_triplets = function(dt) {
  dt[,lapply(.SD,unique_values),by=c("protein","ko")]
}

query_endpoint = function(ids) {
  values = paste(unlist(ids),collapse=" ")
  results = as.data.table(protein_to_pathway(values))
  return(results)
}

shinyServer(function(input, output) {
  
  results = reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Querying the endpoints...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    proteins = unlist(strsplit(input$proteins,"\n"))
    if (length(proteins) >= 1) {
      
      proteins = unlist(lapply(proteins,function(x) paste("uniprotkb:",x,sep="")))
      splitted_ids = split(proteins, ceiling(seq_along(proteins)/200))
      dt <- data.frame(protein=character(),
                       ko=character(),
                       unipathway=character(),
                       pathwayName=character(),
                       biocyc=character(),
                       go=character(),
                       stringsAsFactors=FALSE)
      for(i in 1:length(splitted_ids)) {
        results = query_endpoint(splitted_ids[i])
        dt = rbindlist(list(dt,results))
      }
      if (length(dt) != 0) {
        dt[is.na(dt)] = ""
        dt = create_urls(dt)
        return(collapse_triplets(dt))
      }
    }
    
  })
  
  output$table = DT::renderDataTable({
    datatable(results(),escape = F,colnames=c("Protein ID","KEGG Ortholog","UniPathway IDs","Pathways Names","BioCyc IDs","Gene Ontologies"))
  })
  output$pathway = renderPlot({
    # Get the results
    data = results()
    pathways = as.data.frame(unlist(strsplit(data$pathwayName[!is.na(data$pathwayName)],",")))
    colnames(pathways) = c("pathwayName")
    pathways_freq = data.frame(table(pathways$pathwayName))
    pathways_filtered = subset(pathways_freq,Freq > 1)
    ggplot(data=pathways_filtered, aes(x=Var1,y=Freq)) + 
      geom_bar(stat="identity",fill="royalblue2",colour="black") + 
      coord_flip() + 
      theme_bw() + 
      xlab("Pathways name") +
      ggtitle("Pathways") +
      theme(axis.text.y = element_text(colour="grey20",size=10))
  },height=800,width=600)

  
  output$go = renderPlot({
    # Get the results
    data = results()
    go = as.data.frame(unlist(strsplit(data$go[!is.na(data$go)],",")))
    colnames(go) = c("go")
    go_freq = data.frame(table(go$go))
    go_filtered = subset(go_freq,Freq > 1)
    ggplot(data=go_filtered, aes(x=Var1,y=Freq)) + 
      geom_bar(stat="identity",fill="tomato2",colour="black") + 
      coord_flip() + 
      theme_bw() + 
      xlab("Gene Ontologies") + 
      ggtitle("Gene Ontology") +
      theme(axis.text.y = element_text(colour="grey20",size=10))
  },height=800,width=600)
  
})
