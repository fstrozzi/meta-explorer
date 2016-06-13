
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)

shinyUI(fluidPage(

  # Application title
  titlePanel("MetaExplorer"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      textInput("proteins",
                  "List of proteins:",
                  value = "Protein7660",
                  placeholder = "Paste a list of protein ids")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      DT::dataTableOutput('table')
    )
  )
))
