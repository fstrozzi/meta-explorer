
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)
library(shinydashboard)


dashboardPage(
  dashboardHeader(title="Meta Explorer"),
  dashboardSidebar(),
  dashboardBody(
  # Application title
  
  
  # Sidebar with a slider input for number of bins
  fluidRow(
    column(12,
      tags$textarea(id="proteins", rows=3, cols=40, "P0A6P9")
      # textInput("protein",
      #             "List of proteins:",
      #             value = "P0A6P9",
      #             placeholder = "Paste a list of protein ids")
    )
  ),
  fluidRow(
    # Show a plot of the generated distribution
    column(12,
      DT::dataTableOutput('table')
    ),
    column(6,
      plotOutput('pathway')
    ),
    column(6,
      plotOutput('go')
    )
  ),
  fluidRow(
    column(12,
      plotOutput("biocyc")
    )
  )
))


