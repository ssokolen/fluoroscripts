library(fluoroscripts)
library(shiny)

data(protein_spectra)
proteins <- unique(protein.spectra$protein)

# Logistic function
shinyUI(fluidPage(

  # Application title
  titlePanel("Protein overlap"),

  # Sidebar
  sidebarLayout(
    sidebarPanel(
      selectInput('prot1', 'Protein 1', proteins, selected = NULL, 
                  multiple = FALSE, selectize = TRUE),
      selectInput('prot2', 'Protein 2', proteins, selected = NULL, 
                  multiple = FALSE, selectize = TRUE)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("spectra", width=800, height=500)
    )
  )
))
