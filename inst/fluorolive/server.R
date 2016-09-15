library(dplyr)
library(fluoroscripts)
library(shiny)

js <- require('shinyjs')

#if (!js) {
#  message <- "shinyjs package required to run this app. You can install it
#              directly from github by using install_github("daattali/shinyjs").
#              You'll need to install devtools in order to use install_github."
#  error(gsub('\n\\s+',' ', message))
#}

# Loading protein spectra
data(protein_spectra)

# Define server logic required to draw function
shinyServer(function(input, output) {

  # Expression that generates the plot
  output$spectra <- renderPlot({
    if (!is.null(input$prot1) & !is.null(input$prot2)) {

      d1 <- filter(protein.spectra, protein == input$prot1)
      d2 <- filter(protein.spectra, protein == input$prot2)

      w1 <- d1$wavelength
      ex1 <- d1$excitation/max(d1$excitation)
      em1 <- d1$emission/max(d1$emission)

      w2 <- d2$wavelength
      ex2 <- d2$excitation/max(d2$excitation)
      em2 <- d2$emission/max(d2$emission)

      plot(w1, ex1, col = 'black', type = 'l', lty = 2, lwd = 2)
      lines(w1, em1, col = 'black', lwd = 2)

      lines(w2, ex2, col = 'red', lty = 2, lwd = 2)
      lines(w2, em2, col = 'red', lwd = 2)

      xlab('Wavelength (nm)')
      ylab('Relative intensity')

    } else {
      plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '')
    }
  })
})
