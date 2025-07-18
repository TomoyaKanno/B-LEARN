biologicalFilterUI <- function(id) {
  ns <- NS(id)
  
  radioButtons(ns("filterbutton"), "Biological function:",
      choiceNames = c("Everything", "Kinases", "Phosphatases", 
                    "Transcription Factors", "Epigenetic factors"),
      choiceValues = c("Everything", "kinases", "phosphatases", 
                     "HumanTFs", "epigenetic_factors"),
      selected = "Everything")
}

biologicalFilterServer <- function(id) {
  moduleServer(id, function(input, output, session) {
      return(reactive({ input$filterbutton }))
  })
}