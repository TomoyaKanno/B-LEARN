scientificSliderUI <- function(id, label = "Cutoff value:", min_exp = 3, max_exp = 10) {
    ns <- NS(id)
    
    # Create formatted choices because we had to in order make ShinyWidgets work.
    # I know this looks awful, it is. But it works.
    slider_choices <- paste0("<span>1e-", sprintf("%d", min_exp:max_exp), "</span>")
    
    sliderTextInput(ns("value"), 
                label = label,
                choices = slider_choices,
                selected = slider_choices[1],
                grid = TRUE,
                hide_min_max = TRUE)
}

scientificSliderServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        # Extract numeric value from HTML string
        # Yes this is cursed. Don't blame me. Blame ShinyWidgets.
        # I mean ok you can blame me a little.
        reactive({
            as.numeric(gsub("</?span>", "", input$value))
        })
    })
}