regulatoryNetworkUI <- function(id, app_theme) {
    ns <- NS(id)
    
    tabPanel("Regulatory Network",
        # Controls in themed panel
        fluidRow(style = "margin-top: 20px;",
            column(12,
                wellPanel(
                    style = app_theme$styles$control_panel,
                    h3("Select Genes", style = app_theme$headers$base),
                    selectizeInput(ns("GeneNames2"),
                        label = "Type Gene here",
                        choices = NULL,
                        multiple = TRUE
                    )
                )
            )
        ),        
        fluidRow(
            column(12,
                visNetworkOutput(ns("visNetOut"), width = "100%", height = "700px")
            )
        )
    )
}

regulatoryNetworkServer <- function(id, datalist) {
    moduleServer(id, function(input, output, session) {
        # Initial reactive value
        currentGenes <- reactiveVal(c("TCF3", "ID3"))

        updateSelectizeInput(session, 'GeneNames2',
            choices = datalist$effect,
            server = TRUE,
            selected = c("TCF3", "ID3")
        )

        # Update currentGenes only when there's a valid selection
        observeEvent(input$GeneNames2, {
            if (!is.null(input$GeneNames2) && length(input$GeneNames2) > 0) {
                currentGenes(input$GeneNames2)
            }
        }, ignoreInit = TRUE)
        
        output$visNetOut <- renderVisNetwork({
            plotRegulatoryNetwork(genes = currentGenes(), datalist)
        })
    })
}