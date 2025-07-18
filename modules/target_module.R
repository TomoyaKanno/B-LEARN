targetUI <- function(id, app_theme) {
    ns <- NS(id)
    
    tabPanel("Targets",
        # Gene selection panel
        fluidRow(style = "margin-top: 20px;",
            column(12,
                wellPanel(
                    style = app_theme$styles$control_panel,
                    h3("Gene Selection", style = app_theme$headers$base),
                    selectizeInput(ns("GeneName3"),
                        label = "Type gene of interest to display in which of the 47 screens it was detected as a regulator:",
                        choices = NULL,
                        multiple = FALSE
                    )
                )
            )
        ),
        
        # Analysis panels
        fluidRow(
            # Repressor panel
            column(6,
                wellPanel(
                    style = app_theme$styles$repressor_panel,
                    h3("Repressor Analysis", style = app_theme$headers$repressor),
                    fluidRow(
                        column(12, scientificSliderUI(ns("cutoff3"), "Cutoff value:"))
                    ),
                    uiOutput(ns("DynamicSubHeader3")),
                    div(
                        style = 'max-height:500px; overflow-y: scroll; position: relative',
                        withSpinner(plotlyOutput(ns("ScreenBar3")), 
                                  color = app_theme$colors$repressor,
                                  type = 8)
                    ),
                    br(),
                    h3("STRING Interactome", style = app_theme$headers$repressor),
                    withSpinner(visNetworkOutput(ns("STRINGdb3")),
                              color = app_theme$colors$repressor,
                              type = 8)
                )
            ),
            # Activator panel
            column(6,
                wellPanel(
                    style = app_theme$styles$activator_panel,
                    h3("Activator Analysis", style = app_theme$headers$activator),
                    fluidRow(
                        column(12, scientificSliderUI(ns("cutoff4"), "Cutoff value:"))
                    ),
                    uiOutput(ns("DynamicSubHeader4")),
                    div(
                        style = 'max-height:500px; overflow-y: scroll; position: relative',
                        withSpinner(plotlyOutput(ns("ScreenBar4")), 
                                  color = app_theme$colors$activator,
                                  type = 8)
                    ),
                    br(),
                    h3("STRING Interactome", style = app_theme$headers$activator),
                    withSpinner(visNetworkOutput(ns("STRINGdb4")),
                              color = app_theme$colors$activator,
                              type = 8)
                )
            )
        )
    )
}

targetServer <- function(id, datalist) {
    moduleServer(id, function(input, output, session) {
        # Initialize scientific slider modules
        cutoff3 <- scientificSliderServer("cutoff3")
        cutoff4 <- scientificSliderServer("cutoff4")

        # Create reactive value for stable gene selection
        currentGene <- reactiveVal("PAX5") # initial value
        
        # Prevent blanks from resetting the app
        observeEvent(input$GeneName3, {
            if (!is.null(input$GeneName3) && input$GeneName3 != "") {
                currentGene(input$GeneName3)
            }
        }, ignoreInit = TRUE)

        updateSelectizeInput(session, "GeneName3",
        choices = datalist$effect,
        server = TRUE,
        selected = "PAX5")

        # Create filtered data reactives
        filteredData1 <- reactive({
            req(currentGene())
            filtered <- filterGeneData(
                gene = currentGene(),
                column_name = "effect",  
                TFtype = "Repr",
                cutoff = cutoff3(),
                filtercategory = "Everything", # No filter this tab
                datalist = datalist
            )
            return(filtered)
        })
        
        filteredData2 <- reactive({
            req(currentGene())
            filtered <- filterGeneData(
                gene = currentGene(),
                column_name = "effect",
                TFtype = "Acti",
                cutoff = cutoff4(),
                filtercategory = "Everything", # No filter this tab
                datalist = datalist
            )
            return(filtered)
        })
        
        output$DynamicSubHeader3 <- renderUI({
            req(currentGene())
            h2(paste(currentGene(), "as a Repressor of"))
        })
        
        output$DynamicSubHeader4 <- renderUI({
            req(currentGene())
            h2(paste(currentGene(), "as an Activator of"))
        })
        
        output$ScreenBar3 <- renderPlotly({
            plotGeneScoreBar(
                filteredData1(), 
                yaxisname = "target", 
                colorfill = "#D8222A")
        })
        
        output$ScreenBar4 <- renderPlotly({
            plotGeneScoreBar(
                filteredData2(), 
                yaxisname = "target", 
                colorfill = "#28A1D0")
        })
        
        output$STRINGdb3 <- renderVisNetwork({
            getSTRINGVis(filteredData1(),
                center_gene = currentGene(),
                resulttype = "target",
                colorfill = "#D8222A",
                datalist = datalist)
        })
        
        output$STRINGdb4 <- renderVisNetwork({
            getSTRINGVis(filteredData2(),
                center_gene = currentGene(),
                resulttype = "target",
                colorfill = "#28A1D0",
                datalist = datalist)
        })
    })
} 