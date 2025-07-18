regulatorUI <- function(id, app_theme) {
    ns <- NS(id)
    
    tabPanel("Regulators",
        # Gene selection panel
        fluidRow(style = "margin-top: 20px;",
            column(12,
                wellPanel(
                    style = app_theme$styles$control_panel,
                    h3("Gene Selection", style = app_theme$headers$base),
                    selectizeInput(ns("GeneName1"),
                        label = "Select 1 of 47 screened genes to display its identified regulators:",
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
                        column(6, scientificSliderUI(ns("cutoff1"), "Cutoff value:")),
                        column(6, biologicalFilterUI(ns("filter1")))
                    ),
                    uiOutput(ns("DynamicSubHeader1")),
                    div(
                        style = 'max-height:500px; overflow-y: scroll; position: relative',
                        withSpinner(plotlyOutput(ns("ScreenBar1")), 
                                  color = app_theme$colors$repressor,
                                  type = 8)
                    ),
                    br(),
                    h3("STRING Interactome", style = app_theme$headers$repressor),
                    withSpinner(visNetworkOutput(ns("STRINGdb1")),
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
                        column(6, scientificSliderUI(ns("cutoff2"), "Cutoff value:")),
                        column(6, biologicalFilterUI(ns("filter2")))
                    ),
                    uiOutput(ns("DynamicSubHeader2")),
                    div(
                        style = 'max-height:500px; overflow-y: scroll; position: relative',
                        withSpinner(plotlyOutput(ns("ScreenBar2")), 
                                  color = app_theme$colors$activator,
                                  type = 8)
                    ),
                    br(),
                    h3("STRING Interactome", style = app_theme$headers$activator),
                    withSpinner(visNetworkOutput(ns("STRINGdb2")),
                              color = app_theme$colors$activator,
                              type = 8)
                )
            )
        )
    )
}

regulatorServer <- function(id, datalist) {
    moduleServer(id, function(input, output, session) {
      
        # Custom sci notation sliders
        cutoff1 <- scientificSliderServer("cutoff1")
        cutoff2 <- scientificSliderServer("cutoff2")

        # Custom filter values
        filter1 <- biologicalFilterServer("filter1")
        filter2 <- biologicalFilterServer("filter2")

        updateSelectizeInput(session, "GeneName1", 
                           choices = datalist$target, 
                           server = TRUE, 
                           selected = "LYN")
        
        # Creating a reactive value that only updates with valid selection
        # Doing this to prevent blank outs
        currentGene <- reactiveVal("LYN")
        observeEvent(input$GeneName1, {
            if (!is.null(input$GeneName1) && input$GeneName1 != "") {
                currentGene(input$GeneName1)
            }
        }, ignoreInit = TRUE)  # ignore the initial NULL->CD79A transition


        filteredData1 <- reactive({
            req(currentGene())
            filtered <- filterGeneData(
                gene = currentGene(),
                column_name = "target",
                TFtype = "Repr",
                cutoff = cutoff1(),
                filtercategory = filter1(),
                datalist = datalist
            )
            return(filtered)
        })
        
        filteredData2 <- reactive({
            req(currentGene())
            filtered <- filterGeneData(
                gene = currentGene(),
                column_name = "target",
                TFtype = "Acti",
                cutoff = cutoff2(),
                filtercategory = filter2(),
                datalist = datalist
            )
            return(filtered)
        })
        
        output$DynamicSubHeader1 <- renderUI({
            gene <- currentGene()
            h3(paste("Top Repressors for", gene, 
                    "(n=", nrow(filteredData1()), ")"))
        })
        
        output$DynamicSubHeader2 <- renderUI({
            gene <- currentGene()
            h3(paste("Top Activators for", gene, 
                    "(n=", nrow(filteredData2()), ")"))
        })
        
        output$ScreenBar1 <- renderPlotly({
            plotGeneScoreBar(filteredData1(), yaxisname = "effect", colorfill = "#D8222A")
        })
        
        output$ScreenBar2 <- renderPlotly({
            plotGeneScoreBar(filteredData2(), yaxisname = "effect", colorfill = "#28A1D0")
        })
        
        output$STRINGdb1 <- renderVisNetwork({
            getSTRINGVis(filteredData1(), 
                        center_gene = currentGene(),
                        resulttype = "effect",
                        colorfill = "#D8222A",
                        datalist = datalist)
        })
        
        output$STRINGdb2 <- renderVisNetwork({
            gene <- currentGene()
            getSTRINGVis(filteredData2(), 
                        center_gene = currentGene(),
                        resulttype = "effect",
                        colorfill = "#28A1D0",
                        datalist = datalist)
        })
    })
} 