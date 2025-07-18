coregulatorUI <- function(id, app_theme) {
  ns <- NS(id)
  
  tabPanel("Co-regulation",
      # Gene selection panel
      fluidRow(style = "margin-top: 20px;",
      column(12,
          wellPanel(
              style = app_theme$styles$control_panel,
              h3("Gene Selection", style = app_theme$headers$base),
              fluidRow(
                  column(3, 
                      numericInput(ns("num_genes"), 
                          "Number of genes to compare (max4):", 
                          value = 2, min = 2, max = 4)
                  )
              ),
              uiOutput(ns("gene_selectors"))
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
                      column(6, scientificSliderUI(ns("cutoffRepr"), "Cutoff value:")),
                      column(6, biologicalFilterUI(ns("filterRepr")))
                  ),
                  uiOutput(ns("DynamicSubHeader1")),
                  plotOutput(ns("Venn1"))
              )
          ),
          # Activator panel
          column(6,
              wellPanel(
                  style = app_theme$styles$activator_panel,
                  h3("Activator Analysis", style = app_theme$headers$activator),
                  fluidRow(
                      column(6, scientificSliderUI(ns("cutoffActi"), "Cutoff value:")),
                      column(6, biologicalFilterUI(ns("filterActi")))
                  ),
                  uiOutput(ns("DynamicSubHeader2")),
                  plotOutput(ns("Venn2"))
              )
          )
      ),

        # data table section

        # Add intersection selector
        fluidRow(
          column(12,
              wellPanel(
                  style = app_theme$styles$control_panel,
                  h3("Intersection Selection", style = app_theme$headers$base),
                  uiOutput(ns("intersection_selector"))
              )
          )
      ),

        # Actual data tables
        fluidRow(
          column(6,
              wellPanel(
                  style = app_theme$styles$repressor_panel,
                  h3("Shared Repressors Table", style = app_theme$headers$repressor),
                  DT::dataTableOutput(ns("repressor_table"))
              )
          ),
          column(6,
              wellPanel(
                  style = app_theme$styles$activator_panel,
                  h3("Shared Activators Table", style = app_theme$headers$activator),
                  DT::dataTableOutput(ns("activator_table"))
              )
          )
      )
  )
}

coregulatorServer <- function(id, datalist) {
  moduleServer(id, function(input, output, session) {
      # Initialize scientific slider modules
      cutoffRepr <- scientificSliderServer("cutoffRepr")
      cutoffActi <- scientificSliderServer("cutoffActi")

      # Initialize custom filter modules
      filterRepr <- biologicalFilterServer("filterRepr")
      filterActi <- biologicalFilterServer("filterActi")

      # Dynamic gene selection UI
      output$gene_selectors <- renderUI({
        req(input$num_genes)
        num_genes <- input$num_genes
        
        # Create a list of select inputs
        selectors <- lapply(1:num_genes, function(i) {
            column(12/num_genes,
                selectizeInput(
                    session$ns(paste0("GeneName", i)),
                    label = paste("Gene", i, ":"),
                    choices = datalist$target,
                    selected = ifelse(i == 1, "MAP3K7",
                              ifelse(i == 2, "PARP1", "")),
                    multiple = FALSE
                )
            )
        })
        
        fluidRow(selectors)
    })

        # Reactive values for selected genes
        selected_genes <- reactive({
          req(input$num_genes)
          genes <- character(0)
          for(i in 1:input$num_genes) {
              gene_input <- input[[paste0("GeneName", i)]]
              if(!is.null(gene_input) && gene_input != "") {
                  genes <- c(genes, gene_input)
              }
          }
          return(genes)
      })
      
    output$DynamicSubHeader1 <- renderUI({
      genes <- selected_genes()
      if(length(genes) >= 2) {
          h4(paste("Shared Repressors between", paste(genes, collapse = ", ")))
      }
    })
    
    output$DynamicSubHeader2 <- renderUI({
        genes <- selected_genes()
        if(length(genes) >= 2) {
            h4(paste("Shared Activators between", paste(genes, collapse = ", ")))
        }
    })
      
    output$Venn1 <- renderPlot({
      genes <- selected_genes()
      req(length(genes) >= 2)
      
      getVennFiltered(
          genes = genes,
          TFtype = "Repr", 
          coloring = "#D8222A", 
          cutoff = cutoffRepr(),
          filter = filterRepr(),
          datalist = datalist
      )
    })
    
    output$Venn2 <- renderPlot({
        genes <- selected_genes()
        req(length(genes) >= 2)
        
        getVennFiltered(
            genes = genes,
            TFtype = "Acti", 
            coloring = "#28A1D0", 
            cutoff = cutoffActi(),
            filter = filterActi(),
            datalist = datalist
        )
    })

      # Dynamic intersection selector
      output$intersection_selector <- renderUI({
        genes <- selected_genes()
        req(length(genes) >= 2)
        
        # Generate all possible intersection combinations
        combinations <- list()
        n <- length(genes)
        
        # Add "all genes" option
        combinations[["All genes"]] <- genes

        # Add pairwise combinations
        if(n > 2) {
          for(i in 1:(n-1)) {
              for(j in (i+1):n) {
                  key <- paste(genes[i], "and", genes[j])
                  combinations[[key]] <- genes[c(i,j)]
              }
          }
          
          # Add three-way combinations for 4-bubble Venns
          if(n > 3) {
              for(i in 1:(n-2)) {
                  for(j in (i+1):(n-1)) {
                      for(k in (j+1):n) {
                          key <- paste(genes[i], "and", genes[j], "and", genes[k])
                          combinations[[key]] <- genes[c(i,j,k)]
                      }
                  }
              }
          }
      }
        
        selectInput(
            session$ns("intersection_choice"),
            "Select intersection to view:",
            choices = names(combinations),
            selected = "All genes"
        )
    })


      # Get intersection data to use in data tables
      repressor_data <- reactive({
        get_intersection_data(
            genes = selected_genes(),
            intersection_choice = input$intersection_choice,
            TFtype = "Repr",
            cutoff_value = cutoffRepr(),
            datalist = datalist
        )
    })
    
    activator_data <- reactive({
        get_intersection_data(
            genes = selected_genes(),
            intersection_choice = input$intersection_choice,
            TFtype = "Acti",
            cutoff_value = cutoffActi(),
            datalist = datalist
        )
    })

      output$repressor_table <- DT::renderDataTable({
        req(selected_genes(), input$intersection_choice)
        DT::datatable(
            repressor_data(),
            options = list(
                scrollX = TRUE,  # Enable horizontal scrolling
                autoWidth = TRUE,  # Optimal column width
                pageLength = 10,   # Rows per page
                dom = 'frtip'      # DataTable controls (f=filter, r=processing, t=table, i=info, p=pagination)
            )
        )
    })
    
    output$activator_table <- DT::renderDataTable({
        req(selected_genes(), input$intersection_choice)
        DT::datatable(
            activator_data(),
            options = list(
                scrollX = TRUE,
                autoWidth = TRUE,
                pageLength = 10,
                dom = 'frtip'
            )
        )
    })

  })
} 