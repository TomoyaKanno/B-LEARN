# Load required packages
library(shiny)
library(ggplot2)
library(tidyr)
library(readr)
library(shinyWidgets)
library(rlang)
library(shinycssloaders)
library(yaml)
library(DT)

# Load theme
app_theme <- yaml::yaml.load_file("theme.yml")

# Load data
datalist <- readRDS("B-LEARN_Data.RDS")

# Source all module files
source("modules/regulator_module.R")
source("modules/target_module.R")
source("modules/coregulator_module.R")
source("modules/regulatory_network_module.R")
source("modules/scientific_slider_module.R")
source("modules/biological_filter_module.R")
source("modules/utils.R")

ui <- fluidPage(
    # Custom title bar HTML elements
    tags$head(
      tags$style(HTML("
              .main-title { 
                  font-size: 28px;
                  margin-bottom: 5px;
              }
              .expanded-name {
                  font-size: 16px;
                  color: #666;
                  margin-bottom: 5px;
              }
              .header-container {
                  display: flex;
                  align-items: center;
                  gap: 20px; 
              }
              .logo {
                  height: 100px; 
                  width: auto;
              }
              .text-container {
                  flex: 1; 
              }
          "))
    ),
    
    # Header using native Shiny elements with minimal customization
    div(
      class = "container",  
      style = "padding-bottom: 10px;",
      div(class = "header-container",
          img(src = "Logo3.png", class = "logo"),
          div(class = "text-container",
              h1(class = "main-title", "B-LEARN Portal"),
              p(class = "expanded-name", 
                HTML(paste(
                    "<strong>B</strong>-cell ",
                    "<strong>L</strong>ymphoma Gene ",
                    "<strong>E</strong>xpression ",
                    "<strong>A</strong>ctivator and ",
                    "<strong>R</strong>epressor ",
                    "<strong>N</strong>etwork"
                ))),
              p("Comprehensive regulatory network analysis from 47 CRISPR/Cas9",
                "screens in human B lymphocytes")
          )
      )
    ),
      
    tabsetPanel(
        id = "tabs",
        type = "pills",
        regulatorUI("perGene", app_theme),
        targetUI("acrossGenes", app_theme),
        coregulatorUI("sharedGenes", app_theme),
        regulatoryNetworkUI("regNetwork", app_theme)
    )
)

# Server
server <- function(input, output, session) {
    regulatorServer("perGene", datalist)
    targetServer("acrossGenes", datalist)
    coregulatorServer("sharedGenes", datalist)
    regulatoryNetworkServer("regNetwork", datalist)
}

shinyApp(ui, server)