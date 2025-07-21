# B-LEARN Portal

[![DOI](https://zenodo.org/badge/1022341370.svg)](https://doi.org/10.5281/zenodo.16120941)

**B**-cell **L**ymphoma Gene **E**xpression **A**ctivator and **R**epressor **N**etwork

## Overview

B-LEARN Portal is a comprehensive Shiny web application for analyzing regulatory networks derived from 47 CRISPR/Cas9 screens in human B lymphocytes. The application provides interactive tools for exploring gene regulation patterns, target identification, and regulatory network analysis.

For issues or questions contact tomoya.kanno97@gmail.com

## Getting Started

[work in progress]

## Features

### Core Modules

#### 1. Regulator Analysis
- Select from 47 CRISPR screens to display identified regulators sorted by most significant
- Separate analysis panels for repressors and activators
- MAGeCK cutoff sliders to customize significance score values
- Biological filtering (Everything, Kinases, Phosphatases, Transcription Factors, Epigenetic factors)
- Interactive bar plots showing top regulators
- STRING protein-protein interaction network visualization

#### 2. Target Analysis
- Search any regulator to see which of the 47 screens detected it as a regulator (Opposite direction from Regulator tab)
- Identify targets that a regulator of interest represses or activates
- MAGeCK cutoff controls for effect significance
- Interactive visualizations of regulatory relationships
- STRING interactome networks for target validation

#### 3. Corregulation 
- Compare 2-4 CRISPR screens simultaneously to find shared regulators
- Dynamic Venn diagrams showing regulatory overlaps
- Separate analysis for shared repressors and activators
- Interactive data tables for intersection exploration
- Customizable intersection selection (all genes, pairwise, three-way combinations)
- Biological function filtering and MAGeCK cutoffs

#### 4. Regulatory Network Analysis 
- Multi-regulator network visualization
- Pick any number of regulators of interest to see which shared genes it regulates
- Arrow color and thickness signify activtation/repression and significance of the effect
- Comprehensive view of regulatory relationships
- Visual network analysis with node and edge interactions

## Installation

### Prerequisites
- R (version 4.5 tested.)
- Required R packages:
  ```r
  install.packages(c(
    "shiny", "ggplot2", "tidyr", "readr", 
    "shinyWidgets", "rlang", "shinycssloaders", 
    "yaml", "DT", "plotly", "visNetwork"
  ))
  ```

### Setup
1. Clone this repository
2. Ensure `B-LEARN_Data.RDS` is in the project root (Built using Data_Builder.R)
3. Configure `theme.yml` for custom styling

## Usage

### Running the Application
```r
shiny::runApp()
```
Recommended to use RStudio where native shiny UI is present. Optionally, use RStudio to push the app to shinyapps.io or a local Posit Connect server.

### Data Requirements
- Primary dataset: `B-LEARN_Data.RDS`
- Theme configuration: `theme.yml`
- Logo image: `www/Logo3.png`

### Navigation
The application features a tabbed interface with four main analysis modules:
1. **Regulator**: Individual gene regulatory analysis
2. **Target**: Cross-gene target analysis  
3. **Coregulator**: Shared regulatory mechanisms
4. **Regulatory Network**: Comprehensive network analysis

## Data Source

Based on comprehensive CRISPR/Cas9 screening data from 47 screens in human B lymphocytes, providing insights into B-cell lymphoma gene regulation mechanisms.

## File Structure

```
B-LEARN_Portal/
├── app.R                          # Main application file
├── B-LEARN_Data.RDS              # Primary dataset
├── theme.yml                      # Theme configuration
├── README.md                      
├── www/
│   └── Logo3.png                 # Application logo
└── modules/
    ├── regulator_module.R        # Regulator analysis module
    ├── target_module.R           # Target analysis module
    ├── coregulator_module.R      # Coregulator analysis module
    ├── regulatory_network_module.R # Network analysis module
    ├── scientific_slider_module.R # Custon Shiny UI Utility
    ├── biological_filter_module.R # Biological filtering Utility
    └── utils.R                   # Other Utility functions
```

## Citation

If you use B-LEARN Portal in your research, please cite:

[![DOI](https://zenodo.org/badge/1022341370.svg)](https://doi.org/10.5281/zenodo.16120941)

## Author

Tomoya Kanno (tomoya.kanno97@gmail.com)

## License

[TBD]
