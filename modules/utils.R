library(ggVennDiagram)
library(plotly)
library(visNetwork)
library(dplyr)
library(httr)
library(readr)
library(stats)

# Function for filtering data
filterData <- function(gene, column_name, TFtype, cutoff, filtercategory, datalist) {
    column_sym = as.symbol(column_name)
    
    result <- datalist[[TFtype]] %>%
        filter(!!column_sym == gene) %>%
        filter(score < cutoff) %>%
        mutate(logscore = -log(score)) %>%
        arrange(logscore)
    
    if (filtercategory != "Everything") {
        result <- result %>%
            filter(target %in% datalist[[filtercategory]] | effect %in% datalist[[filtercategory]])
    }
    
    result %>%
        mutate(
            target = factor(target, levels = unique(target)),
            effect = factor(effect, levels = unique(effect)),
            text = paste("score: ", score)
        )
}

# NEW filtering function to cover both barpltos and STRINGvis
filterGeneData <- function(gene, column_name, TFtype, cutoff, filtercategory, datalist) {
    column_sym = as.symbol(column_name)
    
    result <- datalist[[TFtype]] %>%
        filter(!!column_sym == gene) %>%
        filter(score < cutoff)
    
    if (filtercategory != "Everything") {
        result <- result %>%
            filter(effect %in% datalist[[filtercategory]])
    }
    
    return(result)
}

# Function for bar plots version 2
plotGeneScoreBar <- function(filtered_data, yaxisname, colorfill) {

    # Log transform and arrange data
    plot_data <- filtered_data %>%
        mutate(
            logscore = -log(score),
            text = paste("score:", score)
        ) %>%
        arrange(logscore) %>%
        mutate(  # Arranging by size
            target = factor(target, levels = unique(target)),
            effect = factor(effect, levels = unique(effect))
        )

    
    fig <- plot_ly(
        plot_data, 
        x = ~logscore, 
        y = ~plot_data[[yaxisname]], 
        type = 'bar', 
        orientation = 'h', 
        color = I(colorfill), 
        text = ~text, 
        hoverinfo = "text",
        height = if(nrow(plot_data) > 15) 35*nrow(plot_data) else NULL
    )
    fig <- fig %>% layout(showlegend = FALSE,
        yaxis = list(title = ''))
    
    return(fig)
}

# Function for STRING visualization
# Updated getSTRINGVis function with API fallback
getSTRINGVis <- function(filtered_data, center_gene, resulttype, colorfill, datalist) {
    # Extract unique gene to visualize
    result <- filtered_data %>%
        select(as.symbol(resulttype)) %>%
        pull() %>%
        unique()
    
    # Add the gene we're looking at
    result <- c(center_gene, result)
    
    map <- tibble(genes = result) %>%
        mutate(ID_upper = toupper(genes)) %>%
        inner_join(datalist$STRINGIDmap, by = 'ID_upper') %>%
        select(-ID_upper)
    
    # Try to get pre-calculated network first
    network_data <- datalist$string_networks[[center_gene]]
    
    # If no pre-calculated network exists, fall back to API
    if(is.null(network_data)) {
        message(sprintf("No pre-calculated network found for %s. Falling back to API...", gene))
        
        # Add small delay to respect API rate limits
        Sys.sleep(1)
        
        tryCatch({
            base_url <- "https://string-db.org/api"
            query_url_net <- paste0(base_url, "/tsv/network?",
                                  "identifiers=", paste(map$STRING_id, collapse = "%0d"),
                                  "&species=9606",
                                  "&caller_identity=TomoyaScript")
            
            APIres = GET(query_url_net)
            
            if(status_code(APIres) == 200) {
                network_data <- content(APIres, "text", encoding = "UTF-8") %>%
                    read_tsv() %>%
                    select(stringId_A = stringId_A,
                           stringId_B = stringId_B,
                           score = score)
            } else {
                warning(sprintf("API request failed for %s. Status code: %d", 
                              gene, status_code(APIres)))
                return(visNetwork(data.frame(), data.frame()))
            }
        }, error = function(e) {
            warning(sprintf("Error fetching network data for %s: %s", gene, e$message))
            return(visNetwork(data.frame(), data.frame()))
        })
    }
    
    # Process network data (whether from cache or API)
    Net_DF <- network_data %>%
        inner_join(map, by = c("stringId_A" = "STRING_id")) %>%
        rename(nodeA_name = genes) %>%
        inner_join(map, by = c("stringId_B" = "STRING_id")) %>%
        rename(nodeB_name = genes)
    
    # If no interactions found after filtering
    if(nrow(Net_DF) == 0) {
        return(visNetwork(data.frame(), data.frame()))
    }
    
    edges <- Net_DF %>%
        rename(from = stringId_A,
               to = stringId_B) %>%
        mutate(
            value = score,
            title = paste("<b>",nodeA_name, "to", nodeB_name, "</b><br>STRING Score:",score),
            color = colorfill
        )
    
    # Create nodes dataframe with special highlighting for center gene
    nodes <- map %>%
        rename(id = STRING_id) %>%
        mutate(
            label = genes,
            shape = "ellipse",
            group = if_else(genes == center_gene, "center", "others"),
            borderWidth = if_else(genes == center_gene, 30, 1),
            shadow = genes == center_gene  # Add shadow only to center gene
        )
    
    visNetwork(as.data.frame(nodes), edges) %>%
        visGroups(
            groupname = "center", 
            color = colorfill, 
            shape = "ellipse"
        ) %>%
        visGroups(
            groupname = "others", 
            color = colorfill, 
            shape = "ellipse"
        ) %>%
        visNodes(font = list(color = "white")) %>%
        visPhysics(barnesHut = list(springConstant = 0)) %>%
        visEdges(smooth = FALSE) %>%
        visOptions(
            highlightNearest = list(
                enabled = TRUE, 
                degree = 1,
                hover = FALSE,
                algorithm = "hierarchical"
            ),
            nodesIdSelection = TRUE
        ) %>%
        visInteraction(
            hover = TRUE, 
            selectConnectedEdges = FALSE,
            hideEdgesOnDrag = FALSE,
            multiselect = FALSE,
            zoomView = FALSE, # Disable scroll zoom
            navigationButtons = TRUE
        )
}

# Fancy Venn diagrams with external controls
getVennFiltered <- function(genes, TFtype, coloring, cutoff, filter, datalist) {
    # Get filtered data for each gene
    all_sets <- lapply(genes, function(gene) {
        filtered <- filterData(gene, "target", TFtype, cutoff, filter, datalist)
        filtered$effect  # Get just the effects
    })
    
    # Name the sets with gene names
    names(all_sets) <- genes
    
    # Create the Venn diagram with custom styling
    ggVennDiagram(all_sets, 
                 set_color = "black",
                 edge_size = 1) + 
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_gradient(
            low = "white",
            high = coloring,
            # This makes the gradient start further up from zero
            # This limit trick works. I don't know why. Don't change it.
            limits = c(-0.001, NA),  # Slight negative value ensures true zero is white
            guide = "none"
        ) +
        theme(legend.position = "none")
}

# Complex intersection finder for the data table under Venn.
get_intersection_data <- function(genes, intersection_choice, TFtype, cutoff_value, datalist) {
    # Validate inputs
    if(length(genes) < 2) return(data.frame())
    if(is.null(intersection_choice)) return(data.frame())
    
    # Get data source based on type
    data_source <- datalist[[TFtype]]
    
    # Get the genes for the selected intersection
    selected_combination <- if(intersection_choice == "All genes") {
        genes
    } else {
        strsplit(intersection_choice, " and ")[[1]]
    }
    
    # Get data for selected genes
    gene_data <- lapply(selected_combination, function(gene) {
        data_source %>%
            filter(target == gene,
                   score < cutoff_value) %>%
            select(effect, score)
    })
    
    # Find shared effects for selected intersection
    # Using lapply because gene_data is a list of dataframes not a single dataframe
    shared_effects <- Reduce(intersect, 
                           lapply(gene_data, function(x) x$effect))
    
    # Get effects that are NOT in other genes (if not showing all)
    if(intersection_choice != "All genes") {
        other_genes <- setdiff(genes, selected_combination)
        if(length(other_genes) > 0) {
            other_effects <- unique(unlist(lapply(other_genes, function(gene) {
                data_source %>%
                    filter(target == gene,
                           score < cutoff_value) %>%
                    pull(effect)
            })))
            shared_effects <- setdiff(shared_effects, other_effects)
        }
    }
    
    # Create result table
    result <- data.frame(
        Regulator = shared_effects
    )
    
    # Add score columns for the selected intersection genes
    for(gene in selected_combination) {
        scores <- gene_data[[which(selected_combination == gene)]]
        result[[gene]] <- scores$score[match(shared_effects, scores$effect)]
    }
    
    result
}

# Function for network visualization
plotRegulatoryNetwork <- function(genes, datalist) {
    result <- datalist$ActiRepr %>%
        filter(effect %in% genes)
    
    edges <- result %>%
        rename(from = effect) %>%
        rename(to = target) %>%
        mutate(
            value = -log(score),
            arrows = "to",
            color = if_else(type == "Activator",'#28A1D0','#D8222A'),
            title = paste("<b>",from, "to", to, "</b><br>Score:",score,"<br>Type:",type)
        )
    
    nodes <- result %>%
        select(effect, target) %>%
        pivot_longer(everything()) %>%
        distinct(value) %>%
        rename(id = value) %>%
        mutate(
            label = id,
            group = if_else(label %in% genes, "inputs", "others"),
            fixed = if_else(group == "inputs", TRUE,FALSE)
        )
    
    visNetwork(as.data.frame(nodes), edges) %>%
        visGroups(groupname = "inputs", color = "yellow", shape = "elipse") %>%
        visGroups(groupname = "others", color = "lightgrey", shape = 'circle') %>%
        visOptions(
            highlightNearest = list(
                enabled = T, 
                degree = 1,
                hover = FALSE,
                algorithm = "hierarchical")) %>%
        visInteraction(zoomView = FALSE) %>%
        visPhysics(
          solver = "barnesHut",
          barnesHut = list(
            gravitationalConstant = -2000,
            springLength = 400,
            damping = 2
          )
        ) # %>%
      # visEvents(stabilizationIterationsDone = "function() { this.physics.disable(); }")
}