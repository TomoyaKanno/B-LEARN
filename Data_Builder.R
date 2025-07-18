library(readxl)
library(dplyr)
library(readr)
library(STRINGdb)
library(httr)
library(tidyr)

# Original data loading code
acti_excel <- read_excel('data/all_activators NO_cellsize_Capitalized.xlsx')
repr_excel <- read_excel('data/all_repressors_NO_cellsize_Capitalized.xlsx')

Acti <- acti_excel %>%
  select(target, effect, score) %>%
  mutate(type = "Activator")

Repr <- repr_excel %>%
  select(target, effect, score) %>%
  mutate(type = "Repressor")

ActiRepr <- full_join(Acti, Repr)

# Added %>% sort() to make it alphabetical
target <- unique(ActiRepr$target) %>% na.omit() %>% sort()
effect <- unique(ActiRepr$effect) %>% na.omit() %>% sort()

# Import the gene categories lists
epigenetic_factors <- read_excel('data/EpiGenes_main.xlsx')
epigenetic_factors <- epigenetic_factors$...2
epigenetic_factors <- epigenetic_factors[-1]

kinases <- read_excel('data/Kincat_Hsap.08.02.xls', range = cell_cols("A"))
kinases <- kinases %>% 
  slice(-1) %>% 
  pull()

phosphatases <- read_excel('data/Phosphatases from aag1796_Tables S1 to S23.xlsx',
                           range = cell_cols("A"))
phosphatases <- phosphatases %>%
  slice(-(1:2)) %>%
  pull()

HumanTFs <- read_excel('data/human_TFs.xlsx', range = cell_cols("A"))
HumanTFs <- pull(HumanTFs)

# Use local STRINGdb for mapping
AllGenes <- unique(c(target, effect))
AllGenesDF <- data.frame(ID_upper = AllGenes)

string_db <- STRINGdb$new(version="12.0", 
                          species=9606, 
                          score_threshold=200,
                          network_type="full", 
                          link_data="combined_only",
                          input_directory='data')

STRINGIDmap <- string_db$map(AllGenesDF, "ID_upper", removeUnmappedRows = TRUE)

# Function to get STRING interactions via POST
get_string_interactions <- function(string_ids) {
  # Add small delay to respect API rate limits
  Sys.sleep(1)
  
  response <- POST(
    "https://string-db.org/api/tsv/network",
    body = list(
      identifiers = paste(string_ids, collapse = "%0d"),
      species = "9606",
      caller_identity = "TomoyaScript"
    ),
    encode = "form"
  )
  
  if(status_code(response) == 200) {
    content(response, "text", encoding = "UTF-8") %>%
      read_tsv() %>%
      select(stringId_A = stringId_A,
             stringId_B = stringId_B,
             score = score)
  } else {
    warning(sprintf("Failed to get interactions. Status code: %d", status_code(response)))
    NULL
  }
}

# Pre-calculate networks for each gene
message("\nPre-calculating STRING networks...")
string_networks <- list()
total_genes <- length(AllGenes)

for(i in seq_along(AllGenes)) {
  gene <- AllGenes[i]
  message(sprintf("Processing gene %d of %d: %s", i, total_genes, gene))
  
  # Get potential interaction partners from your data
  partners_acti <- ActiRepr %>%
    filter(target == gene | effect == gene) %>%
    select(target, effect) %>%
    unlist() %>%
    unique()
  
  # Get STRING IDs for these genes from our local mapping
  relevant_string_ids <- STRINGIDmap %>%
    filter(ID_upper %in% toupper(c(gene, partners_acti))) %>%
    pull(STRING_id)
  
  # Get interactions
  if(length(relevant_string_ids) >= 2) {
    network <- get_string_interactions(relevant_string_ids)
    if(!is.null(network)) {
      string_networks[[gene]] <- network
      message(sprintf("  Retrieved network with %d interactions", nrow(network)))
    } else {
      warning(sprintf("  Failed to retrieve network for %s", gene))
    }
  } else {
    message(sprintf("  Skipping %s - fewer than 2 mapped IDs", gene))
  }
}

# Add a summary at the end
successful_networks <- sum(!sapply(string_networks, is.null))
message(sprintf("\nSummary:\nProcessed %d genes\nSuccessfully retrieved %d networks", 
                total_genes, successful_networks))

# Create the complete datalist
datalist <- list(
  Acti = Acti,
  Repr = Repr,
  ActiRepr = ActiRepr,
  target = target,
  effect = effect,
  epigenetic_factors = epigenetic_factors,
  kinases = kinases,
  phosphatases = phosphatases,
  HumanTFs = HumanTFs,
  STRINGIDmap = STRINGIDmap,
  string_networks = string_networks
)

saveRDS(datalist, "B-LEARN_Data.RDS")