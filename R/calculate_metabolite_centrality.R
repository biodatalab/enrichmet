#' Calculate Metabolite Centrality
#'
#' Computes relative betweenness centrality for metabolites based on their 
#' co-membership in metabolic pathways.
#'
#' @param PathwayVsMetabolites A data frame containing pathways and their associated 
#'        metabolites. Must include columns 'Pathway' and 'Metabolites'.
#'
#' @return A data frame with metabolite names and their relative betweenness centrality scores.
#'
#' @examples
#' # Create example pathway data with valid KEGG metabolite IDs
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Glycolysis / Gluconeogenesis",
#'               "Citrate cycle (TCA cycle)",
#'               "Pentose phosphate pathway"),
#'   Metabolites = c("C00031,C00022,C00197,C00221,C00631,C01172",
#'                   "C00022,C00036,C00024,C00158,C00149,C00311",
#'                   "C00117,C00257,C00121,C00085,C00118")
#' )
#'
#' # Calculate centrality
#' centrality <- calculate_metabolite_centrality(PathwayVsMetabolites)
#' head(centrality)
#'
#' @importFrom dplyr filter
#' @importFrom tidyr unnest
#' @importFrom igraph graph_from_data_frame betweenness graph.adjacency
#' @export
calculate_metabolite_centrality <- function(PathwayVsMetabolites) {
    # Input validation
    if (!is.data.frame(PathwayVsMetabolites) || 
        !all(c("Pathway", "Metabolites") %in% colnames(PathwayVsMetabolites))) {
        stop("PathwayVsMetabolites must be a data frame with 'Pathway' and 'Metabolites' columns")
    }
    
    # Clean the data first - remove NA and empty values
    PathwayVsMetabolites_clean <- PathwayVsMetabolites %>%
        dplyr::filter(!is.na(Pathway), 
                      !is.na(Metabolites),
                      Pathway != "",
                      Metabolites != "",
                      !grepl("^\\s*$", Metabolites))  # Remove whitespace-only entries
    
    if (nrow(PathwayVsMetabolites_clean) == 0) {
        warning("No valid pathway-metabolite relationships found after cleaning")
        return(data.frame(
            Metabolite = character(0),
            RBC_Metabolite = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    
    # Prepare edge list with additional cleaning
    data <- PathwayVsMetabolites_clean %>%
        dplyr::mutate(
            Metabolites = strsplit(as.character(Metabolites), ",")
        ) %>%
        tidyr::unnest(Metabolites) %>%
        dplyr::mutate(
            Metabolites = trimws(Metabolites)  # Trim whitespace
        ) %>%
        dplyr::filter(
            Metabolites != "",  # Remove empty metabolites
            !is.na(Metabolites)  # Remove NA metabolites
        )
    
    if (nrow(data) == 0) {
        warning("No valid metabolite entries found after cleaning")
        return(data.frame(
            Metabolite = character(0),
            RBC_Metabolite = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    
    # Create edge list
    edge_list_metabolites <- data.frame(
        from = data$Metabolites,
        to = data$Pathway,
        stringsAsFactors = FALSE
    )
    
    # Remove any remaining NA values in edge list
    edge_list_metabolites <- edge_list_metabolites %>%
        dplyr::filter(!is.na(from), !is.na(to), from != "", to != "")
    
    if (nrow(edge_list_metabolites) == 0) {
        warning("No valid edges found for graph construction")
        return(data.frame(
            Metabolite = character(0),
            RBC_Metabolite = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    
    # Create graph and calculate centrality
    tryCatch({
        g_metabolites <- igraph::graph_from_data_frame(
            d = edge_list_metabolites, 
            directed = FALSE
        )
        
        betweenness_metabolites <- igraph::betweenness(
            g_metabolites, 
            directed = FALSE, 
            normalized = TRUE
        )
        
        # Prepare results
        metabolite_centrality <- data.frame(
            Metabolite = names(betweenness_metabolites),
            RBC_Metabolite = betweenness_metabolites,
            stringsAsFactors = FALSE
        ) %>%
            dplyr::arrange(desc(RBC_Metabolite))
        
        # Filter to keep only metabolites (remove pathways if they appear)
        # Metabolites typically start with "C" followed by numbers
        metabolite_centrality <- metabolite_centrality %>%
            dplyr::filter(grepl("^C\\d+", Metabolite))
        
        return(metabolite_centrality)
        
    }, error = function(e) {
        warning("Graph construction or centrality calculation failed: ", e$message)
        # Return empty data frame
        return(data.frame(
            Metabolite = character(0),
            RBC_Metabolite = numeric(0),
            stringsAsFactors = FALSE
        ))
    })
}