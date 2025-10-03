#' Calculate Metabolite Centrality
#'
#' Computes relative betweenness centrality for metabolites in pathway networks.
#'
#' @param PathwayVsMetabolites A data frame containing pathways and their associated 
#'        metabolites. Must include columns 'Pathway' and 'Metabolites'.
#'
#' @return A data frame with metabolite names and their relative betweenness centrality scores.
#'
#' @examples
#' # Create example pathway data
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Pathway1", "Pathway2"),
#'   Metabolites = c("M1,M2,M3", "M2,M4,M5")
#' )
#'
#' # Calculate centrality
#' centrality <- calculate_metabolite_centrality(PathwayVsMetabolites)
#' head(centrality)
#' # See ?enrichmet for complete examples
#' @importFrom dplyr filter
#' @importFrom tidyr unnest
#' @importFrom igraph graph_from_data_frame betweenness
#' @export
calculate_metabolite_centrality <- function(PathwayVsMetabolites) {
    # Input validation
    if (!is.data.frame(PathwayVsMetabolites) || 
        !all(c("Pathway", "Metabolites") %in% colnames(PathwayVsMetabolites))) {
        stop("PathwayVsMetabolites must be a data frame with 'Pathway' and 'Metabolites' columns")
    }
    
    # Prepare edge list - EXACTLY LIKE ORIGINAL
    data <- PathwayVsMetabolites %>%
        dplyr::mutate(Metabolites = strsplit(Metabolites, ",")) %>%
        tidyr::unnest(Metabolites)
    
    edge_list_metabolites <- data.frame(
        from = unlist(data$Metabolites),
        to = rep(data$Pathway, lengths(data$Metabolites))
    )
    
    # Create graph and calculate centrality - EXACTLY LIKE ORIGINAL
    g_metabolites <- igraph::graph_from_data_frame(d = edge_list_metabolites, directed = FALSE)
    betweenness_metabolites <- igraph::betweenness(g_metabolites, directed = FALSE, normalized = TRUE)
    
    # Prepare results - EXACTLY LIKE ORIGINAL (NO FILTERING)
    metabolite_centrality <- data.frame(
        Metabolite = names(betweenness_metabolites),
        RBC_Metabolite = betweenness_metabolites
    )
    
    return(metabolite_centrality)
}
