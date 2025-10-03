#' Create Metabolite-Pathway Network
#'
#' Generates a metabolite-pathway interaction network visualization.
#'
#' @param inputMetabolites A character vector of metabolite IDs.
#' @param PathwayVsMetabolites A data frame with pathway-metabolite associations.
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping.
#'
#' @return A ggraph object showing the metabolite-pathway network.
#'
#' @examples
#' # Create comprehensive example data
#' inputMetabolites <- paste0("C", sprintf("%05d", 1:25))
#' 
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Glycolysis", "TCA cycle", "Pentose phosphate", 
#'               "Amino acid metabolism", "Lipid metabolism",
#'               "Nucleotide metabolism", "Oxidative phosphorylation"),
#'   Metabolites = c("C00031,C00022,C00118,C00197,C00074,C00036,C00186,C00221,C00089,C00103",
#'                   "C00024,C00036,C00042,C00122,C00149,C00158,C00417,C05379",
#'                   "C00117,C00231,C00279,C00345,C01172,C01236",
#'                   "C00025,C00026,C00041,C00049,C00064,C00135",
#'                   "C00083,C00162,C00249,C00422,C01205",
#'                   "C00106,C00112,C00144,C00212,C00262",
#'                   "C00009,C00013,C00080,C01345,C05345")
#' )
#' 
#' # Create KEGG lookup
#' kegg_lookup <- data.frame(
#'   kegg_id = paste0("C", sprintf("%05d", 1:30)),
#'   name = c("Glucose", "Lactate", "Pyruvate", "Alanine", "Glutamate",
#'            "Citrate", "Succinate", "Malate", "Aspartate", "Glutamine",
#'            "Acetyl-CoA", "Oxaloacetate", "Alpha-ketoglutarate", "Fumarate",
#'            "Isocitrate", "Ribose-5P", "Glycerol-3P", "Dihydroxyacetone-P",
#'            "Fructose-6P", "Glucose-6P", "Phosphoenolpyruvate", "ATP", "ADP",
#'            "NAD+", "NADH", "Coenzyme A", "Acetate", "Butyrate", "Propionate", "Glycine")
#' )
#'
#' # Create network plot
#' plot <- create_network_plot(inputMetabolites, PathwayVsMetabolites, kegg_lookup)
#' plot
#'
#' @importFrom dplyr filter select mutate
#' @importFrom tidyr unnest
#' @importFrom igraph graph_from_data_frame betweenness degree vcount ecount
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text theme_graph
#' @export
create_network_plot <- function(inputMetabolites, PathwayVsMetabolites, kegg_lookup = NULL) {
    df <- PathwayVsMetabolites %>%
        dplyr::mutate(Metabolite = strsplit(Metabolites, ",")) %>%
        tidyr::unnest(Metabolite) %>%
        dplyr::filter(Metabolite %in% inputMetabolites)
    
    if (nrow(df) == 0) {
        warning("No matching metabolites found for network plot")
        return(NULL)
    }
    
    edges <- df %>% dplyr::select(Pathway, Metabolite)
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    
    # Calculate centrality for node sizing
    node_centrality <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)
    
    # Set node attributes
    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% df$Pathway, "Pathway", "Metabolite")
    igraph::V(g)$centrality <- node_centrality[igraph::V(g)$name]
    igraph::V(g)$degree <- igraph::degree(g)
    
    # Add metabolite names if available
    if (!is.null(kegg_lookup)) {
        metabolite_names <- kegg_lookup$name
        names(metabolite_names) <- kegg_lookup$kegg_id
        
        igraph::V(g)$display_name <- ifelse(
            igraph::V(g)$type == "Metabolite",
            ifelse(!is.na(metabolite_names[igraph::V(g)$name]), 
                   metabolite_names[igraph::V(g)$name], 
                   igraph::V(g)$name),
            igraph::V(g)$name
        )
    } else {
        igraph::V(g)$display_name <- igraph::V(g)$name
    }
    
    # Create network plot
    network_plot <- ggraph::ggraph(g, layout = "fr") +
        ggraph::geom_edge_link(alpha = 0.3, color = "grey70") +
        ggraph::geom_node_point(aes(color = type, size = degree), alpha = 0.8) +
        ggraph::geom_node_text(
            aes(label = display_name),
            size = 3, repel = TRUE, max.overlaps = 50,
            fontface = ifelse(igraph::V(g)$type == "Pathway", "bold", "plain")
        ) +
        ggplot2::scale_color_manual(
            name = "Node Type",
            values = c("Metabolite" = "red", "Pathway" = "blue"),
            labels = c("Metabolite", "Pathway")
        ) +
        ggplot2::scale_size_continuous(name = "Connections", range = c(2, 8)) +
        ggplot2::labs(
            title = "Metabolite-Pathway Network",
            subtitle = paste("Showing", igraph::vcount(g), "nodes and", igraph::ecount(g), "connections")
        ) +
        ggraph::theme_graph() +
        ggplot2::theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)
        )
    
    return(network_plot)
}