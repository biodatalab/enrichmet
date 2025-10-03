#' Create STITCH Interaction Network
#'
#' Generates a chemical interaction network using STITCH data.
#'
#' @param inputMetabolites A character vector of metabolite IDs.
#' @param mapping_df Data frame mapping KEGG IDs to STITCH IDs and PubChem CIDs.
#' @param stitch_df Data frame containing STITCH interaction data.
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping.
#'
#' @return A ggraph object showing the chemical interaction network.
#'
#' @examples
#' # Create comprehensive example STITCH data
#' set.seed(123)
#' inputMetabolites <- c("C00031", "C00022", "C00074", "C00036", "C00103", 
#'                      "C00197", "C00186", "C00221", "C00024", "C00042")
#'
#' mapping_df <- data.frame(
#'   KEGG_ID = c("C00031", "C00022", "C00074", "C00036", "C00103", 
#'               "C00197", "C00186", "C00221", "C00024", "C00042",
#'               "C00122", "C00149", "C00158", "C00117"),
#'   STITCH_ID = c("CID00000579", "CID00000607", "CID00000586", "CID00000387",
#'                 "CID00000534", "CID00000734", "CID00000596", "CID00000551",
#'                 "CID00000176", "CID00000946", "CID00000330", "CID00001097",
#'                 "CID00000625", "CID00000867"),
#'   PubChem_CID = c(5793, 607, 586, 387, 534, 734, 596, 551, 176, 946, 
#'                   330, 1097, 625, 867)
#' )
#'
#' # Create STITCH interaction data
#' stitch_df <- data.frame(
#'   chemical1 = c("CID00000579", "CID00000579", "CID00000607", "CID00000586",
#'                 "CID00000387", "CID00000534", "CID00000734", "CID00000596",
#'                 "CID00000551", "CID00000176", "CID00000946", "CID00000330",
#'                 "CID00001097", "CID00000625", "CID00000867"),
#'   chemical2 = c("CID00000534", "CID00000734", "CID00000586", "CID00000387",
#'                 "CID00000551", "CID00000946", "CID00000596", "CID00000176",
#'                 "CID00000330", "CID00001097", "CID00000625", "CID00000867",
#'                 "CID00000579", "CID00000607", "CID00000534"),
#'   combined_score = c(850, 720, 680, 790, 810, 730, 690, 760, 820, 710,
#'                      750, 780, 670, 740, 710)
#' )
#'
#' kegg_lookup <- data.frame(
#'   kegg_id = c("C00031", "C00022", "C00074", "C00036", "C00103", 
#'               "C00197", "C00186", "C00221", "C00024", "C00042",
#'               "C00122", "C00149", "C00158", "C00117"),
#'   name = c("D-Glucose", "L-Lactate", "Oxaloacetate", "Citrate", 
#'            "D-Fructose 6-phosphate", "3-Phospho-D-glyceroyl phosphate",
#'            "L-Aspartate", "Oxoglutaric acid", "Acetyl-CoA", 
#'            "Oxaloacetic acid", "L-Glutamate", "Succinyl-CoA",
#'            "L-Malate", "D-Ribose 5-phosphate")
#' )
#'
#' # Create interaction network plot
#' plot <- create_interaction_plot(inputMetabolites, mapping_df, stitch_df, kegg_lookup)
#' plot
#'
#' @importFrom dplyr filter distinct mutate select semi_join
#' @importFrom igraph graph_from_data_frame degree betweenness components vcount ecount
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text theme_graph scale_edge_width
#' @importFrom scales rescale
#' @importFrom stringr str_wrap str_trunc
#' @export
create_interaction_plot <- function(inputMetabolites, mapping_df, stitch_df, kegg_lookup = NULL) {
    if (is.null(mapping_df) || is.null(stitch_df)) {
        warning("STITCH interaction analysis requested but mapping_df or stitch_df not provided.")
        return(NULL)
    }
    
    # Create vertex_df
    vertex_df <- mapping_df %>%
        dplyr::filter(KEGG_ID %in% inputMetabolites) %>%
        dplyr::filter(!is.na(PubChem_CID)) %>%
        dplyr::distinct(KEGG_ID, .keep_all = TRUE)
    
    # KEGG name mapping (optional)
    if (!is.null(kegg_lookup)) {
        vertex_df <- vertex_df %>%
            dplyr::left_join(kegg_lookup, by = c("KEGG_ID" = "kegg_id")) %>%
            dplyr::mutate(display_name = ifelse(!is.na(name), name, KEGG_ID))
    } else {
        vertex_df <- vertex_df %>%
            dplyr::mutate(display_name = KEGG_ID)
    }
    
    vertex_df <- vertex_df %>%
        dplyr::mutate(display_name = stringr::str_trunc(display_name, 25)) %>%
        dplyr::select(STITCH_ID, everything())
    
    # Create edge list
    valid_edges <- stitch_df %>%
        dplyr::filter(combined_score >= 100) %>%
        dplyr::semi_join(vertex_df, by = c("chemical1" = "STITCH_ID")) %>%
        dplyr::semi_join(vertex_df, by = c("chemical2" = "STITCH_ID")) %>%
        dplyr::distinct(chemical1, chemical2, .keep_all = TRUE)
    
    # Proceed only if there are vertices and edges
    if (nrow(vertex_df) == 0 || nrow(valid_edges) == 0) {
        warning("Insufficient STITCH data to build interaction graph.")
        return(NULL)
    }
    
    g <- igraph::graph_from_data_frame(
        d = valid_edges %>% dplyr::mutate(weight = scales::rescale(combined_score, to = c(0.1, 1))),
        directed = FALSE,
        vertices = vertex_df %>%
            dplyr::filter(STITCH_ID %in% c(valid_edges$chemical1, valid_edges$chemical2))
    )
    
    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
        warning("STITCH graph has no vertices or edges.")
        return(NULL)
    }
    
    igraph::V(g)$degree <- igraph::degree(g)
    igraph::V(g)$betweenness <- igraph::betweenness(g)
    igraph::V(g)$component <- igraph::components(g)$membership
    
    # Get actual degree range for meaningful breaks
    degree_range <- range(igraph::V(g)$degree)
    degree_breaks <- seq(degree_range[1], degree_range[2], by = 1)
    
    # If there are too many breaks, use a subset
    if (length(degree_breaks) > 5) {
        degree_breaks <- pretty(degree_range, n = 4)
    }
    
    interaction_plot <- ggraph::ggraph(g, layout = "fr") +
        ggraph::geom_edge_link(aes(width = weight), color = "grey50", alpha = 0.7, show.legend = TRUE) +
        ggraph::geom_node_point(aes(size = degree, color = as.factor(component)), alpha = 0.8) +
        ggraph::geom_node_text(
            aes(label = stringr::str_wrap(display_name, 12)),
            size = 3,
            repel = TRUE,
            max.overlaps = 20,
            family = "sans"
        ) +
        ggplot2::scale_color_discrete(name = "Network Component") +
        ggplot2::scale_size_continuous(
            name = "Degree (Connections)", 
            range = c(3, 10),
            breaks = degree_breaks,
            guide = ggplot2::guide_legend(
                override.aes = list(color = "grey50"),
                nrow = length(degree_breaks)
            )
        ) +
        ggraph::scale_edge_width(
            name = "Interaction Strength",
            range = c(0.5, 2)
        ) +
        ggplot2::labs(
            title = "Metabolite Interaction Network (STITCH)",
            subtitle = paste(
                igraph::vcount(g),
                "compounds with",
                igraph::ecount(g),
                "interactions"
            )
        ) +
        ggraph::theme_graph(base_family = "sans") +
        ggplot2::theme(
            legend.position = "right",
            plot.title = element_text(hjust = 0.5, face = "bold", family = "sans"),
            legend.key.height = unit(0.8, "lines")
        )
    
    return(interaction_plot)
}