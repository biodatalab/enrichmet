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
create_interaction_plot <- function(inputMetabolites, mapping_df, stitch_df, 
                                    kegg_lookup = NULL) {
    if (is.null(mapping_df) || is.null(stitch_df)) {
        warning("STITCH interaction analysis requested but mapping_df or stitch_df not provided.")
        return(NULL)
    }
    
    # Helper function to extract KEGG ID from complex metabolite IDs
    extract_kegg_id <- function(met_id) {
        if (is.na(met_id) || met_id == "") return(NA_character_)
        
        if (grepl("^C\\d{5}$", met_id)) {
            return(met_id)
        }
        
        if (grepl("_C\\d", met_id)) {
            parts <- unlist(strsplit(met_id, "_"))
            kegg_parts <- grep("^C\\d", parts, value = TRUE)
            
            if (length(kegg_parts) > 0) {
                if (grepl("\\|", kegg_parts[1])) {
                    all_kegg <- unlist(strsplit(kegg_parts[1], "\\|"))
                    valid_kegg <- grep("^C\\d", all_kegg, value = TRUE)
                    if (length(valid_kegg) > 0) return(valid_kegg[1])
                } else {
                    return(kegg_parts[1])
                }
            }
        }
        
        return(NA_character_)
    }
    
    # Extract KEGG IDs from input metabolites
    message("Extracting KEGG IDs from input metabolites...")
    extracted_kegg_ids <- vapply(inputMetabolites, extract_kegg_id, 
                                 FUN.VALUE = character(1))
    valid_kegg_ids <- extracted_kegg_ids[!is.na(extracted_kegg_ids) & 
                                             extracted_kegg_ids != ""]
    
    message("Successfully extracted ", length(valid_kegg_ids), 
            " KEGG IDs from ", length(inputMetabolites), " input metabolites")
    
    if (length(valid_kegg_ids) > 0) {
        message("Sample extracted KEGG IDs: ", 
                paste(utils::head(unique(valid_kegg_ids)), collapse = ", "))
    }
    
    if (length(valid_kegg_ids) == 0) {
        warning("No valid KEGG IDs could be extracted from input metabolites")
        return(NULL)
    }
    
    # Create vertex_df using extracted KEGG IDs
    vertex_df <- mapping_df %>%
        dplyr::filter(KEGG_ID %in% valid_kegg_ids) %>%
        dplyr::filter(!is.na(PubChem_CID)) %>%
        dplyr::distinct(KEGG_ID, .keep_all = TRUE)
    
    message("Found ", nrow(vertex_df), 
            " metabolites in mapping_df with valid PubChem CIDs")
    
    if (nrow(vertex_df) == 0) {
        warning("No metabolites found in mapping_df after KEGG ID extraction and filtering")
        return(NULL)
    }
    
    # KEGG name mapping (optional)
    if (!is.null(kegg_lookup)) {
        if (all(c("kegg_id", "name") %in% colnames(kegg_lookup))) {
            vertex_df <- vertex_df %>%
                dplyr::left_join(kegg_lookup, 
                                 by = c("KEGG_ID" = "kegg_id")) %>%
                dplyr::mutate(display_name = ifelse(!is.na(name), 
                                                    name, KEGG_ID))
            message("Applied KEGG pathway name mapping")
        } else {
            warning("kegg_lookup provided but missing required columns 'kegg_id' and 'name'")
            vertex_df <- vertex_df %>%
                dplyr::mutate(display_name = KEGG_ID)
        }
    } else {
        vertex_df <- vertex_df %>%
            dplyr::mutate(display_name = KEGG_ID)
    }
    
    vertex_df <- vertex_df %>%
        dplyr::mutate(display_name = stringr::str_trunc(display_name, 25)) %>%
        dplyr::select(STITCH_ID, dplyr::everything())
    
    # Debug: Show what display names we have
    if (nrow(vertex_df) > 0) {
        message("Display names sample: ", 
                paste(utils::head(vertex_df$display_name), collapse = ", "))
    }
    
    # Create edge list
    valid_edges <- stitch_df %>%
        dplyr::filter(combined_score >= 50) %>%
        dplyr::semi_join(vertex_df, by = c("chemical1" = "STITCH_ID")) %>%
        dplyr::semi_join(vertex_df, by = c("chemical2" = "STITCH_ID")) %>%
        dplyr::distinct(chemical1, chemical2, .keep_all = TRUE)
    
    message("Found ", nrow(valid_edges), " valid interactions between ", 
            nrow(vertex_df), " metabolites")
    
    if (nrow(vertex_df) == 0 || nrow(valid_edges) == 0) {
        warning("Insufficient STITCH data to build interaction graph.")
        return(NULL)
    }
    
    # Create graph - CRITICAL: Make sure vertex_df has the right columns
    # The vertices data frame for igraph::graph_from_data_frame should have 
    # name as first column
    vertices_for_graph <- vertex_df %>%
        dplyr::filter(STITCH_ID %in% c(valid_edges$chemical1, 
                                       valid_edges$chemical2)) %>%
        dplyr::select(name = STITCH_ID, display_name, KEGG_ID, PubChem_CID)  
    
    message("Creating graph with ", nrow(vertices_for_graph), " vertices")
    message("Vertex attributes: ", paste(names(vertices_for_graph), 
                                         collapse = ", "))
    
    g <- igraph::graph_from_data_frame(
        d = valid_edges %>% 
            dplyr::select(chemical1, chemical2, dplyr::everything()) %>%
            dplyr::mutate(weight = scales::rescale(combined_score, 
                                                   to = c(0.1, 1))),
        directed = FALSE,
        vertices = vertices_for_graph  # Use the properly formatted vertices
    )
    
    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
        warning("STITCH graph has no vertices or edges.")
        return(NULL)
    }
    
    # Debug: Check what attributes are in the graph
    message("Graph vertex attributes: ", 
            paste(names(igraph::vertex_attr(g)), collapse = ", "))
    message("Sample vertex display_names: ", 
            paste(utils::head(igraph::V(g)$display_name), collapse = ", "))
    
    # Calculate graph metrics
    igraph::V(g)$degree <- igraph::degree(g)
    igraph::V(g)$betweenness <- igraph::betweenness(g)
    igraph::V(g)$component <- igraph::components(g)$membership
    
    # Get actual degree range for meaningful breaks
    degree_range <- range(igraph::V(g)$degree)
    degree_breaks <- seq(degree_range[1], degree_range[2], by = 1)
    
    if (length(degree_breaks) > 5) {
        degree_breaks <- pretty(degree_range, n = 4)
    }
    
    # Check if graph is disconnected
    comps <- igraph::components(g)
    num_components <- comps$no
    
    if (num_components > 1) {
        message("Graph has ", num_components, 
                " disconnected components - using GEM layout")
        best_layout <- "gem"
    } else {
        # Try multiple layouts for connected graph
        layouts_to_try <- c("fr", "kk", "dh", "gem", "lgl")
        best_layout <- NULL
        best_spacing <- 0
        
        for (layout_name in layouts_to_try) {
            tryCatch({
                layout_pos <- ggraph::create_layout(g, layout = layout_name)
                node_distances <- as.matrix(stats::dist(layout_pos[, 
                                                                   seq_len(2)]))
                diag(node_distances) <- NA
                avg_distance <- mean(node_distances, na.rm = TRUE)
                
                if (avg_distance > best_spacing) {
                    best_spacing <- avg_distance
                    best_layout <- layout_name
                }
            }, error = function(e) NULL)
        }
        
        if (is.null(best_layout)) best_layout <- "fr"
        message("Using layout: ", best_layout, " (spacing score: ", 
                round(best_spacing, 2), ")")
    }
    
    # Create the plot
    interaction_plot <- ggraph::ggraph(g, layout = best_layout) +
        # Draw edges
        ggraph::geom_edge_link(
            aes(width = weight), 
            color = "#606060",
            alpha = 0.4,
            show.legend = TRUE
        ) +
        # Draw nodes
        ggraph::geom_node_point(
            aes(size = degree, color = as.factor(component)), 
            alpha = 0.8,
            stroke = 0.5
        ) +
        # Add labels - use display_name from vertex attributes
        ggraph::geom_node_text(
            aes(label = display_name),  # Just use display_name directly
            size = 3.5,
            repel = TRUE,
            box.padding = 0.8,
            point.padding = 0.5,
            max.overlaps = Inf,
            min.segment.length = 0.2,
            segment.color = "grey30",
            segment.alpha = 0.6,
            segment.size = 0.3,
            family = "sans",
            fontface = "bold"
        ) +
        ggplot2::scale_color_discrete(name = "Network Component") +
        ggplot2::scale_size_continuous(
            name = "Degree (Connections)", 
            range = c(3, 10),
            breaks = degree_breaks,
            guide = ggplot2::guide_legend(
                override.aes = list(color = "grey50"),
                nrow = min(4, length(degree_breaks))
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
                "interactions | Layout:",
                toupper(best_layout)
            )
        ) +
        ggraph::theme_graph(base_family = "sans") +
        ggplot2::theme(
            legend.position = "right",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", 
                                               family = "sans"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, 
                                                  family = "sans"),
            legend.key.height = ggplot2::unit(0.8, "lines")
        )
    
    return(interaction_plot)
}