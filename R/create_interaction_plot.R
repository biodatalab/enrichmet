#' Create Metabolite Interaction Network (Reactome)
#'
#' Generates a metabolite interaction network based on shared Reactome
#' reactions. Two metabolites are connected by an edge if they co-occur in
#' the same Reactome reaction; edge weight reflects the number of reactions
#' they share.
#'
#' @param inputMetabolites A character vector of metabolite IDs (KEGG IDs).
#' @param reactome_df A data frame with the KEGG-to-Reactome reaction
#'        mapping, e.g. the result of joining a KEGG-ChEBI conversion table
#'        to Reactome's ChEBI2ReactomeReactions.txt. Must contain columns
#'        'KEGG' and 'Reaction'. If you want to restrict to a species (e.g.
#'        "Homo sapiens"), filter reactome_df on 'Species' before calling
#'        this function.
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping.
#'        Should contain columns 'kegg_id' and 'name'.
#' @param min_shared_reactions Minimum number of shared Reactome reactions
#'        required for two metabolites to be connected by an edge
#'        (default 1).
#' @param max_metabolites_per_reaction Reactions involving more than this
#'        many of the input metabolites are skipped when building edges.
#'        This avoids hairball edges coming from very generic/large
#'        reactions that touch many unrelated metabolites (default 25).
#'
#' @return A ggraph object showing the Reactome-based metabolite
#'         interaction network, or NULL if insufficient data is available.
#'
#' @examples
#' reactome_df <- data.frame(
#'   KEGG = c("C00031", "C00031", "C00022", "C00022", "C00074", "C00074"),
#'   ChEBI = c("4167", "4167", "16651", "16651", "16452", "16452"),
#'   Reaction = c("R-HSA-1", "R-HSA-2", "R-HSA-1", "R-HSA-3",
#'                "R-HSA-2", "R-HSA-3"),
#'   Species = "Homo sapiens"
#' )
#' plot <- create_interaction_plot(
#'   inputMetabolites = c("C00031", "C00022", "C00074"),
#'   reactome_df = reactome_df
#' )
#' plot
#'
#' @importFrom dplyr filter distinct mutate select group_by summarise
#' @importFrom dplyr n_distinct left_join everything
#' @importFrom igraph graph_from_data_frame degree betweenness components
#' @importFrom igraph vcount ecount
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggraph theme_graph scale_edge_width create_layout
#' @importFrom scales rescale
#' @importFrom stringr str_trunc
#' @importFrom utils combn
#' @export
create_interaction_plot <- function(inputMetabolites,
                                    reactome_df,
                                    kegg_lookup = NULL,
                                    min_shared_reactions = 1,
                                    max_metabolites_per_reaction = 25) {
    
    if (is.null(reactome_df) || nrow(reactome_df) == 0) {
        warning("Reactome interaction analysis requested but reactome_df not provided or empty.")
        return(NULL)
    }
    
    if (!all(c("KEGG", "Reaction") %in% colnames(reactome_df))) {
        warning("reactome_df must contain 'KEGG' and 'Reaction' columns.")
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
    valid_kegg_ids <- unique(extracted_kegg_ids[!is.na(extracted_kegg_ids) &
                                                    extracted_kegg_ids != ""])
    
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
    
    # Restrict the Reactome mapping to metabolites of interest
    reactome_sub <- reactome_df %>%
        dplyr::filter(KEGG %in% valid_kegg_ids) %>%
        dplyr::distinct(KEGG, Reaction, .keep_all = TRUE)
    
    message("Found ", dplyr::n_distinct(reactome_sub$KEGG),
            " metabolites with ", dplyr::n_distinct(reactome_sub$Reaction),
            " associated Reactome reactions")
    
    if (nrow(reactome_sub) == 0) {
        warning("No Reactome reactions found for the supplied metabolites")
        return(NULL)
    }
    
    # Build edges: two metabolites are connected if they co-occur in the
    # same Reactome reaction. Weight = number of shared reactions.
    message("Building metabolite co-occurrence edges from shared reactions...")
    reactions_split <- split(reactome_sub$KEGG, reactome_sub$Reaction)
    
    edge_list <- lapply(reactions_split, function(mets_in_reaction) {
        mets_in_reaction <- unique(mets_in_reaction)
        n_mets <- length(mets_in_reaction)
        if (n_mets < 2 || n_mets > max_metabolites_per_reaction) return(NULL)
        
        pairs <- utils::combn(sort(mets_in_reaction), 2, simplify = FALSE)
        do.call(rbind, lapply(pairs, function(p) {
            data.frame(from = p[1], to = p[2], stringsAsFactors = FALSE)
        }))
    })
    edge_list <- edge_list[!vapply(edge_list, is.null, logical(1))]
    
    if (length(edge_list) == 0) {
        warning("No shared reactions found between any pair of input metabolites")
        return(NULL)
    }
    
    all_edges <- do.call(rbind, edge_list)
    
    valid_edges <- all_edges %>%
        dplyr::group_by(from, to) %>%
        dplyr::summarise(shared_reactions = dplyr::n(), .groups = "drop") %>%
        dplyr::filter(shared_reactions >= min_shared_reactions)
    
    message("Found ", nrow(valid_edges), " valid metabolite pairs sharing >= ",
            min_shared_reactions, " reaction(s)")
    
    if (nrow(valid_edges) == 0) {
        warning("Insufficient Reactome data to build interaction graph after filtering.")
        return(NULL)
    }
    
    # Create vertex_df using metabolites present in the valid edges
    vertex_ids <- unique(c(valid_edges$from, valid_edges$to))
    vertex_df <- data.frame(KEGG_ID = vertex_ids, stringsAsFactors = FALSE)
    
    # KEGG name mapping (optional)
    if (!is.null(kegg_lookup)) {
        if (all(c("kegg_id", "name") %in% colnames(kegg_lookup))) {
            vertex_df <- vertex_df %>%
                dplyr::left_join(kegg_lookup,
                                 by = c("KEGG_ID" = "kegg_id")) %>%
                dplyr::mutate(display_name = ifelse(!is.na(name),
                                                    name, KEGG_ID)) %>%
                dplyr::select(-name)
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
        dplyr::select(name = KEGG_ID, display_name, dplyr::everything())
    
    # Debug: Show what display names we have
    if (nrow(vertex_df) > 0) {
        message("Display names sample: ",
                paste(utils::head(vertex_df$display_name), collapse = ", "))
    }
    
    message("Creating graph with ", nrow(vertex_df), " vertices")
    message("Vertex attributes: ", paste(names(vertex_df), collapse = ", "))
    
    # Create graph
    g <- igraph::graph_from_data_frame(
        d = valid_edges %>%
            dplyr::mutate(weight = scales::rescale(shared_reactions,
                                                   to = c(0.1, 1))),
        directed = FALSE,
        vertices = vertex_df
    )
    
    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
        warning("Reactome graph has no vertices or edges.")
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
            aes(label = display_name),
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
            name = "Shared Reactions",
            range = c(0.5, 2)
        ) +
        ggplot2::labs(
            title = "Metabolite Interaction Network (Reactome)",
            subtitle = paste(
                igraph::vcount(g),
                "compounds with",
                igraph::ecount(g),
                "shared-reaction edges"
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
