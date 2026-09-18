#' Create Metabolite Pathway Network
#'
#' Generates a metabolite pathway interaction network visualization.
#'
#' @param inputMetabolites Character vector of metabolite IDs
#' @param PathwayVsMetabolites Data frame with pathway metabolite associations
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping
#' @param top_n Number of top metabolites to display based on centrality
#' @param font_family Font family for text elements
#'
#' @return A ggraph object showing the metabolite pathway network
#'
#' @examples
#' 
#' # Always-runnable minimal example (no network)
#' pw <- data.frame(
#'     Pathway = c("Glycolysis", "TCA cycle"),
#'     Metabolites = c("C00031,C00022,C00074", "C00036,C00042,C00026"),
#'     stringsAsFactors = FALSE
#' )
#' mets <- c("C00031", "C00022", "C00036")
#' create_network_plot(mets, pw, kegg_lookup = NULL, top_n = 2)
#'
#'
#' PathwayVsMetabolites <- fetch_kegg_pathway_metabolites(organism = "hsa")
#' kegg_lookup <- fetch_kegg_compound_lookup()
#'
#' # Metabolomics matrix (bundled example data)
#' metabolomics_path <- system.file(
#'     "extdata", "example_data.csv",
#'     package = "enrichmet"
#' )
#'
#' if (metabolomics_path == "") {
#'     stop("Example file 'example_data.csv' not found in inst/extdata/")
#' }
#'
#' metabolomics_mat <- read.csv(
#'     metabolomics_path,
#'     row.names = 1,
#'     check.names = FALSE
#' )
#' da_out <- run_de(
#'     metabolomics_mat, "TK-CMV", "K-CMV",
#'     fc_threshold = 1, pval_threshold = 0.05
#' )
#' query <- da_out$kegg_ready$kegg_id[
#'     da_out$kegg_ready$Significant != "Not significant"
#' ]
#' query <- query[!is.na(query) & query != ""]
#'
#' create_network_plot(
#'     query,
#'     PathwayVsMetabolites,
#'     kegg_lookup = kegg_lookup,
#'     top_n = 10
#' )
#' 
#'
#' @export
create_network_plot <- function(
        inputMetabolites,
        PathwayVsMetabolites,
        kegg_lookup = NULL,
        top_n = 20,
        font_family = "sans"
) {
    
    if (!is.null(top_n) && top_n < length(inputMetabolites)) {
        
        all_centrality <- calculate_metabolite_centrality(PathwayVsMetabolites)
        
        selected_metabolites <- all_centrality |>
            dplyr::filter(Metabolite %in% inputMetabolites) |>
            dplyr::arrange(desc(RBC_Metabolite)) |>
            head(top_n) |>
            dplyr::pull(Metabolite)
        
        message("Using top ", top_n, " metabolites by centrality")
    } else {
        selected_metabolites <- inputMetabolites
    }
    
    df <- PathwayVsMetabolites |>
        dplyr::mutate(Metabolite = strsplit(Metabolites, ",")) |>
        tidyr::unnest(Metabolite) |>
        dplyr::filter(Metabolite %in% selected_metabolites)
    
    if (nrow(df) == 0) {
        warning("No matching metabolites found for network plot")
        return(NULL)
    }
    
    edges <- df |> dplyr::select(Pathway, Metabolite)
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    
    igraph::V(g)$type <- ifelse(
        igraph::V(g)$name %in% df$Pathway,
        "Pathway",
        "Metabolite"
    )
    
    igraph::V(g)$degree <- igraph::degree(g)
    
    if (!is.null(kegg_lookup)) {
        metabolite_names <- kegg_lookup$name
        names(metabolite_names) <- kegg_lookup$kegg_id
        
        igraph::V(g)$display_name <- ifelse(
            igraph::V(g)$type == "Metabolite",
            ifelse(
                !is.na(metabolite_names[igraph::V(g)$name]),
                metabolite_names[igraph::V(g)$name],
                igraph::V(g)$name
            ),
            igraph::V(g)$name
        )
    } else {
        igraph::V(g)$display_name <- igraph::V(g)$name
    }
    
    network_plot <- ggraph::ggraph(g, layout = "fr") +
        ggraph::geom_edge_link(
            alpha = 0.3,
            color = "#808080"
        ) +
        ggraph::geom_node_point(
            ggplot2::aes(
                color = type,
                size = degree
            ),
            alpha = 0.85
        ) +
        ggraph::geom_node_text(
            ggplot2::aes(label = display_name),
            repel = TRUE,
            size = 3,
            family = font_family,
            fontface = ifelse(
                igraph::V(g)$type == "Pathway",
                "bold",
                "plain"
            ),
            max.overlaps = 50
        ) +
        ggplot2::scale_color_manual(
            name = "Node type",
            values = c(
                Metabolite = "#1F78B4",
                Pathway = "#33A02C"
            )
        ) +
        ggplot2::scale_size_continuous(
            name = "Number of connections",
            range = c(2, 8)
        ) +
        ggplot2::labs(
            title = paste(
                "Metabolite Pathway Network Top",
                top_n,
                "metabolites"
            ),
            subtitle = paste(
                igraph::vcount(g),
                "nodes and",
                igraph::ecount(g),
                "edges"
            )
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold",
                hjust = 0.5,
                size = 14,
                family = font_family
            ),
            plot.subtitle = ggplot2::element_text(
                hjust = 0.5,
                size = 11,
                family = font_family
            ),
            legend.position = "right",
            legend.title = ggplot2::element_text(
                face = "bold",
                family = font_family
            ),
            legend.text = ggplot2::element_text(
                family = font_family
            )
        )
    
    return(network_plot)
}
