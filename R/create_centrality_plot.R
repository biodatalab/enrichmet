#' Create Centrality Plot
#'
#' Generates a visualization of metabolite centrality results (relative
#' betweenness centrality, RBC).
#'
#' @param centrality_results Data frame from
#'   \code{calculate_metabolite_centrality()}, typically filtered to the query
#'   metabolites. Must contain \code{Metabolite} and \code{RBC_Metabolite};
#'   may already contain \code{Display_Name}.
#' @param top_n Number of top metabolites to display (default = 20).
#' @param kegg_lookup Optional data frame with columns \code{kegg_id} and
#'   \code{name} for axis labels.
#'
#' @return A \code{ggplot} object (horizontal bar chart of RBC). Returns
#'   \code{NULL} if \code{centrality_results} has zero rows.
#'
#' @details
#' Bars show relative betweenness centrality from the pathway–metabolite
#' annotation graph only. Fold change and \eqn{p}-values are not used in the
#' underlying centrality scores.
#'
#' @importFrom dplyr left_join mutate
#' @importFrom ggplot2 ggplot aes geom_col labs theme_minimal theme
#'   element_text element_rect coord_flip scale_y_continuous
#' @importFrom stats reorder
#' @importFrom utils head
#'
#' @examples
#' # Always-runnable minimal example (no network)
#' cen <- data.frame(
#'     Metabolite = c("C00031", "C00022", "C00074"),
#'     RBC_Metabolite = c(0.12, 0.08, 0.05),
#'     stringsAsFactors = FALSE
#' )
#' create_centrality_plot(cen, top_n = 3)
#'
#' 
#' PathwayVsMetabolites <- fetch_kegg_pathway_metabolites(organism = "hsa")
#' kegg_lookup <- fetch_kegg_compound_lookup()
#'
#' metabolomics_path <- get_cached_file(
#'     "https://zenodo.org/api/records/17819145/files/example_data.csv/content"
#' )
#' metabolomics_mat <- read.csv(
#'     metabolomics_path, row.names = 1, check.names = FALSE
#' )
#' da_out <- run_de(
#'     metabolomics_mat, "TK-CMV", "K-CMV",
#'     fc_threshold = 1, pval_threshold = 0.05
#' )
#'
#' query <- da_out$kegg_ready$kegg_id[
#'     da_out$kegg_ready$Significant != "Not significant"
#' ]
#' query <- query[!is.na(query) & query != ""]
#'
#' cen <- calculate_metabolite_centrality(PathwayVsMetabolites)
#' cen <- cen[cen$Metabolite %in% query, ]
#'
#' create_centrality_plot(cen, top_n = 15, kegg_lookup = kegg_lookup)
#' 
#'
#' @export
create_centrality_plot <- function(centrality_results, top_n = 20, kegg_lookup = NULL) {
    
    if (nrow(centrality_results) == 0) {
        warning("No centrality results to plot")
        return(NULL)
    }
    
    # Prepare top metabolites
    top_central <- head(centrality_results, top_n)
    
    # Add display names if lookup provided
    if (!is.null(kegg_lookup)) {
        top_central <- top_central %>%
            dplyr::left_join(kegg_lookup, by = c("Metabolite" = "kegg_id")) %>%
            dplyr::mutate(
                Display_Name = ifelse(!is.na(name), name, Metabolite)
            )
    } else {
        top_central$Display_Name <- top_central$Metabolite
    }
    
    # Axis limits
    y_min <- 0
    y_max <- max(top_central$RBC_Metabolite, na.rm = TRUE) * 1.05
    
    ggplot2::ggplot(
        top_central,
        ggplot2::aes(
            x = reorder(Display_Name, RBC_Metabolite),
            y = RBC_Metabolite
        )
    ) +
        ggplot2::geom_col(
            alpha = 0.85,
            fill = "#1F78B4"
        ) +
        ggplot2::labs(
            title = paste("Top", top_n, "Metabolites by Relative Betweenness Centrality"),
            subtitle = "Measures connectivity importance in the metabolic network",
            x = "Metabolite",
            y = "Relative Betweenness Centrality"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            axis.text.y = element_text(size = 12, color = "black"),
            axis.text.x = element_text(size = 12, color = "black"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            panel.border = ggplot2::element_rect(
                color = "black",
                fill = NA,
                linewidth = 0.8
            ),
            legend.position = "none"
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(
            expand = c(0, 0),
            limits = c(y_min, y_max)
        )
}
