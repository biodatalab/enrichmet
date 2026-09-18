#' Create Impact vs Significance Plot
#'
#' Generates a visualization of pathway impact versus statistical significance.
#'
#' @param significant_results_df Data frame from \code{perform_enrichment_analysis()}
#'   (or \code{enrichmet()}$pathway_enrichment_results). Must contain columns
#'   \code{Impact}, \code{Log_P_value}, and \code{Pathway}.
#'
#' @return A \code{ggplot} object showing pathway impact (x) versus
#'   \eqn{-\log_{10}} \eqn{p}-value (y). Returns \code{NULL} if
#'   \code{significant_results_df} has zero rows.
#'
#' @details
#' Each point is a pathway. Impact is a topology-based score; the vertical
#' axis is enrichment significance. Pathways that are both statistically
#' supported and structurally central stand out toward the upper right.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal theme
#'   element_text element_rect
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' 
#' PathwayVsMetabolites <- fetch_kegg_pathway_metabolites(organism = "hsa")
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
#' enr <- perform_enrichment_analysis(
#'     inputMetabolites = da_out$kegg_ready$kegg_id[
#'         da_out$kegg_ready$Significant != "Not significant"
#'     ],
#'     PathwayVsMetabolites = PathwayVsMetabolites,
#'     top_n = 15,
#'     p_value_cutoff = 1
#' )
#'
#' create_impact_plot(enr)
#' 
create_impact_plot <- function(significant_results_df) {
    
    if (nrow(significant_results_df) == 0) {
        warning("No pathways passed the p value cutoff.")
        return(NULL)
    }
    
    ggplot2::ggplot(
        significant_results_df,
        ggplot2::aes(
            x = Impact,
            y = Log_P_value,
            label = Pathway
        )
    ) +
        ggplot2::geom_point(
            alpha = 0.85,
            size = 4,
            color = "#1F78B4"
        ) +
        ggrepel::geom_text_repel(
            size = 4.5,
            max.overlaps = 20,
            box.padding = 0.35,
            point.padding = 0.3,
            min.segment.length = 0.2,
            show.legend = FALSE
        ) +
        ggplot2::labs(
            title = "Pathway Impact vs Statistical Significance",
            x = "Pathway Impact",
            y = expression(-log[10](italic(P)~value))
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            legend.position = "none",
            panel.border = ggplot2::element_rect(
                color = "black",
                fill = NA,
                linewidth = 0.8
            ),
            plot.title = ggplot2::element_text(
                face = "bold",
                hjust = 0.5,
                size = 14
            )
        )
}
create_impact_plot <- function(significant_results_df) {
    
    if (nrow(significant_results_df) == 0) {
        warning("No pathways passed the p value cutoff.")
        return(NULL)
    }
    
    ggplot2::ggplot(
        significant_results_df,
        ggplot2::aes(
            x = Impact,
            y = Log_P_value,
            label = Pathway
        )
    ) +
        ggplot2::geom_point(
            alpha = 0.85,
            size = 4,
            color = "#1F78B4"
        ) +
        ggrepel::geom_text_repel(
            size = 4.5,
            max.overlaps = 20,
            box.padding = 0.35,
            point.padding = 0.3,
            min.segment.length = 0.2,
            show.legend = FALSE
        ) +
        ggplot2::labs(
            title = "Pathway Impact vs Statistical Significance",
            x = "Pathway Impact",
            y = "-log10(P value)"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            legend.position = "none",
            panel.border = ggplot2::element_rect(
                color = "black",
                fill = NA,
                linewidth = 0.8
            ),
            plot.title = element_text(
                face = "bold",
                hjust = 0.5,
                size = 14
            )
        )
}
