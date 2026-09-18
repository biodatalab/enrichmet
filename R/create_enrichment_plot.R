#' Create Pathway Enrichment Plot
#'
#' Generates a visualization of pathway enrichment results.
#'
#' @param enrichment_results Data frame with enrichment analysis results.
#'   Must contain columns \code{Pathway}, \code{Log_P_value}, and \code{Count}.
#'
#' @return A \code{ggplot} object showing enrichment results.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous
#' @importFrom ggplot2 labs theme_minimal theme element_text element_rect coord_flip
#' @importFrom stats reorder
#'
#' @examples
#' 
#' # Always-runnable minimal example (no network, no dependencies)
#' enr <- data.frame(
#'     Pathway = c("Glycolysis", "TCA cycle", "Pentose phosphate pathway"),
#'     Log_P_value = c(3.2, 2.1, 1.4),
#'     Count = c(8, 5, 3),
#'     stringsAsFactors = FALSE
#' )
#' create_enrichment_plot(enr)
#'
#' 
#' PathwayVsMetabolites <- fetch_kegg_pathway_metabolites(organism = "hsa")
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
#' enr_real <- perform_enrichment_analysis(
#'     inputMetabolites = query,
#'     PathwayVsMetabolites = PathwayVsMetabolites,
#'     top_n = 15,
#'     p_value_cutoff = 1
#' )
#'
#' create_enrichment_plot(enr_real)
#' 
#'
#' @export
create_enrichment_plot <- function(enrichment_results) {
    
    if (is.null(enrichment_results) || nrow(enrichment_results) == 0) {
        stop("No enrichment results provided for plotting.")
    }
    
    # Ensure Count is numeric and rounded
    enrichment_results$Count <- as.numeric(enrichment_results$Count)
    enrichment_results$Count <- round(enrichment_results$Count)
    
    # Determine size scale breaks from observed values
    unique_counts <- sort(unique(enrichment_results$Count))
    count_breaks <- unique_counts
    
    if (length(count_breaks) > 6) {
        idx <- round(seq(1, length(unique_counts), length.out = 6))
        count_breaks <- unique_counts[idx]
    }
    
    p <- ggplot2::ggplot(
        enrichment_results,
        ggplot2::aes(
            x = reorder(Pathway, Log_P_value),
            y = Log_P_value,
            size = Count
        )
    ) +
        ggplot2::geom_point(
            alpha = 0.8,
            color = "#1F78B4"
        ) +
        ggplot2::scale_size_continuous(
            range = c(3, 12),
            breaks = count_breaks,
            name = "Count"
        ) +
        ggplot2::labs(
            title = "Pathway Enrichment Analysis",
            x = "Pathway",
            y = "-log10(P-value)"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            axis.text.y = element_text(size = 12, color = "black"),
            axis.text.x = element_text(size = 12, color = "black"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            legend.position = "right",
            panel.border = ggplot2::element_rect(
                color = "black",
                fill = NA,
                linewidth = 0.8
            )
        ) +
        ggplot2::coord_flip()
    
    return(p)
}
