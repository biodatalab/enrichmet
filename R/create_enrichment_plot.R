
#' Create Pathway Enrichment Plot
#'
#' Generates a visualization of pathway enrichment results.
#'
#' @param significant_results_df Data frame from perform_enrichment_analysis().
#'
#' @return A ggplot object showing pathway enrichment results.
#'
#' Create Enrichment Plot
#'
#' Generates a visualization of pathway enrichment results.
#'
#' @param enrichment_results Data frame with enrichment analysis results.
#'
#' @return A ggplot object showing enrichment results.
#'
#' @examples
#' # Create comprehensive example enrichment results
#' set.seed(123)
#' enrichment_results <- data.frame(
#'   Pathway = c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)",
#'               "Pentose phosphate pathway", "Pyruvate metabolism",
#'               "Amino sugar and nucleotide sugar metabolism",
#'               "Glycine, serine and threonine metabolism",
#'               "Cysteine and methionine metabolism",
#'               "Valine, leucine and isoleucine degradation",
#'               "Fatty acid degradation", "Oxidative phosphorylation",
#'               "Purine metabolism", "Pyrimidine metabolism",
#'               "Alanine, aspartate and glutamate metabolism",
#'               "Butanoate metabolism", "Propanoate metabolism"),
#'   Log_P_value = -log10(c(1.2e-8, 3.4e-6, 2.1e-5, 7.8e-5, 1.5e-4,
#'                         3.2e-4, 8.7e-4, 0.0012, 0.0023, 0.0041,
#'                         0.0067, 0.0089, 0.012, 0.018, 0.025)),
#'   Adjusted_P_value = c(1.2e-6, 1.7e-4, 7.0e-4, 0.0013, 0.0024,
#'                       0.0040, 0.0087, 0.012, 0.018, 0.027,
#'                       0.035, 0.042, 0.048, 0.060, 0.075),
#'   Count = c(15, 12, 10, 9, 8, 7, 7, 6, 6, 5, 8, 7, 6, 5, 4)
#' )
#'
#' # Create enrichment plot
#' plot <- create_enrichment_plot(enrichment_results)
#' plot
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_gradient scale_size_continuous 
#' @importFrom ggplot2 labs theme_minimal theme element_text coord_flip
#' @importFrom stats reorder
#' @export
create_enrichment_plot <- function(enrichment_results) {
    if (is.null(enrichment_results) || nrow(enrichment_results) == 0) {
        stop("No enrichment results provided for plotting.")
    }
    
    # Ensure Count is numeric
    enrichment_results$Count <- as.numeric(enrichment_results$Count)
    
    # Detect min and max automatically
    min_count <- min(enrichment_results$Count, na.rm = TRUE)
    max_count <- max(enrichment_results$Count, na.rm = TRUE)
    
    # Plot
    p <- ggplot2::ggplot(
        enrichment_results,
        ggplot2::aes(
            x = reorder(Pathway, Log_P_value),
            y = Log_P_value,
            size = Count,
            color = Adjusted_P_value
        )
    ) +
        ggplot2::geom_point(alpha = 0.8) +
        ggplot2::scale_size_continuous(
            range = c(3, 12),
            breaks = pretty(c(min_count, max_count), n = 4),
            name = "Metabolite Count"
        ) +
        ggplot2::scale_color_gradient(
            low = "#D32F2F",        # small adj. p-values (more significant)
            high = "#1976D2",       # large adj. p-values (less significant)
            name = "Adj. p-value",
            trans = "log10"
        ) +
        ggplot2::labs(
            title = "Pathway Enrichment Analysis",
            x = "Pathway",
            y = "-log10(P-value)"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            axis.text.y = element_text(size = 10, color = "black"),
            axis.text.x = element_text(size = 10, color = "black"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            legend.position = "right"
        ) +
        ggplot2::coord_flip()
    
    return(p)
}