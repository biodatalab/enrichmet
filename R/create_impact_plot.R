#' Create Impact vs Significance Plot
#'
#' Generates a visualization of pathway impact versus statistical significance.
#'
#' @param significant_results_df Data frame from perform_enrichment_analysis().
#'
#' @return A ggplot object showing impact vs significance.
#'
#' @examples
#' # Create comprehensive example enrichment results
#' set.seed(123)
#' significant_results_df <- data.frame(
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
#'   Impact = c(0.85, 0.72, 0.68, 0.61, 0.55, 0.48, 0.45, 0.42, 0.38, 0.35,
#'              0.52, 0.47, 0.41, 0.36, 0.32),
#'   Log_P_value = -log10(c(1.2e-8, 3.4e-6, 2.1e-5, 7.8e-5, 1.5e-4,
#'                         3.2e-4, 8.7e-4, 0.0012, 0.0023, 0.0041,
#'                         0.0067, 0.0089, 0.012, 0.018, 0.025)),
#'   P_value = c(1.2e-8, 3.4e-6, 2.1e-5, 7.8e-5, 1.5e-4, 3.2e-4, 8.7e-4,
#'               0.0012, 0.0023, 0.0041, 0.0067, 0.0089, 0.012, 0.018, 0.025),
#'   Adjusted_P_value = c(1.2e-6, 1.7e-4, 7.0e-4, 0.0013, 0.0024, 0.0040,
#'                       0.0087, 0.012, 0.018, 0.027, 0.035, 0.042, 0.048, 0.060, 0.075)
#' )
#'
#' # Create impact plot
#' plot <- create_impact_plot(significant_results_df)
#' plot
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient scale_size_continuous 
#' @importFrom ggplot2 labs theme_minimal theme
#' @importFrom ggrepel geom_text_repel
#' @export
create_impact_plot <- function(significant_results_df) {
    if (nrow(significant_results_df) == 0) {
        warning("No pathways passed the p-value cutoff.")
        return(NULL)
    }
    
    ggplot2::ggplot(
        significant_results_df,
        aes(
            x = Impact,
            y = Log_P_value,
            size = Impact,
            color = Adjusted_P_value,
            label = Pathway
        )
    ) +
        ggplot2::geom_point(alpha = 0.85) +
        ggrepel::geom_text_repel(
            size = 3.2, 
            max.overlaps = 20, 
            box.padding = 0.35,
            point.padding = 0.3,
            min.segment.length = 0.2,
            show.legend = FALSE
        ) +
        ggplot2::scale_size_continuous(
            range = c(3, 10),
            name = "Pathway Impact"
        ) +
        ggplot2::scale_color_gradient(
            low = "#D32F2F", 
            high = "#1976D2", 
            name = "Adj. p-value",
            trans = "log10"
        ) +
        ggplot2::labs(
            title = "Pathway Impact vs Statistical Significance",
            x = "Pathway Impact", 
            y = "-log10(P-value)"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "right",
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
        )
}