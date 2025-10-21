#' Create GSEA Plot
#'
#' Generates a visualization of GSEA results.
#'
#' @param gsea_results Data frame from perform_gsea_analysis().
#' @param top_n Number of top pathways to display (default = 20).
#'
#' @return A ggplot object showing GSEA results.
#'
#' @examples
#' # Create comprehensive example GSEA results
#' set.seed(123)
#' gsea_results <- data.frame(
#'   pathway = c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)",
#'               "Pentose phosphate pathway", "Pyruvate metabolism",
#'               "Amino sugar and nucleotide sugar metabolism",
#'               "Glycine, serine and threonine metabolism",
#'               "Cysteine and methionine metabolism",
#'               "Valine, leucine and isoleucine degradation",
#'               "Fatty acid degradation", "Oxidative phosphorylation",
#'               "Purine metabolism", "Pyrimidine metabolism",
#'               "Alanine, aspartate and glutamate metabolism",
#'               "Butanoate metabolism", "Propanoate metabolism",
#'               "Starch and sucrose metabolism", "Galactose metabolism",
#'               "Fructose and mannose metabolism", "Glycerolipid metabolism",
#'               "Glycerophospholipid metabolism"),
#'   pval = c(1.2e-8, 3.4e-6, 2.1e-5, 7.8e-5, 1.5e-4, 3.2e-4, 8.7e-4,
#'            0.0012, 0.0023, 0.0041, 0.0067, 0.0089, 0.012, 0.018,
#'            0.025, 0.032, 0.041, 0.055, 0.067, 0.078),
#'   NES = c(2.3, 2.1, 1.9, 1.8, 1.7, 1.6, 1.5, -1.4, -1.3, 1.2,
#'           -1.1, 1.0, 0.9, -0.8, 0.7, 0.6, -0.5, 0.4, 0.3, -0.2),
#'   input_count = c(15, 12, 10, 9, 8, 7, 7, 6, 6, 5, 8, 7, 6, 5, 4, 7, 5, 6, 4, 3)
#' )
#'
#' # Create GSEA plot
#' plot <- create_gsea_plot(gsea_results, top_n = 15)
#' plot
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient2 
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @importFrom utils head
#' @export
create_gsea_plot <- function(gsea_results, top_n = 20) {
    if (nrow(gsea_results) == 0) {
        warning("No GSEA results to plot")
        return(NULL)
    }
    
    # Sort and prepare data
    MSEAres <- gsea_results %>%
        dplyr::arrange(pval) %>%
        head(top_n)
    
    MSEAres$pathway <- factor(MSEAres$pathway, levels = rev(MSEAres$pathway))
    
    ggplot2::ggplot(MSEAres,
                    ggplot2::aes(
                        x = -log10(pval),
                        y = pathway,
                        size = input_count,
                        color = NES
                    )) +
        ggplot2::geom_point(alpha = 0.8) +
        ggplot2::labs(
            title = "Metabolite Set Enrichment Analysis",
            x = "-log10(p-value)", 
            y = "Pathway",
            size = "Metabolite Count",
            color = "NES"
        ) +
        ggplot2::scale_color_gradient2(
            low = "#1E88E5", mid = "white", high = "#D81B60", 
            midpoint = 0,
            name = "NES"
        ) +
        ggplot2::scale_size_continuous(
            range = c(3, 10),
            name = "Metabolite Count"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            axis.text.y = element_text(size = 12, color = "black"),
            axis.text.x = element_text(size = 12, color = "black"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            legend.position = "right",   panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8)
        )
}