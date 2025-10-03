#' Create Centrality Plot
#'
#' Generates a visualization of metabolite centrality results.
#'
#' @param centrality_results Data frame from calculate_metabolite_centrality().
#' @param top_n Number of top metabolites to display (default = 20).
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping.
#'
#' @return A ggplot object showing metabolite centrality.
#'
#' @examples
#' # Create comprehensive example centrality results
#' set.seed(123)
#' centrality_results <- data.frame(
#'   Metabolite = paste0("C", sprintf("%05d", 1:50)),
#'   RBC_Metabolite = c(runif(30, 0.7, 0.95), runif(20, 0.3, 0.6))
#' )
#' 
#' # Create KEGG lookup for metabolite names
#' kegg_lookup <- data.frame(
#'   kegg_id = paste0("C", sprintf("%05d", 1:50)),
#'   name = c("Glucose", "Lactate", "Pyruvate", "Alanine", "Glutamate",
#'            "Citrate", "Succinate", "Malate", "Aspartate", "Glutamine",
#'            "Acetyl-CoA", "Oxaloacetate", "Alpha-ketoglutarate", "Fumarate",
#'            "Isocitrate", "Ribose-5P", "Glycerol-3P", "Dihydroxyacetone-P",
#'            "Fructose-6P", "Glucose-6P", "Phosphoenolpyruvate", "2-Phosphoglycerate",
#'            "3-Phosphoglycerate", "Glyceraldehyde-3P", "Sedioheptulose-7P",
#'            "Xylulose-5P", "Ribulose-5P", "Erythrose-4P", "Acetate", "Butyrate",
#'            paste0("Metabolite_", 31:50))
#' )
#'
#' # Create centrality plot
#' plot <- create_centrality_plot(centrality_results, top_n = 15, kegg_lookup = kegg_lookup)
#' plot
#'
#' @importFrom dplyr left_join mutate
#' @importFrom ggplot2 ggplot aes geom_col labs theme_minimal theme element_text coord_flip
#' @importFrom stats reorder
#' @importFrom utils head
#' @export
create_centrality_plot <- function(centrality_results, top_n = 20, kegg_lookup = NULL) {
    if (nrow(centrality_results) == 0) {
        warning("No centrality results to plot")
        return(NULL)
    }
    
    # Prepare data for plotting
    top_central <- head(centrality_results, top_n)
    
    # Add display names if available
    if (!is.null(kegg_lookup)) {
        top_central <- top_central %>%
            dplyr::left_join(kegg_lookup, by = c("Metabolite" = "kegg_id")) %>%
            dplyr::mutate(Display_Name = ifelse(!is.na(name), name, Metabolite))
    } else {
        top_central$Display_Name <- top_central$Metabolite
    }
    
    ggplot2::ggplot(
        top_central,
        aes(x = reorder(Display_Name, RBC_Metabolite), y = RBC_Metabolite)
    ) +
        ggplot2::geom_col(fill = "#2E86AB", alpha = 0.85) +
        ggplot2::labs(
            title = paste("Top", top_n, "Metabolites by Relative Betweenness Centrality"),
            subtitle = "Measures connectivity importance in metabolic network",
            x = "Metabolite", 
            y = "Relative Betweenness Centrality"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            axis.text.y = element_text(size = 10, color = "black"),
            axis.text.x = element_text(size = 10, color = "black"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11)
        ) +
        ggplot2::coord_flip()
}