#' Create GSEA Plot
#'
#' Generates a visualization of GSEA results.
#'
#' @param gsea_results Data frame from perform_gsea_analysis().
#' @param top_n Number of top pathways to display (default = 20).
#' @param kegg_lookup Optional data frame for mapping pathway IDs to names.
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
#' @importFrom ggplot2 labs theme_minimal theme element_text scale_size_continuous
#' @importFrom dplyr arrange left_join
#' @importFrom utils head
#' @export
create_gsea_plot <- function(gsea_results, top_n = 20, kegg_lookup = NULL) {
    if (nrow(gsea_results) == 0) {
        warning("No GSEA results to plot")
        return(NULL)
    }
    
    # Handle both old and new column names for compatibility
    # New GSEA output uses 'pval', old example used 'pval'
    if (!"pval" %in% colnames(gsea_results) && "pval" %in% colnames(gsea_results)) {
        gsea_results$pval <- gsea_results$pval
    }
    
    # Ensure required columns exist
    required_cols <- c("pathway", "pval", "NES")
    missing_cols <- setdiff(required_cols, colnames(gsea_results))
    if (length(missing_cols) > 0) {
        stop("Missing required columns in GSEA results: ", paste(missing_cols, collapse = ", "))
    }
    
    # Handle input_count column - use leadingEdge count if input_count doesn't exist
    if (!"input_count" %in% colnames(gsea_results)) {
        if ("leadingEdge" %in% colnames(gsea_results)) {
            gsea_results$input_count <- vapply(
                gsea_results$leadingEdge, 
                length, 
                FUN.VALUE = integer(1)
            )
        } else {
            # If neither exists, create a dummy column with value 1
            gsea_results$input_count <- 1
            warning("No 'input_count' or 'leadingEdge' column found. Using default size of 1 for all points.")
        }
    }
    
    # Apply KEGG lookup if provided - FIXED COLUMN NAMES
    if (!is.null(kegg_lookup)) {
        # Try different possible column name combinations
        if (all(c("pathway_id", "pathway_name") %in% colnames(kegg_lookup))) {
            # Map pathway IDs to human-readable names
            gsea_results <- gsea_results %>%
                dplyr::left_join(kegg_lookup, by = c("pathway" = "pathway_id")) %>%
                dplyr::mutate(
                    pathway = ifelse(!is.na(pathway_name), pathway_name, pathway)
                )
            message("Applied KEGG pathway name mapping using 'pathway_id' and 'pathway_name' columns")
        } else if (all(c("kegg_id", "name") %in% colnames(kegg_lookup))) {
            # Alternative column names that match your other functions
            gsea_results <- gsea_results %>%
                dplyr::left_join(kegg_lookup, by = c("pathway" = "kegg_id")) %>%
                dplyr::mutate(
                    pathway = ifelse(!is.na(name), name, pathway)
                )
            message("Applied KEGG pathway name mapping using 'kegg_id' and 'name' columns")
        } else if ("pathway" %in% colnames(kegg_lookup) && "name" %in% colnames(kegg_lookup)) {
            # Another possible column name combination
            gsea_results <- gsea_results %>%
                dplyr::left_join(kegg_lookup, by = "pathway") %>%
                dplyr::mutate(
                    pathway = ifelse(!is.na(name), name, pathway)
                )
            message("Applied KEGG pathway name mapping using 'pathway' and 'name' columns")
        } else {
            warning("kegg_lookup provided but missing required columns. Expected one of:\n",
                    "  - 'pathway_id' and 'pathway_name'\n", 
                    "  - 'kegg_id' and 'name'\n",
                    "  - 'pathway' and 'name'\n",
                    "Available columns in kegg_lookup: ", paste(colnames(kegg_lookup), collapse = ", "))
            
            # Debug information
            message("First few pathway names in GSEA results: ", paste(head(gsea_results$pathway), collapse = ", "))
            if (nrow(kegg_lookup) > 0) {
                message("First few entries in kegg_lookup: ", paste(head(kegg_lookup[,1]), collapse = ", "))
            }
        }
    }
    
    # Sort and prepare data - maintaining original format
    MSEAres <- gsea_results %>%
        dplyr::arrange(pval) %>%
        utils::head(top_n)
    
    # Check if we have any results after filtering
    if (nrow(MSEAres) == 0) {
        warning("No pathways to plot after filtering to top_n = ", top_n)
        return(NULL)
    }
    
    MSEAres$pathway <- factor(MSEAres$pathway, levels = rev(MSEAres$pathway))
    
    # Create the exact same plot as original
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
            legend.position = "right",   
            panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8)
        )
}