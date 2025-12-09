#' Create Enrichment Heatmap
#'
#' Generates a heatmap showing metabolite-pathway enrichment significance.
#'
#' @param enrichment_results Data frame from perform_enrichment_analysis().
#' @param PathwayVsMetabolites A data frame with pathway-metabolite associations.
#' @param inputMetabolites A character vector of metabolite IDs.
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping.
#'
#' @return A ComplexHeatmap object showing enrichment significance.
#'
#' @examples
#' # Create Enrichment Heatmap
#'
#' # Generates a heatmap showing metabolite-pathway enrichment significance.
#'
#' @param enrichment_results Data frame from perform_enrichment_analysis().
#' @param PathwayVsMetabolites A data frame with pathway-metabolite associations.
#' @param inputMetabolites A character vector of metabolite IDs.
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping.
#' @param min_pathways Minimum number of metabolites a pathway must have to be included.
#' @param min_metabolites Minimum number of pathways a metabolite must be in to be included.
#' @param base_font_size Base font size for the heatmap labels.
#' @param top_n Number of top metabolites to include (by significance).
#'
#' @return A ComplexHeatmap object showing enrichment significance.
#'
#' @examples
#' # Create comprehensive example data
#' set.seed(123)
#' enrichment_results <- data.frame(
#'   Pathway = c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)",
#'               "Pentose phosphate pathway", "Pyruvate metabolism",
#'               "Amino sugar and nucleotide sugar metabolism",
#'               "Glycine, serine and threonine metabolism",
#'               "Cysteine and methionine metabolism"),
#'   Log_P_value = c(7.2, 5.8, 4.9, 4.3, 3.8, 3.5, 3.2),
#'   P_value = 10^(-c(7.2, 5.8, 4.9, 4.3, 3.8, 3.5, 3.2)),
#'   Adjusted_P_value = 10^(-c(6.5, 5.2, 4.4, 3.9, 3.5, 3.2, 2.9))
#' )
#'
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)",
#'               "Pentose phosphate pathway", "Pyruvate metabolism",
#'               "Amino sugar and nucleotide sugar metabolism",
#'               "Glycine, serine and threonine metabolism",
#'               "Cysteine and methionine metabolism"),
#'   Metabolites = c("C00031,C00022,C00074,C00036,C00103,C00197,C00186,C00221,C01172",
#'                   "C00024,C00036,C00042,C00122,C00149,C00158,C00417,C05379",
#'                   "C00117,C00231,C00279,C00345,C01172,C01236,C01218",
#'                   "C00022,C00024,C00036,C00074,C00122,C00149,C00232",
#'                   "C00031,C00095,C00103,C00121,C00140,C00208,C00267",
#'                   "C00037,C00048,C00065,C00143,C00258,C00300,C00497",
#'                   "C00021,C00073,C00094,C00155,C00283,C00542,C00957")
#' )
#'
#' inputMetabolites <- c("C00031", "C00022", "C00074", "C00036", "C00103", "C00197",
#'                      "C00186", "C00221", "C00024", "C00042", "C00122", "C00149",
#'                      "C00158", "C00117", "C00231", "C00279", "C00037", "C00048")
#'
#' kegg_lookup <- data.frame(
#'   kegg_id = c("C00031", "C00022", "C00074", "C00036", "C00103", "C00197",
#'               "C00186", "C00221", "C00024", "C00042", "C00122", "C00149",
#'               "C00158", "C00117", "C00231", "C00279", "C00037", "C00048",
#'               "C00417", "C05379", "C00345", "C01172", "C01236", "C00065",
#'               "C00143", "C00258", "C00021", "C00073", "C00094"),
#'   name = c("D-Glucose", "L-Lactate", "Oxaloacetate", "Citrate",
#'            "D-Fructose 6-phosphate", "3-Phospho-D-glyceroyl phosphate",
#'            "L-Aspartate", "Oxoglutaric acid", "Acetyl-CoA",
#'            "Oxaloacetic acid", "L-Glutamate", "Succinyl-CoA",
#'            "L-Malate", "D-Ribose 5-phosphate", "D-Xylulose 5-phosphate",
#'            "Sedoheptulose 7-phosphate", "Glycine", "L-Serine",
#'            "Succinate", "Isocitrate", "Glyceraldehyde 3-phosphate",
#'            "D-Fructose 1,6-bisphosphate", "D-Erythrose 4-phosphate",
#'            "L-Threonine", "L-Cysteine", "O-Acetyl-L-serine",
#'            "L-Cystathionine", "L-Homocysteine", "S-Adenosyl-L-methionine")
#' )
#'
#' # Create enrichment heatmap with all parameters
#' heatmap_plot <- create_heatmap_plot(
#'   enrichment_results = enrichment_results,
#'   PathwayVsMetabolites = PathwayVsMetabolites,
#'   inputMetabolites = inputMetabolites,
#'   kegg_lookup = kegg_lookup,
#'   top_n = 10,
#'   min_pathways = 1,
#'   min_metabolites = 1
#' )
#'
#' # Note: ComplexHeatmap objects need to be drawn explicitly
#' if (!is.null(heatmap_plot)) {
#'   ComplexHeatmap::draw(heatmap_plot)
#' }
#'
#'
#' @importFrom dplyr filter mutate
#' @importFrom tidyr unnest
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @export
create_heatmap_plot <- function(enrichment_results, PathwayVsMetabolites, inputMetabolites, 
                                kegg_lookup = NULL, min_pathways = 2, min_metabolites = 1,
                                base_font_size = 12, top_n = NULL) {
    
    if (nrow(enrichment_results) == 0) {
        warning("No significant pathways found for heatmap")
        return(NULL)
    }
    
    data_filtered <- PathwayVsMetabolites %>%
        dplyr::mutate(Metabolites = strsplit(Metabolites, ",")) %>%
        tidyr::unnest(Metabolites) %>%
        dplyr::filter(Pathway %in% enrichment_results$Pathway,
                      Metabolites %in% inputMetabolites)
    
    if (nrow(data_filtered) == 0) {
        warning("No matching metabolites found for the significant pathways")
        return(NULL)
    }
    
    # Select top N metabolites based on enrichment significance
    if (!is.null(top_n)) {
        metabolite_significance <- data_filtered %>%
            dplyr::left_join(enrichment_results %>% 
                                 dplyr::select(Pathway, Adjusted_P_value), 
                             by = "Pathway") %>%
            dplyr::group_by(Metabolites) %>%
            dplyr::summarise(
                Min_Adjusted_P = min(Adjusted_P_value, na.rm = TRUE),
                Mean_NegLog_P = mean(-log10(Adjusted_P_value), na.rm = TRUE),
                .groups = "drop"
            ) %>%
            dplyr::arrange(Min_Adjusted_P) %>%
            head(top_n)
        
        data_filtered <- data_filtered %>%
            dplyr::filter(Metabolites %in% metabolite_significance$Metabolites)
        
        message("Selected top ", nrow(metabolite_significance), 
                " metabolites by enrichment significance for heatmap")
    }
    
    if (!is.null(kegg_lookup)) {
        data_filtered <- data_filtered %>%
            dplyr::left_join(kegg_lookup, by = c("Metabolites" = "kegg_id")) %>%
            dplyr::mutate(Metabolites = ifelse(!is.na(name), name, Metabolites))
    }
    
    heatmap_matrix <- table(data_filtered$Metabolites, data_filtered$Pathway)
    heatmap_matrix <- as.matrix(heatmap_matrix)
    
    # Apply filtering based on user input
    pathway_counts <- colSums(heatmap_matrix)
    metabolite_counts <- rowSums(heatmap_matrix)
    
    # Filter pathways and metabolites
    keep_pathways <- pathway_counts >= min_pathways
    keep_metabolites <- metabolite_counts >= min_metabolites
    
    if (sum(keep_pathways) == 0 || sum(keep_metabolites) == 0) {
        warning("No pathways or metabolites meet the filtering criteria for heatmap")
        return(NULL)
    }
    
    filtered_matrix <- heatmap_matrix[keep_metabolites, keep_pathways, drop = FALSE]
    
    logp_vec <- enrichment_results$Log_P_value
    names(logp_vec) <- enrichment_results$Pathway
    
    common_pathways <- intersect(colnames(filtered_matrix), names(logp_vec))
    if (length(common_pathways) == 0) {
        warning("No common pathways found between enrichment results and pathway data")
        return(NULL)
    }
    
    filtered_matrix <- filtered_matrix[, common_pathways, drop = FALSE]
    logp_vec <- logp_vec[common_pathways]
    
    heatmap_values <- sweep(filtered_matrix, 2, logp_vec, `*`)
    
    value_range <- range(heatmap_values[heatmap_values > 0])
    
    # Color gradient
    if (length(value_range) > 1 && !any(is.na(value_range)) && value_range[2] > value_range[2]) {
        color_breaks <- seq(value_range[1], value_range[2], length.out = 3)
        white_blue_red_gradient <- c("#E0F4FF","#FFA07A", "#B2182B")
    } else {
        color_breaks <- c(0, 1.5, 3)
        white_blue_red_gradient <- c("#E0F4FF","#FFA07A", "#B2182B")
    }
    
    # Calculate dynamic font sizes
    n_rows <- nrow(heatmap_values)
    n_cols <- ncol(heatmap_values)
    
    row_fontsize <- ifelse(n_rows > 30, base_font_size - 4, 
                           ifelse(n_rows > 15, base_font_size - 2, base_font_size))
    col_fontsize <- ifelse(n_cols > 20, base_font_size - 4, 
                           ifelse(n_cols > 10, base_font_size - 2, base_font_size))
    
    # Title now matches actual filtering
    plot_title <- ifelse(!is.null(top_n), 
                         paste("Top", top_n, "Metabolites"))
    
    heatmap_plot <- ComplexHeatmap::Heatmap(
        heatmap_values,
        name = "-log10(P-value)",
        col = circlize::colorRamp2(color_breaks, white_blue_red_gradient),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_side = "right",
        column_names_rot = 90,
        column_names_side = "bottom",
        border = TRUE,
        rect_gp = grid::gpar(col = "white", lwd = 0.5),
        row_names_gp = grid::gpar(fontsize = row_fontsize, fontface = "plain"),
        column_names_gp = grid::gpar(fontsize = col_fontsize, fontface = "plain"),
        column_title = plot_title,
        column_title_gp = grid::gpar(fontsize = base_font_size + 2, fontface = "bold"),
        row_title = "Metabolites",
        row_title_gp = grid::gpar(fontsize = base_font_size + 2, fontface = "bold"),
        heatmap_legend_param = list(
            title = "-log10(P-value)",
            title_gp = grid::gpar(fontsize = base_font_size, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = base_font_size - 2),
            legend_height = grid::unit(4, "cm")
        ),
        # Key parameters to make dendrogram touch the heatmap
        row_dend_width = grid::unit(0.5, "cm"),
        column_dend_height = grid::unit(0.5, "cm"),
        row_dend_side = "left",
        column_dend_side = "top",
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        gap = grid::unit(0, "mm")
    )
    
    return(heatmap_plot)
}