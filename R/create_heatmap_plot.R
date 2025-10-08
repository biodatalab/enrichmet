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
#'   name = c("D-Glucose", "L-Lactate", "Oxaloacetate", "Citrate", "D-Fructose 6-phosphate",
#'            "3-Phospho-D-glyceroyl phosphate", "L-Aspartate", "Oxoglutaric acid",
#'            "Acetyl-CoA", "Oxaloacetic acid", "L-Glutamate", "Succinyl-CoA",
#'            "L-Malate", "D-Ribose 5-phosphate", "D-Xylulose 5-phosphate",
#'            "Sedoheptulose 7-phosphate", "Glycine", "L-Serine",
#'            "Succinate", "Isocitrate", "Glyceraldehyde 3-phosphate",
#'            "D-Fructose 1,6-bisphosphate", "D-Erythrose 4-phosphate",
#'            "L-Threonine", "L-Cysteine", "O-Acetyl-L-serine",
#'            "L-Cystathionine", "L-Homocysteine", "S-Adenosyl-L-methionine")
#' )
#'
#' # Create enrichment heatmap
#' heatmap_plot <- create_heatmap_plot(enrichment_results, PathwayVsMetabolites, 
#'                                    inputMetabolites, kegg_lookup)
#' # Note: ComplexHeatmap objects need to be drawn explicitly
#' if (!is.null(heatmap_plot)) {
#'   ComplexHeatmap::draw(heatmap_plot)
#' }
#'
#' @importFrom dplyr filter mutate
#' @importFrom tidyr unnest
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @export
create_heatmap_plot <- function(enrichment_results, PathwayVsMetabolites, inputMetabolites, kegg_lookup = NULL) {
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
    
    if (!is.null(kegg_lookup)) {
        data_filtered <- data_filtered %>%
            dplyr::left_join(kegg_lookup, by = c("Metabolites" = "kegg_id")) %>%
            dplyr::mutate(Metabolites = ifelse(!is.na(name), name, Metabolites))
    }
    
    heatmap_matrix <- table(data_filtered$Metabolites, data_filtered$Pathway)
    heatmap_matrix <- as.matrix(heatmap_matrix)
    
    logp_vec <- enrichment_results$Log_P_value
    names(logp_vec) <- enrichment_results$Pathway
    
    common_pathways <- intersect(colnames(heatmap_matrix), names(logp_vec))
    if (length(common_pathways) == 0) {
        warning("No common pathways found between enrichment results and pathway data")
        return(NULL)
    }
    
    heatmap_matrix <- heatmap_matrix[, common_pathways, drop = FALSE]
    logp_vec <- logp_vec[common_pathways]
    
    heatmap_values <- sweep(heatmap_matrix, 2, logp_vec, `*`)
    
    value_range <- range(heatmap_values[heatmap_values > 0])
    color_breaks <- seq(0, ceiling(max(value_range)), length.out = 4)
    
    white_blue_red_gradient <- c("#FFFFFF", "#6BAED6", "#EF3B2C", "#67000D")
    
    # Calculate dynamic font sizes based on data size
    n_rows <- nrow(heatmap_values)
    n_cols <- ncol(heatmap_values)
    
    # Adjust font sizes based on matrix size
    row_fontsize <- ifelse(n_rows > 30, 8, ifelse(n_rows > 15, 10, 12))
    col_fontsize <- ifelse(n_cols > 20, 8, ifelse(n_cols > 10, 10, 12))
    
    heatmap_plot <- ComplexHeatmap::Heatmap(
        heatmap_values,
        name = "-log10(P-value)",
        col = circlize::colorRamp2(color_breaks, white_blue_red_gradient),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_side = "left",
        column_names_rot = 45,
        border = TRUE,
        rect_gp = grid::gpar(col = "white", lwd = 0.5),
        row_names_gp = grid::gpar(fontsize = row_fontsize),
        column_names_gp = grid::gpar(fontsize = col_fontsize),
        column_title = "Metabolite-Pathway Enrichment Significance",
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
        row_title = "Metabolites",
        row_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
        heatmap_legend_param = list(
            title = "Enrichment\nScore",
            title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 10),
            legend_height = grid::unit(4, "cm")
        )
        # REMOVED: width and height parameters - let it auto-size
    )
    
    return(heatmap_plot)
}