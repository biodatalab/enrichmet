#' Create Pathway Membership Plot
#'
#' Generates a heatmap showing pathway membership for metabolites.
#'
#' @param PathwayVsMetabolites A data frame with pathway-metabolite associations.
#' @param inputMetabolites A character vector of metabolite IDs.
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping.
#' @param top_n Number of top metabolites to display (by centrality).
#' @param min_pathway_occurrence Minimum number of metabolites a pathway must have to be included.
#' @param min_metabolite_occurrence Minimum number of pathways a metabolite must be in to be included.
#'
#' @return A ComplexHeatmap object showing pathway membership.
#'
#' @examples
#' # Create comprehensive example data with matching metabolite IDs
#' inputMetabolites <- c("C00031", "C00022", "C00074", "C00036", "C00103", 
#'                      "C00197", "C00186", "C00221", "C00024", "C00042",
#'                      "C00122", "C00149", "C00158", "C00117", "C00231",
#'                      "C00279", "C00037", "C00048", "C00065", "C00143")
#' 
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Glycolysis/Gluconeogenesis", "Citrate cycle (TCA cycle)", 
#'               "Pentose phosphate pathway", "Pyruvate metabolism",
#'               "Amino sugar and nucleotide sugar metabolism",
#'               "Glycine, serine and threonine metabolism",
#'               "Cysteine and methionine metabolism",
#'               "Valine, leucine and isoleucine degradation"),
#'   Metabolites = c("C00031,C00022,C00074,C00036,C00103,C00197,C00186,C00221",
#'                   "C00024,C00036,C00042,C00122,C00149,C00158,C00417",
#'                   "C00117,C00231,C00279,C00345,C01172",
#'                   "C00022,C00024,C00036,C00074,C00122",
#'                   "C00031,C00095,C00103,C00121,C00140",
#'                   "C00037,C00048,C00065,C00143,C00258",
#'                   "C00021,C00073,C00094,C00155,C00283",
#'                   "C00123,C00183,C00233,C00407,C00979")
#' )
#' 
#' # Create KEGG lookup with matching IDs
#' kegg_lookup <- data.frame(
#'   kegg_id = c("C00031", "C00022", "C00074", "C00036", "C00103", 
#'               "C00197", "C00186", "C00221", "C00024", "C00042",
#'               "C00122", "C00149", "C00158", "C00117", "C00231",
#'               "C00279", "C00037", "C00048", "C00065", "C00143",
#'               "C00258", "C00021", "C00073", "C00094", "C00155"),
#'   name = c("D-Glucose", "L-Lactate", "Oxaloacetate", "Citrate", 
#'            "D-Fructose 6-phosphate", "3-Phospho-D-glyceroyl phosphate",
#'            "L-Aspartate", "Oxoglutaric acid", "Acetyl-CoA", 
#'            "Oxaloacetic acid", "L-Glutamate", "Succinyl-CoA",
#'            "L-Malate", "D-Ribose 5-phosphate", "D-Xylulose 5-phosphate",
#'            "Sedoheptulose 7-phosphate", "Glycine", "L-Serine",
#'            "L-Threonine", "L-Cysteine", "O-Acetyl-L-serine",
#'            "L-Cystathionine", "L-Homocysteine", "S-Adenosyl-L-methionine", "Coenzyme A")
#' )
#'
#' # Create membership plot
#' plot <- create_membership_plot(PathwayVsMetabolites, inputMetabolites, kegg_lookup, top_n = 10)
#' # Note: ComplexHeatmap objects need to be drawn explicitly
#' if (!is.null(plot)) {
#'   ComplexHeatmap::draw(plot)
#' }
#' @importFrom dplyr filter rename
#' @importFrom tidyr separate_rows
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid unit gpar
#' @export
create_membership_plot <- function(PathwayVsMetabolites, inputMetabolites, kegg_lookup = NULL, 
                                   top_n = NULL, min_pathway_occurrence = 1, min_metabolite_occurrence = 1) {
    
    message("=== MEMBERSHIP PLOT DEBUG ===")
    message("Input metabolites: ", length(inputMetabolites))
    message("Top_n parameter: ", top_n)
    
    # Prepare long format
    long_pathway_df <- PathwayVsMetabolites %>%
        tidyr::separate_rows(Metabolites, sep = ",") %>%
        dplyr::rename(pathway_name = Pathway, metabolite = Metabolites)
    
    # Filter for input metabolites
    matched_df <- long_pathway_df %>%
        dplyr::filter(metabolite %in% inputMetabolites)
    
    message("Initial matched metabolites: ", length(unique(matched_df$metabolite)))
    
    if (nrow(matched_df) == 0) {
        warning("No matching metabolites found for membership plot")
        return(NULL)
    }
    
    # Filter to top N metabolites by centrality if specified
    if (!is.null(top_n)) {
        message("Filtering to top ", top_n, " metabolites by centrality...")
        
        # Calculate centrality to identify top metabolites
        all_centrality <- calculate_metabolite_centrality(PathwayVsMetabolites)
        message("Total metabolites with centrality: ", nrow(all_centrality))
        
        top_metabolites <- all_centrality %>%
            dplyr::filter(Metabolite %in% inputMetabolites) %>%
            dplyr::arrange(desc(RBC_Metabolite)) %>%
            head(top_n) %>%
            dplyr::pull(Metabolite)
        
        message("Top metabolites selected: ", length(top_metabolites))
        message("Sample top metabolites: ", paste(head(top_metabolites), collapse = ", "))
        
        matched_df <- matched_df %>%
            dplyr::filter(metabolite %in% top_metabolites)
        
        message("After centrality filtering - unique metabolites: ", length(unique(matched_df$metabolite)))
    }
    
    # Map KEGG IDs to names if lookup provided
    if (!is.null(kegg_lookup)) {
        matched_df <- matched_df %>%
            dplyr::left_join(kegg_lookup, by = c("metabolite" = "kegg_id")) %>%
            dplyr::mutate(name = ifelse(is.na(name), metabolite, name))
    } else {
        matched_df <- matched_df %>% dplyr::mutate(name = metabolite)
    }
    
    # Create membership matrix
    membership_matrix <- table(matched_df$name, matched_df$pathway_name)
    membership_matrix <- as.matrix(membership_matrix)
    membership_matrix[membership_matrix > 1] <- 1
    
    message("Matrix before occurrence filtering: ", nrow(membership_matrix), " x ", ncol(membership_matrix))
    
    if (nrow(membership_matrix) == 0 || ncol(membership_matrix) == 0) {
        warning("Membership matrix is empty")
        return(NULL)
    }
    
    # Apply pathway and metabolite occurrence filtering
    pathway_counts <- colSums(membership_matrix)
    metabolite_counts <- rowSums(membership_matrix)
    
    message("Pathway counts range: ", min(pathway_counts), " to ", max(pathway_counts))
    message("Metabolite counts range: ", min(metabolite_counts), " to ", max(metabolite_counts))
    
    # Filter pathways based on minimum occurrence
    keep_pathways <- pathway_counts >= min_pathway_occurrence
    keep_metabolites <- metabolite_counts >= min_metabolite_occurrence
    
    message("Pathways to keep: ", sum(keep_pathways), "/", length(keep_pathways))
    message("Metabolites to keep: ", sum(keep_metabolites), "/", length(keep_metabolites))
    
    if (sum(keep_pathways) == 0 || sum(keep_metabolites) == 0) {
        warning("No pathways or metabolites meet the filtering criteria for membership plot")
        return(NULL)
    }
    
    filtered_matrix <- membership_matrix[keep_metabolites, keep_pathways, drop = FALSE]
    
    message("Final matrix dimensions: ", nrow(filtered_matrix), " x ", ncol(filtered_matrix))
    message("=== END DEBUG ===")
    
    # Rest of the function remains the same...
    # Update title to show filtering information
    plot_title <- "Pathway Membership"
    if (!is.null(top_n)) {
        plot_title <- paste(plot_title, "(Top", top_n, "Metabolites by Centrality)")
    }
    if (min_pathway_occurrence > 1 || min_metabolite_occurrence > 1) {
        plot_title <- paste0(plot_title, "\nPathways with >=", min_pathway_occurrence, 
                             " metabolites | Metabolites in >=", min_metabolite_occurrence, " pathways")
    }
    
    # Calculate dynamic font sizes for better visualization
    n_rows <- nrow(filtered_matrix)
    n_cols <- ncol(filtered_matrix)
    
    # Adjust font sizes based on matrix size
    row_fontsize <- ifelse(n_rows > 30, 8, 
                           ifelse(n_rows > 15, 9, 10))
    col_fontsize <- ifelse(n_cols > 20, 8, 
                           ifelse(n_cols > 10, 9, 10))
    
    # Create the heatmap object
    membership_plot <- ComplexHeatmap::Heatmap(
        filtered_matrix,
        name = "Membership",
        col = c("0" = "white", "1" = "#2E86AB"),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = grid::gpar(fontsize = row_fontsize),
        column_names_gp = grid::gpar(fontsize = col_fontsize, rot = 45),
        row_title = "Metabolites",
        row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        column_title = plot_title,
        column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(
            title = "Member",
            title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 10),
            legend_direction = "vertical"
        ),
        border = TRUE,
        rect_gp = grid::gpar(col = "white", lwd = 1)
    )
    
    return(membership_plot)
}