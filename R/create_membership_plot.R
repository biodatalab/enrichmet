#' Create Pathway Membership Plot
#'
#' Generates a heatmap showing pathway membership for metabolites.
#'
#' @param PathwayVsMetabolites A data frame with pathway-metabolite associations.
#' @param inputMetabolites A character vector of metabolite IDs.
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping.
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
#' plot <- create_membership_plot(PathwayVsMetabolites, inputMetabolites, kegg_lookup)
#' # Note: ComplexHeatmap objects need to be drawn explicitly
#' if (!is.null(plot)) {
#'   ComplexHeatmap::draw(plot)
#' }
#' @importFrom dplyr filter rename
#' @importFrom tidyr separate_rows
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid unit gpar
#' @export
create_membership_plot <- function(PathwayVsMetabolites, inputMetabolites, kegg_lookup = NULL) {
    
    # Prepare long format
    long_pathway_df <- PathwayVsMetabolites %>%
        tidyr::separate_rows(Metabolites, sep = ",") %>%
        dplyr::rename(pathway_name = Pathway, metabolite = Metabolites)
    
    # Filter for input metabolites
    matched_df <- long_pathway_df %>%
        dplyr::filter(metabolite %in% inputMetabolites)
    
    if (nrow(matched_df) == 0) {
        warning("No matching metabolites found for membership plot")
        return(NULL)
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
    
    if (nrow(membership_matrix) == 0 || ncol(membership_matrix) == 0) {
        warning("Membership matrix is empty")
        return(NULL)
    }
    
    # Create the heatmap object
    membership_plot <- ComplexHeatmap::Heatmap(
        membership_matrix,
        name = "Membership",
        col = c("0" = "white", "1" = "#2E86AB"),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10, rot = 45),
        row_title = "Metabolites",
        row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        column_title = "Pathway Membership",
        column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(
            title = "Member",
            title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 9),
            legend_direction = "vertical"
        ),
        border = TRUE,
        rect_gp = grid::gpar(col = "white", lwd = 1)
    )
    
    return(membership_plot)
}