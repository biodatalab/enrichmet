#' Prepare GMT Data for Pathway Analysis
#'
#' Converts pathway-metabolite data to GMT format required for GSEA analysis.
#' Handles various metabolite ID formats and provides detailed filtering information.
#'
#' @param PathwayVsMetabolites A data frame containing pathways and their associated 
#'        metabolites. Must include columns 'Pathway' and 'Metabolites'.
#' @param metabolites_in_data Optional character vector of metabolite IDs to filter pathways.
#' @param min_pathway_size Minimum number of metabolites required in a pathway (default = 5).
#'
#' @return A list of pathways with their associated metabolites in GMT format.
#'
#' @examples
#' # Example 1: Basic usage with KEGG metabolites
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Glycolysis", "TCA Cycle", "Pentose Phosphate"),
#'   Metabolites = c(
#'     "C00031,C00022,C00036,C00103,C00111,C00118,C00186,C00197,C00221,C00236",
#'     "C00024,C00026,C00042,C00122,C00149,C00158,C00311,C00417,C05125",
#'     "C00117,C00118,C00121,C00231,C00279,C00345,C00620,C01151,C01236"
#'   )
#' )
#'
#' gmt_data <- prepare_gmt_data(PathwayVsMetabolites)
#' str(gmt_data)
#'
#' # Example 2: Filtering with specific KEGG IDs
#' my_metabolites <- c(
#'   "C00031", "C00022", "C00036", "C00103", "C00024", "C00026"
#' )
#' filtered_gmt <- prepare_gmt_data(PathwayVsMetabolites, my_metabolites)
#' str(filtered_gmt)
#'
#' # Example 3: Complex IDs with multiple separators
#' complex_pathways <- data.frame(
#'   Pathway = c("PathwayA", "PathwayB"),
#'   Metabolites = c(
#'     "C00031|C00022,C00036; C00103", 
#'     "C00024|C00026; C00042,C00149"
#'   )
#' )
#' complex_gmt <- prepare_gmt_data(complex_pathways)
#' str(complex_gmt)
#'
#' @examples
#' # Example 4: Realistic human metabolic pathways
#' human_pathways <- data.frame(
#'   Pathway = c(
#'     "Glycolysis / Gluconeogenesis", 
#'     "Citrate cycle (TCA cycle)",
#'     "Pentose phosphate pathway",
#'     "Fatty acid metabolism"
#'   ),
#'   Metabolites = c(
#'     paste0(
#'       "C00031,C00022,C00036,C00103,C00111,C00118,C00186,C00197,",
#'       "C00221,C00236,C00267,C00354,C00631,C00668,C01172"
#'     ),
#'     "C00024,C00026,C00036,C00042,C00122,C00149,C00158,C00311,C00417,C05125",
#'     paste0(
#'       "C00117,C00118,C00121,C00231,C00279,C00345,C00620,C01151,",
#'       "C01236,C01801,C05378,C05382"
#'     ),
#'     "C00083,C00162,C00233,C00487,C00558,C01205,C02571,C05263,C06104,C16255"
#'   )
#' )
#' human_gmt <- prepare_gmt_data(human_pathways)
#' str(human_gmt)
#'
#' # Example 5: Adjusting minimum pathway size
#' # For smaller datasets, you might want to reduce the minimum size
#' small_gmt <- prepare_gmt_data(
#'   PathwayVsMetabolites, 
#'   min_pathway_size = 3
#' )
#' str(small_gmt)
#'
#' @export
prepare_gmt_data <- function(PathwayVsMetabolites, 
                             metabolites_in_data = NULL, 
                             min_pathway_size = 5) {
    # Input validation
    if (!is.data.frame(PathwayVsMetabolites) || 
        !all(c("Pathway", "Metabolites") %in% colnames(PathwayVsMetabolites))) {
        stop("PathwayVsMetabolites must be a data frame with 'Pathway' and 'Metabolites' columns")
    }
    
    message("Preparing GMT data from ", nrow(PathwayVsMetabolites), " pathways")
    
    pathway_list <- list()
    pathway_stats <- data.frame(
        pathway = character(),
        original_size = integer(),
        filtered_size = integer(),
        stringsAsFactors = FALSE
    )
    
    # Use seq_len instead of 1:nrow()
    for (i in seq_len(nrow(PathwayVsMetabolites))) {
        pathway_name <- as.character(PathwayVsMetabolites$Pathway[i])
        metabolites_str <- as.character(PathwayVsMetabolites$Metabolites[i])
        
        # Parse metabolites - handle various separators (, | ;) and clean up
        metabolites <- unlist(strsplit(metabolites_str, "[,|;]"))
        metabolites <- trimws(metabolites)
        metabolites <- metabolites[metabolites != "" & !is.na(metabolites)]
        
        original_size <- length(metabolites)
        
        # Filter to metabolites in data if provided
        if (!is.null(metabolites_in_data)) {
            metabolites <- metabolites[metabolites %in% metabolites_in_data]
        }
        
        filtered_size <- length(metabolites)
        
        # Only include pathways with sufficient metabolites
        if (filtered_size >= min_pathway_size) {
            pathway_list[[pathway_name]] <- unique(metabolites)
        }
        
        # Store statistics
        pathway_stats <- rbind(pathway_stats, data.frame(
            pathway = pathway_name,
            original_size = original_size,
            filtered_size = filtered_size,
            stringsAsFactors = FALSE
        ))
    }
    
    # Print summary statistics
    total_pathways <- nrow(pathway_stats)
    pathways_kept <- length(pathway_list)
    pathways_removed <- total_pathways - pathways_kept
    
    message("GMT preparation summary:")
    message("  Total pathways processed: ", total_pathways)
    message("  Pathways kept (size >= ", min_pathway_size, "): ", pathways_kept)
    message("  Pathways removed (size < ", min_pathway_size, "): ", pathways_removed)
    
    if (pathways_kept > 0) {
        # Use vapply instead of sapply
        pathway_sizes <- vapply(pathway_list, length, FUN.VALUE = integer(1))
        message("  Pathway size distribution:")
        message("    Min: ", min(pathway_sizes))
        message("    Median: ", stats::median(pathway_sizes))
        message("    Max: ", max(pathway_sizes))
        message("    Mean: ", round(mean(pathway_sizes), 1))
    }
    
    return(pathway_list)
}