#' Perform Metabolite Set Enrichment Analysis (GSEA)
#'
#' Implements pre-ranked GSEA using fgsea algorithm for metabolomics data.
#'
#' @param example_data A data frame containing metabolite-level data for GSEA analysis. 
#'        Should include columns "met_id", "pval", and "log2fc".
#' @param PathwayVsMetabolites A data frame containing pathways and their associated 
#'        metabolites. Must include columns 'Pathway' and 'Metabolites'.
#'
#' @return A data frame of GSEA results including normalized enrichment scores (NES) 
#'         and p-values.
#'
#' @examples
#' # Create example data
#' example_data <- data.frame(
#'   met_id = c("M1", "M2", "M3", "M4", "M5"),
#'   pval = c(0.001, 0.01, 0.05, 0.1, 0.5),
#'   log2fc = c(2.0, 1.5, -1.0, 0.5, -0.5)
#' )
#'
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Pathway1", "Pathway2"),
#'   Metabolites = c("M1,M2,M3", "M3,M4,M5")
#' )
#'
#' # Perform GSEA analysis
#' gsea_results <- perform_gsea_analysis(example_data, PathwayVsMetabolites)
#' head(gsea_results)
#'# See ?enrichmet for complete examples
#' @importFrom dplyr arrange distinct filter
#' @importFrom fgsea fgsea
#' @export
perform_gsea_analysis <- function(example_data, PathwayVsMetabolites) {
    # Input validation
    if (!is.data.frame(example_data) || 
        !all(c("met_id", "pval", "log2fc") %in% colnames(example_data))) {
        stop("example_data must be a data frame with 'met_id', 'pval', and 'log2fc' columns")
    }
    
    if (!is.data.frame(PathwayVsMetabolites) || 
        !all(c("Pathway", "Metabolites") %in% colnames(PathwayVsMetabolites))) {
        stop("PathwayVsMetabolites must be a data frame with 'Pathway' and 'Metabolites' columns")
    }
    
    # Prepare metabolite rankings
    example_filtered_data <- example_data %>%
        dplyr::arrange(pval) %>%
        dplyr::distinct(met_id, .keep_all = TRUE) %>%
        dplyr::filter(met_id != "No Metabolites found")
    
    meta <- example_filtered_data$met_id
    bg_metabolites <- prepare_gmt_data(PathwayVsMetabolites, meta)
    
    rankings <- sign(example_filtered_data$log2fc) * (-log10(as.numeric(example_filtered_data$pval)))
    names(rankings) <- example_filtered_data$met_id
    rankings <- sort(rankings, decreasing = TRUE)
    
    # Perform GSEA
    MSEAres <- fgsea::fgsea(
        pathways = bg_metabolites,
        stats = rankings,
        scoreType = 'std',
        minSize = 10,
        maxSize = 500,
        nproc = 1
    )
    
    MSEAres$input_count <- vapply(MSEAres$leadingEdge, length, integer(1))
    MSEAres <- MSEAres %>% dplyr::arrange(pval)
    
    return(MSEAres)
}
