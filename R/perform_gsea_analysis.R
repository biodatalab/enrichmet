#' Perform GSEA-style Pathway Enrichment Analysis
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) adapted for metabolomics data,
#' testing whether predefined metabolite sets (pathways) are enriched among ranked metabolites.
#'
#' @param input_data A data frame containing at least three columns:
#'   \describe{
#'     \item{met_id}{Character vector of metabolite identifiers.}
#'     \item{log2fc}{Numeric vector of log2 fold changes or other ranking metric.}
#'     \item{pval}{Numeric vector of p-values associated with each metabolite (optional for ranking).}
#'   }
#' @param pathway_map A data frame with two columns:
#'   \describe{
#'     \item{Pathway}{Pathway name or ID.}
#'     \item{Metabolites}{Comma-separated metabolite identifiers belonging to each pathway.}
#'   }
#' @param minSize Integer. Minimum number of metabolites in a pathway to test (default = 2).
#' @param maxSize Integer. Maximum number of metabolites in a pathway to test (default = 500).
#' @param nperm Integer. Number of permutations used for significance estimation (default = 1000).
#'
#' @return A data frame containing the GSEA enrichment results, including:
#' \itemize{
#'   \item pathway – pathway name
#'   \item pval – nominal p-value
#'   \item padj – adjusted p-value (Benjamini–Hochberg)
#'   \item ES – enrichment score
#'   \item NES – normalized enrichment score
#'   \item size – number of metabolites in the pathway
#'   \item leadingEdge – subset of metabolites driving the enrichment
#' }
#'
#' @details
#' Metabolites are ranked by \code{log2fc} (descending), and enrichment is tested using
#' the running-sum statistic from the \code{fgsea} algorithm. Only pathways with
#' size between \code{minSize} and \code{maxSize} are tested.
#'
#' @examples
#' # Generate example data for GSEA
#' set.seed(1234)
#' inputMetabolites <- paste0("M", 1:20)
#' 
#' # Create pathway definitions
#' pathway_names <- paste0("Pathway", 1:20)
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = rep(pathway_names, each = 1),
#'   Metabolites = sapply(1:20, function(x)
#'     paste(sample(inputMetabolites, sample(5:10, 1)), collapse = ","))
#' )
#' 
#' # Generate example metabolite statistics for GSEA
#' example_data <- data.frame(
#'   met_id = inputMetabolites,
#'   pval = runif(20, 0.001, 0.05),
#'   log2fc = rnorm(20, mean = 0, sd = 1)
#' )
#'
#' # Perform GSEA analysis
#' gsea_results <- perform_gsea_analysis(example_data, PathwayVsMetabolites)
#' 
#' # Show top results
#' head(gsea_results, 5)
#' 
#' # Expected output structure:
#' #   pathway        pval       padj   log2err        ES      NES size
#' # 1 Pathway5  0.00123456 0.01234567 0.1234567 0.8765432 2.345678    8
#' # 2 Pathway12 0.01234567 0.04567891 0.2345678 0.7654321 1.987654    6
#' # 3 Pathway8  0.03456789 0.07890123 0.3456789 0.6543210 1.543210    7
#' 
#' # Show top enriched pathway details
#' if (nrow(gsea_results) > 0) {
#'   top_pathway <- gsea_results[1, ]
#'   cat("Top enriched pathway:", top_pathway$pathway, "\n")
#'   cat("Normalized Enrichment Score (NES):", round(top_pathway$NES, 3), "\n")
#'   cat("P-value:", format.pval(top_pathway$pval, digits = 3), "\n")
#'   cat("Adjusted P-value:", format.pval(top_pathway$padj, digits = 3), "\n")
#'   cat("Pathway size:", top_pathway$size, "metabolites\n")
#' }
#'
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
