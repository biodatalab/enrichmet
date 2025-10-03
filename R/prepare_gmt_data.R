#' Prepare GMT Data for Pathway Analysis
#'
#' Converts pathway-metabolite data to GMT format required for GSEA analysis.
#'
#' @param PathwayVsMetabolites A data frame containing pathways and their associated 
#'        metabolites. Must include columns 'Pathway' and 'Metabolites'.
#' @param metabolites_in_data Optional character vector of metabolite IDs to filter pathways.
#'
#' @return A list of pathways with their associated metabolites in GMT format.
#'
#' @examples
#' # Create example pathway data
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Pathway1", "Pathway2"),
#'   Metabolites = c("M1,M2,M3", "M2,M4,M5")
#' )
#'
#' # Convert to GMT format
#' gmt_data <- prepare_gmt_data(PathwayVsMetabolites)
#' str(gmt_data)
#'# See ?enrichmet for complete examples
#' @importFrom fgsea gmtPathways
#' @export
prepare_gmt_data <- function(PathwayVsMetabolites, metabolites_in_data = NULL) {
    # Input validation
    if (!is.data.frame(PathwayVsMetabolites) || 
        !all(c("Pathway", "Metabolites") %in% colnames(PathwayVsMetabolites))) {
        stop("PathwayVsMetabolites must be a data frame with 'Pathway' and 'Metabolites' columns")
    }
    
    # Convert to GMT format
    PathwayVsMetabolites$description <- "https://www.genome.jp/kegg/pathway.html#metabolism"
    
    convert_to_gmt <- function(pathway, description, metabolites) {
        pathway_underscore <- gsub(" ", "_", pathway)
        metabolites_vector <- unlist(strsplit(metabolites, ","))
        gmt_line <- paste(
            pathway_underscore,
            description,
            paste(metabolites_vector, collapse = "\t"),
            sep = "\t"
        )
        return(gmt_line)
    }
    
    gmt_data <- mapply(
        convert_to_gmt,
        PathwayVsMetabolites$Pathway,
        PathwayVsMetabolites$description,
        PathwayVsMetabolites$Metabolites,
        SIMPLIFY = TRUE
    )
    
    # Write GMT to temporary file
    gmt_file <- tempfile(pattern = "pathways_", fileext = ".gmt")
    writeLines(gmt_data, gmt_file)
    
    # Prepare GMT pathways
    gmt <- fgsea::gmtPathways(gmt_file)
    hidden <- unique(unlist(gmt))
    mat <- matrix(
        NA,
        dimnames = list(hidden, names(gmt)),
        nrow = length(hidden),
        ncol = length(gmt)
    )
    
    for (i in seq_len(dim(mat)[2])) {
        mat[, i] <- as.numeric(hidden %in% gmt[[i]])
    }
    
    # Filter if metabolites_in_data provided
    if (!is.null(metabolites_in_data)) {
        hidden1 <- intersect(metabolites_in_data, hidden)
        mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1, ]) > 5)]]
    }
    
    # Helper function to convert matrix to list
    matrix_to_list <- function(pws) {
        pws.l <- list()
        for (pw in colnames(pws)) {
            pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
        }
        return(pws.l)
    }
    
    final_list <- matrix_to_list(mat)
    
    # Clean up temp file
    unlink(gmt_file)
    
    return(final_list)
}
