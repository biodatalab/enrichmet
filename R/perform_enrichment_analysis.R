#' Perform Pathway Enrichment Analysis
#'
#' Uses Fisher's exact test to identify enriched pathways for a set of metabolites.
#'
#' @param inputMetabolites A character vector of metabolite IDs for which pathway 
#'        enrichment analysis is to be performed.
#' @param PathwayVsMetabolites A data frame containing pathways and their associated 
#'        metabolites. Must include columns 'Pathway' and 'Metabolites'.
#' @param top_n An integer specifying the number of top significant pathways to return (default is 100).
#' @param p_value_cutoff A numeric value for adjusting the p-value threshold for 
#'        filtering significant pathways (default is 1, no filtering).
#'
#' @return A data frame of pathway enrichment results including p-values, impact scores, 
#'         and coverage statistics.
#'
#' @examples
#' # Create example data
#' inputMetabolites <- c("M1", "M2", "M3")
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Pathway1", "Pathway2"),
#'   Metabolites = c("M1,M2,M3", "M2,M4,M5")
#' )
#'
#' # Perform enrichment analysis
#' results <- perform_enrichment_analysis(inputMetabolites, PathwayVsMetabolites)
#' head(results)
#'# See ?enrichmet for complete examples
#' @importFrom dplyr filter arrange mutate
#' @importFrom tidyr unnest
#' @importFrom stats fisher.test p.adjust
#' @importFrom utils head
#' @export
perform_enrichment_analysis <- function(inputMetabolites, PathwayVsMetabolites, 
                                        top_n = 100, p_value_cutoff = 1) {
    # Load required package
    if (!requireNamespace("qvalue", quietly = TRUE)) {
        stop("Please install the 'qvalue' package: install.packages('qvalue')")
    }
    
    # Input validation
    if (!is.character(inputMetabolites) || length(inputMetabolites) == 0) {
        stop("inputMetabolites must be a non-empty character vector")
    }
    
    if (!is.data.frame(PathwayVsMetabolites) || 
        !all(c("Pathway", "Metabolites") %in% colnames(PathwayVsMetabolites))) {
        stop("PathwayVsMetabolites must be a data frame with 'Pathway' and 'Metabolites' columns")
    }
    
    # Expand pathway-metabolite relationships
    data <- PathwayVsMetabolites %>%
        dplyr::mutate(Metabolites = strsplit(Metabolites, ",")) %>%
        tidyr::unnest(Metabolites)
    
    allMetabolitesSet <- unique(data$Metabolites)
    
    # Centrality calculation
    metabolite_centrality <- calculate_metabolite_centrality(PathwayVsMetabolites)
    
    results <- list()
    
    for (i in seq_len(nrow(PathwayVsMetabolites))) {
        row <- PathwayVsMetabolites[i, ]
        pathway <- row$Pathway
        pathwayMetabolites <- unlist(strsplit(row$Metabolites, ","))
        matchedMet <- intersect(pathwayMetabolites, inputMetabolites)
        
        # Skip Fisher's test if no matches
        if (length(matchedMet) == 0) {
            results[[i]] <- data.frame(
                Pathway = pathway,
                P_value = 1.0,
                Log_P_value = 0,
                Impact = 0,
                Coverage = 0,
                Count = 0
            )
            next
        }
        
        # Fisher's exact test
        a <- length(matchedMet)
        b <- length(setdiff(inputMetabolites, pathwayMetabolites))
        c <- length(setdiff(pathwayMetabolites, inputMetabolites))
        d <- length(setdiff(allMetabolitesSet, union(inputMetabolites, pathwayMetabolites)))
        
        contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
        fisher_test_result <- fisher.test(contingency_table, alternative = "two.sided")
        
        # Impact & Coverage
        matched_centrality <- metabolite_centrality %>%
            dplyr::filter(Metabolite %in% matchedMet)
        
        all_centrality <- metabolite_centrality %>%
            dplyr::filter(Metabolite %in% pathwayMetabolites)
        
        impact <- ifelse(
            nrow(all_centrality) > 0,
            sum(matched_centrality$RBC_Metabolite) / sum(all_centrality$RBC_Metabolite),
            0
        )
        coverage <- length(matchedMet) / length(pathwayMetabolites)
        count <- length(matchedMet)   
        
        # Save results
        results[[i]] <- data.frame(
            Pathway = pathway,
            P_value = fisher_test_result$p.value,
            Log_P_value = -log10(fisher_test_result$p.value),
            Impact = impact,
            Coverage = coverage,
            Count = count
        )
    }
    
    # Combine and adjust
    results_combined <- do.call(rbind, lapply(results, as.data.frame))
    results_combined$Adjusted_P_value <- p.adjust(results_combined$P_value, method = "BH")
    
    # Clean p-values for qvalue calculation
    clean_p_values <- results_combined$P_value
    
    # Handle edge cases for qvalue
    if (all(is.na(clean_p_values))) {
        results_combined$Q_value <- NA
    } else {
        # Remove NA/NaN values and ensure p-values are in [0, 1]
        clean_p_values <- clean_p_values[!is.na(clean_p_values) & !is.nan(clean_p_values)]
        clean_p_values <- pmax(pmin(clean_p_values, 1), 0)  # Clamp to [0, 1]
        
        # Only calculate q-values if there are valid p-values
        if (length(clean_p_values) > 0 && any(clean_p_values > 0 & clean_p_values < 1)) {
            qobj <- qvalue::qvalue(clean_p_values)
            # Create full Q_value vector with NA for invalid entries
            full_qvalues <- rep(NA, nrow(results_combined))
            valid_indices <- which(!is.na(results_combined$P_value) & 
                                       !is.nan(results_combined$P_value) &
                                       results_combined$P_value >= 0 & 
                                       results_combined$P_value <= 1)
            full_qvalues[valid_indices] <- qobj$qvalues
            results_combined$Q_value <- full_qvalues
        } else {
            results_combined$Q_value <- NA
        }
    }
    
    # Filter & sort
    significant_results_df <- results_combined %>%
        dplyr::filter(Adjusted_P_value < p_value_cutoff) %>%
        dplyr::arrange(desc(Log_P_value))
    
    if (!is.null(top_n)) {
        significant_results_df <- head(significant_results_df, top_n)
    }
    
    return(significant_results_df)
}