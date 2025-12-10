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
#' # Create meaningful pathway data with real KEGG IDs
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Glycolysis / Gluconeogenesis",
#'               "Citrate cycle (TCA cycle)",
#'               "Pentose phosphate pathway",
#'               "Amino sugar and nucleotide sugar metabolism",
#'               "Pyruvate metabolism"),
#'   Metabolites = c("C00031,C00022,C00197,C00221,C00631,C01172",
#'                   "C00022,C00036,C00024,C00158,C00149,C00311",
#'                   "C00117,C00257,C00121,C00085,C00118",
#'                   "C00031,C00095,C00103,C00185",
#'                   "C00022,C00024,C00036,C00158")
#' )
#'
#' # Example: Input metabolites enriched in Glycolysis
#' inputMetabolites <- c("C00031", "C00022", "C00197", "C00221", "C00631", "C01172")
#'
#' # Perform enrichment analysis
#' results <- perform_enrichment_analysis(
#'   inputMetabolites = inputMetabolites,
#'   PathwayVsMetabolites = PathwayVsMetabolites,
#'   p_value_cutoff = 0.05
#' )
#' head(results)
#'
#' @importFrom dplyr filter arrange mutate
#' @importFrom tidyr unnest
#' @importFrom stats fisher.test p.adjust
#' @importFrom utils head
#' @export
perform_enrichment_analysis <- function(inputMetabolites, 
                                        PathwayVsMetabolites, 
                                        top_n = 100, 
                                        p_value_cutoff = 1) {
    # Load required packages
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
    
    # Remove any NA or empty values from PathwayVsMetabolites
    PathwayVsMetabolites_clean <- PathwayVsMetabolites %>%
        dplyr::filter(!is.na(Pathway), 
                      !is.na(Metabolites),
                      Pathway != "",
                      Metabolites != "",
                      !grepl("^\\s*$", Metabolites))  # Remove whitespace-only entries
    
    if (nrow(PathwayVsMetabolites_clean) == 0) {
        stop("PathwayVsMetabolites contains no valid data after cleaning NA/empty values")
    }
    
    # Expand pathway-metabolite relationships with proper cleaning
    data <- PathwayVsMetabolites_clean %>%
        dplyr::mutate(Metabolites = strsplit(as.character(Metabolites), ",")) %>%
        tidyr::unnest(Metabolites) %>%
        dplyr::mutate(Metabolites = trimws(Metabolites)) %>%  # Trim whitespace
        dplyr::filter(Metabolites != "", !is.na(Metabolites))  # Remove empty metabolites
    
    if (nrow(data) == 0) {
        stop("No valid metabolite-pathway relationships found after cleaning")
    }
    
    allMetabolitesSet <- unique(data$Metabolites)
    
    # Keep only input metabolites that exist in universe
    inputMetabolites <- inputMetabolites[inputMetabolites %in% allMetabolitesSet]
    if (length(inputMetabolites) == 0) {
        stop("No input metabolites found in the pathway database")
    }
    
    # Calculate centrality with error handling
    metabolite_centrality <- tryCatch({
        calculate_metabolite_centrality(PathwayVsMetabolites_clean)
    }, error = function(e) {
        warning("Centrality calculation failed: ", e$message, 
                "\nProceeding without centrality-based impact scores")
        return(data.frame(Metabolite = character(0), 
                          RBC_Metabolite = numeric(0)))
    })
    
    # Prepare pathway list for iteration
    pathway_list <- unique(PathwayVsMetabolites_clean$Pathway)
    results <- vector("list", length(pathway_list))
    
    for (i in seq_along(pathway_list)) {
        pathway <- pathway_list[i]
        
        # Get metabolites for this pathway
        pathway_rows <- PathwayVsMetabolites_clean %>%
            dplyr::filter(Pathway == pathway)
        
        pathwayMetabolites <- unique(unlist(strsplit(as.character(pathway_rows$Metabolites), ",")))
        pathwayMetabolites <- trimws(pathwayMetabolites)
        pathwayMetabolites <- pathwayMetabolites[pathwayMetabolites != ""]
        
        if (length(pathwayMetabolites) == 0) {
            next  # Skip pathways with no valid metabolites
        }
        
        # Fisher's exact test: correct contingency table
        a <- sum(inputMetabolites %in% pathwayMetabolites)  # In input AND pathway
        b <- length(inputMetabolites) - a                   # In input BUT NOT pathway
        c <- length(pathwayMetabolites) - a                # In pathway BUT NOT input
        d <- length(allMetabolitesSet) - a - b - c         # In neither
        
        # Skip pathways with no overlap
        if (a == 0) next
        
        # Ensure valid contingency table
        if (any(c(a, b, c, d) < 0) || a + b + c + d != length(allMetabolitesSet)) {
            next
        }
        
        contingency_table <- matrix(c(a, c, b, d), nrow = 2, byrow = TRUE)
        
        fisher_test_result <- tryCatch({
            fisher.test(contingency_table, alternative = "greater")
        }, error = function(e) {
            list(p.value = NA)
        })
        
        # Calculate impact if centrality data is available
        impact <- 0
        if (nrow(metabolite_centrality) > 0) {
            matched_centrality <- metabolite_centrality %>%
                dplyr::filter(Metabolite %in% intersect(inputMetabolites, pathwayMetabolites))
            
            all_pathway_centrality <- metabolite_centrality %>%
                dplyr::filter(Metabolite %in% pathwayMetabolites)
            
            if (nrow(all_pathway_centrality) > 0 && 
                sum(all_pathway_centrality$RBC_Metabolite, na.rm = TRUE) > 0) {
                impact <- sum(matched_centrality$RBC_Metabolite, na.rm = TRUE) / 
                    sum(all_pathway_centrality$RBC_Metabolite, na.rm = TRUE)
            }
        }
        
        # Calculate coverage
        coverage <- ifelse(length(pathwayMetabolites) > 0,
                           length(intersect(inputMetabolites, pathwayMetabolites)) / 
                               length(pathwayMetabolites),
                           0)
        
        count <- length(intersect(inputMetabolites, pathwayMetabolites))
        
        results[[i]] <- data.frame(
            Pathway = pathway,
            P_value = ifelse(is.na(fisher_test_result$p.value), 1, fisher_test_result$p.value),
            Log_P_value = ifelse(is.na(fisher_test_result$p.value), 0, 
                                 -log10(fisher_test_result$p.value)),
            Impact = impact,
            Coverage = coverage,
            Count = count,
            Pathway_Size = length(pathwayMetabolites),
            Input_Size = length(inputMetabolites),
            stringsAsFactors = FALSE
        )
    }
    
    # Combine results and remove NULL entries using vapply()
    is_not_null <- vapply(results, function(x) !is.null(x), FUN.VALUE = logical(1))
    results <- results[is_not_null]
    
    if (length(results) == 0) {
        warning("No significant pathway enrichments found")
        return(data.frame(
            Pathway = character(0),
            P_value = numeric(0),
            Log_P_value = numeric(0),
            Impact = numeric(0),
            Coverage = numeric(0),
            Count = integer(0),
            Pathway_Size = integer(0),
            Input_Size = integer(0),
            Adjusted_P_value = numeric(0),
            Q_value = numeric(0),
            Enrichment_Ratio = numeric(0),
            Metabolite_List = character(0),
            stringsAsFactors = FALSE
        ))
    }
    
    results_combined <- do.call(rbind, results)
    
    # BH-adjusted p-values
    results_combined$Adjusted_P_value <- p.adjust(results_combined$P_value, method = "BH")
    
    # Q-value calculation (relaxed method, ignore P=1)
    clean_p_values <- results_combined$P_value
    valid_idx <- which(clean_p_values > 0 & clean_p_values < 1)
    
    if (length(valid_idx) >= 10) {
        qobj <- tryCatch(qvalue::qvalue(clean_p_values[valid_idx]), 
                         error = function(e) NULL)
        full_qvalues <- rep(NA, length(clean_p_values))
        if (!is.null(qobj)) {
            full_qvalues[valid_idx] <- qobj$qvalues
        } else {
            full_qvalues[valid_idx] <- results_combined$Adjusted_P_value[valid_idx]
        }
        results_combined$Q_value <- full_qvalues
    } else {
        results_combined$Q_value <- results_combined$Adjusted_P_value
    }
    
    # Enrichment ratio (with vectorized ifelse instead of if)
    results_combined$Enrichment_Ratio <- with(results_combined, {
        pathway_prop <- Pathway_Size / length(allMetabolitesSet)
        # Use ifelse for vectorized conditional
        ifelse(pathway_prop == 0, 
               NA_real_,  # Return NA if pathway_prop is 0
               (Count / Input_Size) / pathway_prop)
    })
    
    # Filter by p-value cutoff
    results_filtered <- results_combined %>%
        dplyr::filter(P_value <= p_value_cutoff)
    
    if (nrow(results_filtered) == 0) {
        warning(sprintf("No pathways pass the p-value cutoff of %f", p_value_cutoff))
    }
    
    # Apply top_n
    if (!is.null(top_n) && nrow(results_filtered) > 0) {
        if (top_n > 0 && nrow(results_filtered) > top_n) {
            results_topN <- results_filtered %>%
                dplyr::arrange(P_value) %>%
                head(top_n)
        } else {
            results_topN <- results_filtered %>%
                dplyr::arrange(P_value)
        }
    } else {
        results_topN <- results_filtered
    }
    
    # Add metabolite list for each pathway using vapply()
    if (nrow(results_topN) > 0) {
        results_topN$Metabolite_List <- vapply(results_topN$Pathway, function(p) {
            pathway_rows <- PathwayVsMetabolites_clean %>%
                dplyr::filter(Pathway == p)
            pathway_mets <- unique(unlist(strsplit(as.character(pathway_rows$Metabolites), ",")))
            pathway_mets <- trimws(pathway_mets)
            overlapping_mets <- intersect(pathway_mets, inputMetabolites)
            paste(overlapping_mets, collapse = ",")
        }, FUN.VALUE = character(1))
    }
    
    message(sprintf("Enrichment analysis completed: %d pathways tested,",
                          "%d pathways passed filtering (p <= %f)"),
                    nrow(results_combined), 
                    nrow(results_topN), 
                    p_value_cutoff)
    
    return(results_topN)
}