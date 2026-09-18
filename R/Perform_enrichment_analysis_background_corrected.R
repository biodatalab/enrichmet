#' Extract measured KEGG metabolite IDs from a summary statistics table
#'
#' Turns a per-feature summary-statistics table (typical untargeted LC-MS
#' differential analysis output) into a vector of KEGG compound IDs. This is
#' the set of metabolites that were measurable/testable on the platform and is
#' the appropriate universe for background-corrected enrichment, as opposed to
#' every compound annotated in the pathway database.
#'
#' @param summary_stats A data frame with a \code{met_id} column whose values
#'   embed KEGG IDs (pattern \code{"C"} + 5 digits), e.g. \code{"neg_00021_C00245"},
#'   or ambiguous annotations such as \code{"neg_00096_C00025|C00979"} (both
#'   candidate IDs are kept).
#'
#' @return A character vector of unique KEGG IDs found in \code{summary_stats}.
#'
#' @seealso \code{\link{perform_enrichment_analysis_bg}}
#'
#' @examples
#' 
#' example_path <- get_cached_file(
#'     "https://zenodo.org/api/records/17819145/files/summary_stat.csv/content"
#' )
#' example_data <- read.csv(example_path, stringsAsFactors = FALSE)
#' measured <- extract_measured_kegg_ids(example_data)
#' length(measured)
#' head(measured)
#' 
#'
#' @export
extract_measured_kegg_ids <- function(summary_stats) {
    if (!is.data.frame(summary_stats) || !"met_id" %in% colnames(summary_stats)) {
        stop("summary_stats must be a data frame with a 'met_id' column")
    }
    
    met_ids <- unique(summary_stats$met_id[!is.na(summary_stats$met_id) &
                                               summary_stats$met_id != ""])
    
    kegg_id_list <- regmatches(met_ids, gregexpr("C[0-9]{5}", met_ids))
    kegg_ids <- unique(unlist(kegg_id_list))
    kegg_ids <- kegg_ids[nzchar(kegg_ids)]
    
    if (length(kegg_ids) == 0) {
        stop("No KEGG IDs (pattern 'C' + 5 digits) could be parsed from summary_stats$met_id")
    }
    
    kegg_ids
}

#' Perform pathway enrichment with optional background correction
#'
#' Fisher's exact test for pathway over-representation, optionally restricted
#' to a measured-metabolite universe (background correction).
#'
#' @param inputMetabolites Character vector of query metabolite IDs (e.g. KEGG).
#' @param PathwayVsMetabolites Data frame with columns \code{Pathway} and
#'   \code{Metabolites} (comma-separated IDs).
#' @param top_n Integer; number of top pathways to return (default 100).
#' @param p_value_cutoff Numeric; filter on raw \eqn{p}-value (default 1).
#' @param backgroundMetabolites Optional character vector defining the
#'   statistical universe. If \code{NULL}, all metabolites in
#'   \code{PathwayVsMetabolites} are used. If supplied (e.g. from
#'   \code{\link{extract_measured_kegg_ids}}), the universe and each pathway's
#'   membership are intersected with this set before testing.
#'
#' @return A data frame of enrichment results (p-values, impact, coverage, etc.).
#'
#' @details
#' When \code{backgroundMetabolites} is provided, enrichment is tested only
#' among metabolites that were measured, which avoids inflating the background
#' with database compounds that were never observed.
#'
#' @importFrom dplyr filter arrange mutate
#' @importFrom tidyr unnest
#' @importFrom stats fisher.test p.adjust
#' @importFrom utils head
#'
#' @seealso \code{\link{extract_measured_kegg_ids}},
#'   \code{\link{fetch_kegg_pathway_metabolites}}
#'
#' @examples
#' # Always-runnable minimal example (no network)
#' pw <- data.frame(
#'     Pathway = c("Glycolysis", "TCA cycle"),
#'     Metabolites = c("C00031,C00022,C00074", "C00036,C00042,C00026"),
#'     stringsAsFactors = FALSE
#' )
#' query <- c("C00031", "C00022", "C00036")
#' enr <- perform_enrichment_analysis_bg(
#'     inputMetabolites = query,
#'     PathwayVsMetabolites = pw,
#'     top_n = 5,
#'     p_value_cutoff = 1
#' )
#' head(enr)
#'
#'
#' PathwayVsMetabolites <- fetch_kegg_pathway_metabolites(organism = "hsa")
#'
#' example_path <- get_cached_file(
#'     "https://zenodo.org/api/records/17819145/files/summary_stat.csv/content"
#' )
#' example_data <- read.csv(example_path, stringsAsFactors = FALSE)
#' measured <- extract_measured_kegg_ids(example_data)
#'
#' metabolomics_path <- get_cached_file(
#'     "https://zenodo.org/api/records/17819145/files/example_data.csv/content"
#' )
#' metabolomics_mat <- read.csv(
#'     metabolomics_path, row.names = 1, check.names = FALSE
#' )
#' da_out <- run_de(
#'     metabolomics_mat, "TK-CMV", "K-CMV",
#'     fc_threshold = 1, pval_threshold = 0.05
#' )
#' query <- da_out$kegg_ready$kegg_id[
#'     da_out$kegg_ready$Significant != "Not significant"
#' ]
#' query <- query[!is.na(query) & query != ""]
#'
#' enr <- perform_enrichment_analysis_bg(
#'     inputMetabolites = query,
#'     PathwayVsMetabolites = PathwayVsMetabolites,
#'     top_n = 10,
#'     p_value_cutoff = 1
#' )
#' head(enr[, c("Pathway", "P_value", "Adjusted_P_value")], 3)
#'
#' enr_bg <- perform_enrichment_analysis_bg(
#'     inputMetabolites = query,
#'     PathwayVsMetabolites = PathwayVsMetabolites,
#'     top_n = 10,
#'     p_value_cutoff = 1,
#'     backgroundMetabolites = measured
#' )
#' head(enr_bg[, c("Pathway", "P_value", "Adjusted_P_value")], 3)
#' 
#'
#' @export
perform_enrichment_analysis_bg <- function(inputMetabolites,
                                           PathwayVsMetabolites,
                                           top_n = 100,
                                           p_value_cutoff = 1,
                                           backgroundMetabolites = NULL) {
    if (!requireNamespace("qvalue", quietly = TRUE)) {
        stop("Package 'qvalue' is required. Install with BiocManager::install(\"qvalue\").")
    }
    
    if (!is.character(inputMetabolites) || length(inputMetabolites) == 0) {
        stop("inputMetabolites must be a non-empty character vector")
    }
    
    if (!is.data.frame(PathwayVsMetabolites) ||
        !all(c("Pathway", "Metabolites") %in% colnames(PathwayVsMetabolites))) {
        stop("PathwayVsMetabolites must be a data frame with 'Pathway' and 'Metabolites' columns")
    }
    
    if (!is.null(backgroundMetabolites) &&
        (!is.character(backgroundMetabolites) || length(backgroundMetabolites) == 0)) {
        stop("backgroundMetabolites must be NULL or a non-empty character vector of metabolite IDs")
    }
    
    PathwayVsMetabolites_clean <- PathwayVsMetabolites %>%
        dplyr::filter(!is.na(Pathway),
                      !is.na(Metabolites),
                      Pathway != "",
                      Metabolites != "",
                      !grepl("^\\s*$", Metabolites))
    
    if (nrow(PathwayVsMetabolites_clean) == 0) {
        stop("PathwayVsMetabolites contains no valid data after cleaning NA/empty values")
    }
    
    data <- PathwayVsMetabolites_clean %>%
        dplyr::mutate(Metabolites = strsplit(as.character(Metabolites), ",")) %>%
        tidyr::unnest(Metabolites) %>%
        dplyr::mutate(Metabolites = trimws(Metabolites)) %>%
        dplyr::filter(Metabolites != "", !is.na(Metabolites))
    
    if (nrow(data) == 0) {
        stop("No valid metabolite-pathway relationships found after cleaning")
    }
    
    allMetabolitesSet <- unique(data$Metabolites)
    
    if (!is.null(backgroundMetabolites)) {
        original_universe_size <- length(allMetabolitesSet)
        allMetabolitesSet <- intersect(allMetabolitesSet, unique(backgroundMetabolites))
        
        if (length(allMetabolitesSet) == 0) {
            stop("No overlap between PathwayVsMetabolites and backgroundMetabolites; check ID formats (e.g. KEGG IDs vs. other identifiers)")
        }
        
        message(sprintf(
            "Background correction applied: universe restricted from %d (all pathway-annotated metabolites) to %d (measured metabolites also present in the pathway database)",
            original_universe_size, length(allMetabolitesSet)
        ))
    }
    
    inputMetabolites <- inputMetabolites[inputMetabolites %in% allMetabolitesSet]
    if (length(inputMetabolites) == 0) {
        stop("No input metabolites found in the pathway database (after applying backgroundMetabolites, if provided)")
    }
    
    metabolite_centrality <- tryCatch({
        calculate_metabolite_centrality(PathwayVsMetabolites_clean)
    }, error = function(e) {
        warning("Centrality calculation failed: ", e$message,
                "\nProceeding without centrality-based impact scores")
        return(data.frame(Metabolite = character(0),
                          RBC_Metabolite = numeric(0)))
    })
    
    pathway_list <- unique(PathwayVsMetabolites_clean$Pathway)
    results <- vector("list", length(pathway_list))
    
    for (i in seq_along(pathway_list)) {
        pathway <- pathway_list[i]
        
        pathway_rows <- PathwayVsMetabolites_clean %>%
            dplyr::filter(Pathway == pathway)
        
        pathwayMetabolites <- unique(unlist(strsplit(as.character(pathway_rows$Metabolites), ",")))
        pathwayMetabolites <- trimws(pathwayMetabolites)
        pathwayMetabolites <- pathwayMetabolites[pathwayMetabolites != ""]
        
        if (!is.null(backgroundMetabolites)) {
            pathwayMetabolites <- intersect(pathwayMetabolites, allMetabolitesSet)
        }
        
        if (length(pathwayMetabolites) == 0) {
            next
        }
        
        a <- sum(inputMetabolites %in% pathwayMetabolites)
        b <- length(inputMetabolites) - a
        c <- length(pathwayMetabolites) - a
        d <- length(allMetabolitesSet) - a - b - c
        
        if (a == 0) next
        
        if (any(c(a, b, c, d) < 0) || a + b + c + d != length(allMetabolitesSet)) {
            next
        }
        
        contingency_table <- matrix(c(a, c, b, d), nrow = 2, byrow = TRUE)
        
        fisher_test_result <- tryCatch({
            fisher.test(contingency_table, alternative = "greater")
        }, error = function(e) {
            list(p.value = NA)
        })
        
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
    
    is_not_null <- vapply(results, function(x) !is.null(x), FUN.VALUE = logical(1))
    results <- results[is_not_null]
    
    if (length(results) == 0) {
        warning("No significant pathway enrichments found")
        return(data.frame(
            Pathway = character(0), P_value = numeric(0), Log_P_value = numeric(0),
            Impact = numeric(0), Coverage = numeric(0), Count = integer(0),
            Pathway_Size = integer(0), Input_Size = integer(0),
            Adjusted_P_value = numeric(0), Q_value = numeric(0),
            Enrichment_Ratio = numeric(0), Metabolite_List = character(0),
            stringsAsFactors = FALSE
        ))
    }
    
    results_combined <- do.call(rbind, results)
    results_combined$Adjusted_P_value <- p.adjust(results_combined$P_value, method = "BH")
    
    clean_p_values <- results_combined$P_value
    valid_idx <- which(clean_p_values > 0 & clean_p_values < 1)
    
    if (length(valid_idx) >= 10) {
        qobj <- tryCatch(qvalue::qvalue(clean_p_values[valid_idx]), error = function(e) NULL)
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
    
    results_combined$Enrichment_Ratio <- with(results_combined, {
        pathway_prop <- Pathway_Size / length(allMetabolitesSet)
        ifelse(pathway_prop == 0, NA_real_, (Count / Input_Size) / pathway_prop)
    })
    
    results_filtered <- results_combined %>%
        dplyr::filter(P_value <= p_value_cutoff)
    
    if (nrow(results_filtered) == 0) {
        warning(sprintf("No pathways pass the p-value cutoff of %f", p_value_cutoff))
    }
    
    if (!is.null(top_n) && nrow(results_filtered) > 0) {
        if (top_n > 0 && nrow(results_filtered) > top_n) {
            results_topN <- results_filtered %>% dplyr::arrange(P_value) %>% head(top_n)
        } else {
            results_topN <- results_filtered %>% dplyr::arrange(P_value)
        }
    } else {
        results_topN <- results_filtered
    }
    
    if (nrow(results_topN) > 0) {
        results_topN$Metabolite_List <- vapply(results_topN$Pathway, function(p) {
            pathway_rows <- PathwayVsMetabolites_clean %>% dplyr::filter(Pathway == p)
            pathway_mets <- unique(unlist(strsplit(as.character(pathway_rows$Metabolites), ",")))
            pathway_mets <- trimws(pathway_mets)
            overlapping_mets <- intersect(pathway_mets, inputMetabolites)
            paste(overlapping_mets, collapse = ",")
        }, FUN.VALUE = character(1))
    }
    
    message(sprintf(
        "Enrichment analysis completed: %d pathways tested, %d pathways passed filtering (p <= %f)",
        nrow(results_combined), nrow(results_topN), p_value_cutoff
    ))
    
    return(results_topN)
}