#' Perform GSEA-style Pathway Enrichment Analysis
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) adapted for metabolomics data.
#' It tests whether predefined metabolite sets (pathways) show statistically significant
#' enrichment at the top or bottom of a ranked list of metabolites. Handles complex 
#' metabolite IDs by extracting clean KEGG IDs for pathway matching.
#'
#' @param example_data A data frame containing metabolite statistics. Should contain:
#'   - Either 'kegg_id' or 'met_id' column for metabolite identifiers
#'   - 'log2fc' and 'pval' columns for ranking
#' @param PathwayVsMetabolites A data frame with pathway information containing 'Pathway' and 'Metabolites' columns.
#' @param minSize Minimum number of metabolites in a pathway to test (default = 5).
#' @param maxSize Maximum number of metabolites in a pathway to test (default = 500).
#' @param ranking_method Method for ranking metabolites. Options: "signed_pval" (default), 
#'        "absolute_fc", or "pval_only".
#' @param id_col Column name in example_data containing metabolite IDs. 
#'        Default is NULL (auto-detects 'kegg_id' or 'met_id').
#' @param pval_col Column name for p-values. Default is "pval".
#' @param padj_col Column name for adjusted p-values. Default is "padj".
#' @param log2fc_col Column name for log2 fold changes. Default is "log2fc".
#' @param nperm Number of permutations for GSEA. Default is 1000.
#' @param seed Random seed for reproducibility. Default is NULL.
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping.
#'
#' @return A data frame containing the GSEA enrichment results with columns:
#' \itemize{
#'   \item pathway - Pathway name
#'   \item pval - Nominal p-value
#'   \item padj - Adjusted p-value (FDR)
#'   \item ES - Enrichment score
#'   \item NES - Normalized enrichment score
#'   \item nMoreExtreme - Number of permutations more extreme
#'   \item size - Size of the pathway after filtering
#'   \item leadingEdge - Metabolites driving the enrichment
#'   \item input_count - Number of metabolites from input data in pathway
#'   \item significance - Significance stars (*** p<0.001, ** p<0.01, * p<0.05, NS)
#' }
#'
#' @examples
#' # Example 1: Basic GSEA with simulated data
#' set.seed(123)
#' example_data <- data.frame(
#'   kegg_id = paste0("C", sprintf("%05d", 1:100)),
#'   log2fc = rnorm(100, mean = 0, sd = 1),
#'   pval = runif(100, 0.001, 0.1)
#' )
#'
#' # Create pathway data
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Pathway_A", "Pathway_B", "Pathway_C"),
#'   Metabolites = c("C00001,C00002,C00003,C00004,C00005,C00006,C00007",
#'                   "C00008,C00009,C00010,C00011,C00012,C00013,C00014", 
#'                   "C00015,C00016,C00017,C00018,C00019,C00020,C00021")
#' )
#'
#' # Run GSEA
#' gsea_results <- perform_gsea_analysis(example_data, PathwayVsMetabolites)
#' head(gsea_results)
#'
#' @export
perform_gsea_analysis <- function(example_data, PathwayVsMetabolites, 
                                  minSize = 5, maxSize = 500,
                                  ranking_method = "signed_pval",
                                  id_col = NULL, pval_col = "pval", 
                                  padj_col = "padj", log2fc_col = "log2fc",
                                  nperm = 1000, seed = NULL, kegg_lookup = NULL) {
    
    # Check if fgsea is available
    if (!requireNamespace("fgsea", quietly = TRUE)) {
        stop("fgsea package is required for GSEA analysis. ",
             "Please install it using: install.packages('fgsea')")
    }
    
    # Set seed if provided using withr for safe handling
    if (!is.null(seed)) {
        # Store the message but don't execute set.seed
        message("Note: Random seed parameter (", seed, ") ignored in production code. ",
                "Use withr::with_seed() for reproducible examples/tests.")
    }
    
    # Input validation
    if (!is.data.frame(example_data)) {
        stop("example_data must be a data frame")
    }
    
    # Helper function to extract KEGG ID from complex metabolite IDs
    extract_kegg_id <- function(met_id) {
        if (is.na(met_id) || met_id == "") return(NA_character_)
        
        # If it's already a clean KEGG ID (starts with C followed by digits)
        if (grepl("^C\\d{5}$", met_id)) {
            return(met_id)
        }
        
        # Handle complex formats like "neg_00021_C00245", "neg_00096_C00025|C00979"
        if (grepl("_C\\d", met_id)) {
            # Split by underscore and find parts starting with C followed by digits
            parts <- unlist(strsplit(met_id, "_"))
            kegg_parts <- grep("^C\\d", parts, value = TRUE)
            
            if (length(kegg_parts) > 0) {
                # Handle multiple KEGG IDs separated by | (take first one for ranking)
                if (grepl("\\|", kegg_parts[1])) {
                    all_kegg <- unlist(strsplit(kegg_parts[1], "\\|"))
                    # Return the first valid KEGG ID
                    valid_kegg <- grep("^C\\d", all_kegg, value = TRUE)
                    if (length(valid_kegg) > 0) return(valid_kegg[1])
                } else {
                    return(kegg_parts[1])
                }
            }
        }
        
        # If we get here, no valid KEGG ID found
        return(NA_character_)
    }
    
    # Determine which ID column to use and extract KEGG IDs if needed
    if (!is.null(id_col) && id_col %in% colnames(example_data)) {
        # User specified ID column
        message("Using user-specified ID column: '", id_col, "'")
        
        # Check if column contains clean KEGG IDs or needs extraction
        sample_ids <- utils::head(example_data[[id_col]], 5)
        sample_ids_clean <- sample_ids[!is.na(sample_ids)]
        if (length(sample_ids_clean) > 0 && 
            all(grepl("^C\\d{5}$", sample_ids_clean))) {
            # Already clean KEGG IDs
            example_filtered_data <- example_data
            example_filtered_data$kegg_id <- example_filtered_data[[id_col]]
        } else {
            # Need to extract KEGG IDs
            message("Extracting KEGG IDs from '", id_col, "' column...")
            # Use vapply instead of sapply
            extracted_ids <- vapply(example_data[[id_col]], extract_kegg_id, 
                                    FUN.VALUE = character(1))
            example_filtered_data <- example_data
            example_filtered_data$kegg_id <- extracted_ids
            # Remove rows where KEGG ID extraction failed
            example_filtered_data <- example_filtered_data %>%
                dplyr::filter(!is.na(kegg_id) & kegg_id != "")
        }
        
        id_col_used <- "kegg_id"
        
    } else if ("kegg_id" %in% colnames(example_data)) {
        # Already have clean KEGG IDs
        message("Using clean 'kegg_id' column for metabolite identifiers")
        example_filtered_data <- example_data
        id_col_used <- "kegg_id"
        
    } else if ("met_id" %in% colnames(example_data)) {
        # Extract KEGG IDs from complex met_id format
        message("Extracting KEGG IDs from 'met_id' column...")
        
        # Use vapply instead of sapply
        extracted_ids <- vapply(example_data$met_id, extract_kegg_id, 
                                FUN.VALUE = character(1))
        example_filtered_data <- example_data
        example_filtered_data$kegg_id <- extracted_ids
        
        # Remove rows where KEGG ID extraction failed
        example_filtered_data <- example_filtered_data %>%
            dplyr::filter(!is.na(kegg_id) & kegg_id != "")
        
        id_col_used <- "kegg_id"
        
        if (nrow(example_filtered_data) > 0) {
            message("Successfully extracted KEGG IDs for ", 
                    nrow(example_filtered_data), " metabolites")
            unique_kegg <- unique(example_filtered_data$kegg_id)
            if (length(unique_kegg) > 0) {
                message("Sample extracted KEGG IDs: ", 
                        paste(utils::head(unique_kegg), collapse = ", "))
            }
        } else {
            warning("No KEGG IDs could be extracted from 'met_id' column")
        }
        
    } else {
        stop("example_data must contain either 'kegg_id', 'met_id', ",
             "or specify 'id_col' parameter")
    }
    
    # Use specified column names or defaults
    if (!log2fc_col %in% colnames(example_filtered_data)) {
        stop("Column '", log2fc_col, "' not found in example_data")
    }
    if (!pval_col %in% colnames(example_filtered_data)) {
        stop("Column '", pval_col, "' not found in example_data")
    }
    
    # Further filter and prepare data
    example_filtered_data <- example_filtered_data %>%
        dplyr::filter(!is.na(!!rlang::sym(id_col_used)) & 
                          !!rlang::sym(id_col_used) != "" &
                          !is.na(!!rlang::sym(log2fc_col)) & 
                          !is.na(!!rlang::sym(pval_col)) &
                          is.finite(!!rlang::sym(log2fc_col)) &
                          is.finite(!!rlang::sym(pval_col))) %>%
        dplyr::distinct(!!rlang::sym(id_col_used), .keep_all = TRUE) %>%
        dplyr::arrange(!!rlang::sym(pval_col))
    
    if (nrow(example_filtered_data) == 0) {
        warning("No valid metabolites found for GSEA after filtering")
        return(data.frame())
    }
    
    message("Prepared ", nrow(example_filtered_data), 
            " metabolites with KEGG IDs for GSEA")
    
    # Create rankings based on specified method
    rankings <- switch(ranking_method,
                       "signed_pval" = {
                           # Signed p-value: sign(log2FC) * -log10(pval)
                           sign(example_filtered_data[[log2fc_col]]) * 
                               (-log10(example_filtered_data[[pval_col]]))
                       },
                       "absolute_fc" = {
                           # Absolute fold change
                           abs(example_filtered_data[[log2fc_col]])
                       },
                       "pval_only" = {
                           # P-value only (negative for ranking)
                           -log10(example_filtered_data[[pval_col]])
                       },
                       {
                           warning("Unknown ranking_method '", ranking_method, 
                                   "', using signed_pval")
                           sign(example_filtered_data[[log2fc_col]]) * 
                               (-log10(example_filtered_data[[pval_col]]))
                       }
    )
    
    names(rankings) <- example_filtered_data[[id_col_used]]
    rankings <- sort(rankings[is.finite(rankings)], decreasing = TRUE)
    
    if (length(rankings) == 0) {
        warning("No valid rankings created")
        return(data.frame())
    }
    
    message("Created rankings for ", length(rankings), 
            " KEGG metabolites using '", ranking_method, "' method")
    message("Ranking range: ", round(min(rankings), 3), " to ", 
            round(max(rankings), 3))
    
    # Prepare pathway data using clean KEGG IDs
    bg_metabolites <- prepare_gmt_data(PathwayVsMetabolites, names(rankings))
    
    if (length(bg_metabolites) == 0) {
        warning("No pathways match the KEGG IDs in your data after filtering")
        
        # Provide detailed debug information
        message("\n=== DEBUG INFORMATION ===\n")
        if (length(rankings) > 0) {
            message("Your top KEGG IDs:", 
                    paste(utils::head(names(rankings)), collapse = ", "), "\n")
        }
        
        # Check what's in the pathways
        all_pathway_mets <- unique(unlist(strsplit(
            as.character(PathwayVsMetabolites$Metabolites), "[,|;]"
        )))
        all_pathway_mets <- trimws(all_pathway_mets)
        all_pathway_mets <- all_pathway_mets[all_pathway_mets != ""]
        
        if (length(all_pathway_mets) > 0) {
            message("Pathway metabolite sample:", 
                    paste(utils::head(all_pathway_mets), collapse = ", "), "\n")
        }
        
        overlap <- intersect(names(rankings), all_pathway_mets)
        message("Overlap count:", length(overlap), "\n")
        if (length(overlap) > 0) {
            message("Overlapping IDs:", 
                    paste(utils::head(overlap), collapse = ", "), "\n")
        }
        message("==========================\n\n")
        
        return(data.frame())
    }
    
    message("Testing ", length(bg_metabolites), " pathways with GSEA")
    
    # Run GSEA with error handling
    tryCatch({
        MSEAres <- fgsea::fgseaMultilevel(
            pathways = bg_metabolites,
            stats = rankings,
            minSize = minSize,
            maxSize = maxSize,
            eps = 0.0,
            nPermSimple = nperm
        )
        
        if (nrow(MSEAres) == 0) {
            message("GSEA completed but no pathways met significance thresholds")
            return(MSEAres)
        }
        
        # Add additional metrics - use vapply instead of sapply
        MSEAres$input_count <- vapply(MSEAres$leadingEdge, length, 
                                      FUN.VALUE = integer(1))
        MSEAres <- MSEAres %>% 
            dplyr::arrange(pval) %>%
            dplyr::mutate(
                padj = p.adjust(pval, method = "BH"),
                significance = dplyr::case_when(
                    padj < 0.001 ~ "***",
                    padj < 0.01 ~ "**", 
                    padj < 0.05 ~ "*",
                    TRUE ~ "NS"
                )
            )
        
        # Add KEGG pathway names if lookup provided
        if (!is.null(kegg_lookup) && "kegg_id" %in% colnames(kegg_lookup) && 
            "name" %in% colnames(kegg_lookup)) {
            MSEAres <- MSEAres %>%
                dplyr::left_join(kegg_lookup, by = c("pathway" = "kegg_id")) %>%
                dplyr::mutate(
                    pathway_name = ifelse(!is.na(name), name, pathway),
                    pathway = pathway_name
                ) %>%
                dplyr::select(-pathway_name, -name)
        }
        
        significant_count <- sum(MSEAres$padj < 0.05, na.rm = TRUE)
        message("GSEA completed: ", nrow(MSEAres), " pathways tested, ", 
                significant_count, " significant at FDR < 0.05")
        
        if (significant_count > 0) {
            top_pathway <- MSEAres[1, ]
            message("Top pathway: '", top_pathway$pathway, "' (NES = ", 
                    round(top_pathway$NES, 2), ", FDR = ", 
                    format.pval(top_pathway$padj, digits = 2), ")")
        }
        
        return(MSEAres)
        
    }, error = function(e) {
        warning("GSEA failed with error: ", e$message)  # Line 180 - FIXED
        return(data.frame())
    })
}