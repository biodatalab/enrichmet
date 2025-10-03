# utils.R or anywhere near the top of your R script
utils::globalVariables(
    c(
        "%>%",
        "Metabolites",
        "Metabolite",
        "desc",
        "RBC_Metabolite",
        "Adjusted_P_value",
        "Log_P_value",
        "Pathway",
        "Impact",
        "P_value",
        "pval",
        "met_id",
        "input_count",
        "NES",
        "name",
        "V",
        "centrality",
        "metabolite",
        "KEGG_ID",
        "PubChem_CID",
        "display_name",
        "STITCH_ID",
        "everything",
        "combined_score",
        "chemical1",
        "chemical2",
        "weight",
        "degree",
        "component",
        "Matched_Metabolites",
        "Total_Pathway_Metabolites",
        "Coverage",
        "Display_Name",
        "type",
        "pathway",
        "similarity",
        "experimental",
        "database",
        "textmining",
        "from",
        "to",
        "edge_alpha",
        "str_pad",
        "slice_sample",
        "expand.grid",
        "pathway_name",
        "membership_matrix",
        "heatmap_values",
        "logp_vec"
    )
)
# enrichmet: Pathway enrichment and visualization for metabolomics
# Copyright (C) 2025 Yonatan Ayalew Mekonnen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#' ENRICHMET: Comprehensive Pathway Analysis for Metabolomics Data
#'
#' A modular tool for metabolite pathway enrichment analysis that can perform multiple
#' analytical tasks including pathway enrichment, metabolite set enrichment analysis (GSEA),
#' network centrality analysis, and interaction network visualization. The function can
#' be run as a complete workflow or using individual analysis steps.
#'
#' @param inputMetabolites A character vector of metabolite IDs (KEGG IDs recommended) 
#'        for which pathway enrichment and centrality analysis are to be performed.
#' @param PathwayVsMetabolites A data frame containing pathways and their associated 
#'        metabolites. Must include columns 'Pathway' and 'Metabolites'.
#' @param example_data A data frame containing metabolite-level data for GSEA analysis. 
#'        Should include columns "met_id", "pval", and "log2fc". Required for GSEA analysis.
#' @param top_n An integer specifying the number of top pathways to include in the 
#'        pathway enrichment results (default is 100). Use NULL to return all pathways.
#' @param p_value_cutoff A numeric value for adjusting the p-value threshold for 
#'        filtering significant pathways (default is 1, no filtering).
#' @param kegg_lookup Optional data frame for KEGG ID to name mapping. Should contain
#'        columns 'kegg_id' and 'name'.
#' @param mapping_df Optional data frame containing the mapping of metabolite IDs 
#'        (KEGG_ID) to their corresponding STITCH IDs and PubChem CIDs.
#' @param stitch_df Optional data frame containing STITCH interaction data, including 
#'        STITCH IDs and chemical-chemical interaction information downloaded from 
#'        the STITCH database.
#' @param output_dir Optional directory path for saving output files. If NULL (default),
#'        no files are written.
#' @param save_excel Logical indicating whether to save results as Excel files 
#'        (default = FALSE).
#' @param analysis_type Character vector specifying which analyses to run. Options include:
#'        "enrichment" (pathway enrichment), "gsea" (metabolite set enrichment), 
#'        "centrality" (network centrality), "network" (pathway-metabolite network),
#'        "interaction" (chemical interaction network), "heatmap" (enrichment heatmap),
#'        "membership" (pathway membership matrix). Default: all analyses.
#' @param run_plots Logical indicating whether to generate visualization plots 
#'        (default = TRUE).
#'
#' @return A list containing results from the specified analyses. Possible components:
#' \itemize{
#'   \item \code{pathway_enrichment_results} - Data frame of pathway enrichment results
#'   \item \code{gsea_results} - Data frame of GSEA analysis results
#'   \item \code{metabolite_centrality} - Data frame of centrality analysis results
#'   \item \code{pathway_plot} - ggplot object for pathway enrichment visualization
#'   \item \code{impact_plot} - ggplot object for impact vs significance
#'   \item \code{gsea_plot} - ggplot object for GSEA results
#'   \item \code{rbc_plot} - ggplot object for relative betweenness centrality
#'   \item \code{network_plot} - ggraph object for metabolite-pathway network
#'   \item \code{heatmap_plot} - ComplexHeatmap object for enrichment significance
#'   \item \code{membership_plot} - ComplexHeatmap object for pathway membership
#'   \item \code{interaction_plot} - ggraph object for chemical interaction network
#' }
#'
#' @details
#' This function performs comprehensive pathway analysis for metabolomics data.
#' Users can run the complete workflow or call individual modular functions for
#' specific analyses. See the individual function documentation for more details.
#'
#' @seealso
#' Individual analysis functions:
#' \code{\link{prepare_gmt_data}} for GMT format conversion,
#' \code{\link{calculate_metabolite_centrality}} for network centrality,
#' \code{\link{perform_enrichment_analysis}} for pathway enrichment,
#' \code{\link{perform_gsea_analysis}} for metabolite set enrichment,
#' \code{\link{create_enrichment_plot}} for pathway visualization,
#' \code{\link{create_impact_plot}} for impact vs significance,
#' \code{\link{create_gsea_plot}} for GSEA results visualization,
#' \code{\link{create_centrality_plot}} for centrality visualization
#'
#' @examples
#' #' ## ** Complete workflow example
#'
#' # Generate example data
#' set.seed(1234)
#' inputMetabolites <- paste0("M", 1:20)
#'
#' pathway_names <- paste0("Pathway", 1:50)
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = rep(pathway_names, each = 1),
#'   Metabolites = sapply(1:50, function(x)
#'     paste(sample(inputMetabolites, sample(5:15, 1)), collapse = ","))
#' )
#'
#' # Add specific pathway examples
#' new_rows <- data.frame(
#'   Pathway = c("Pathway101", "Pathway102", "Pathway103", "Pathway104", "Pathway105"),
#'   Metabolites = c(
#'     "M12,M13,M14,M15,M16,M1,M18,M3,M29,M6,M16,M4",
#'     "M6,M7,M8,M9,M10,M11,M9,M29,M6,M6,M16,M4",
#'     "M24,M25,M26,M27,M28,M29,M30,M29,M26,M5",
#'     "M13,M14,M15,M16,M17,M24,M27,M14",
#'     "M15,M16,M17,M18,M19,M20,M21,M4,M8,M10"
#'   )
#' )
#' PathwayVsMetabolites <- rbind(PathwayVsMetabolites, new_rows)
#'
#' # Generate example metabolite statistics for GSEA
#' example_data <- data.frame(
#'   met_id = inputMetabolites,
#'   pval = runif(20, 0.001, 0.05),
#'   log2fc = rnorm(20, mean = 0, sd = 1)
#' )
#'
#' # Create mapping data for STITCH network
#' set.seed(42)
#' mapping_df <- data.frame(
#'   KEGG_ID = inputMetabolites,
#'   PubChem_CID = as.character(sample(10000:99999, length(inputMetabolites))),
#'   STITCH_ID = paste0("CIDs", stringr::str_pad(sample(1000:9999, length(inputMetabolites)), 8, pad = "0"))
#' )
#'
#' # Create synthetic STITCH interaction data
#' stitch_ids <- mapping_df$STITCH_ID
#' stitch_pairs <- expand.grid(chemical1 = stitch_ids, chemical2 = stitch_ids) %>%
#'   dplyr::filter(chemical1 != chemical2)
#'
#' set.seed(123)
#' stitch_df <- stitch_pairs %>%
#'   dplyr::slice_sample(n = 200) %>%
#'   dplyr::mutate(
#'     similarity = runif(dplyr::n(), 0, 1),
#'     experimental = sample(0:500, dplyr::n(), replace = TRUE),
#'     database = sample(c(0, 300, 600, 900), dplyr::n(), replace = TRUE),
#'     textmining = sample(0:1000, dplyr::n(), replace = TRUE),
#'     combined_score = similarity * 200 + experimental + database + textmining
#'   ) %>%
#'   tibble::as_tibble()
#'
#' # Run complete analysis with all visualizations
#' results <- enrichmet(
#'   inputMetabolites = inputMetabolites,
#'   PathwayVsMetabolites = PathwayVsMetabolites,
#'   example_data = example_data,
#'   top_n = 10,
#'   mapping_df = mapping_df,
#'   stitch_df = stitch_df,
#'   analysis_type = c("enrichment", "gsea", "centrality", "network", "heatmap", "membership", "interaction")
#' )
#'
#' # Access individual results
#' enrichment_results <- results$pathway_enrichment_results
#' head(enrichment_results)
#'
#' # View all plots (these will display in the help)
#' results$pathway_plot
#' results$impact_plot
#' results$gsea_plot
#' results$rbc_plot
#' results$network_plot
#' results$heatmap_plot
#' results$membership_plot
#' results$interaction_plot
#' @export
enrichmet <- function(inputMetabolites,
                      PathwayVsMetabolites,
                      example_data = NULL,
                      top_n = 100,
                      p_value_cutoff = 1,
                      kegg_lookup = NULL,
                      mapping_df = NULL,
                      stitch_df = NULL,
                      output_dir = NULL,
                      save_excel = FALSE,
                      analysis_type = c("enrichment", "gsea", "centrality", "network", "heatmap", "membership", "interaction"),
                      run_plots = TRUE) {
    
    # Input validation
    if (!is.character(inputMetabolites) || length(inputMetabolites) == 0) {
        stop("inputMetabolites must be a non-empty character vector")
    }
    
    if (!is.data.frame(PathwayVsMetabolites) || 
        !all(c("Pathway", "Metabolites") %in% colnames(PathwayVsMetabolites))) {
        stop("PathwayVsMetabolites must be a data frame with 'Pathway' and 'Metabolites' columns")
    }
    
    if (!is.null(output_dir) && !dir.exists(output_dir)) {
        stop("output_dir must be NULL or an existing directory")
    }
    
    valid_analyses <- c("enrichment", "gsea", "centrality", "network", "heatmap", "membership", "interaction")
    if (!all(analysis_type %in% valid_analyses)) {
        stop("analysis_type must contain only: ", paste(valid_analyses, collapse = ", "))
    }
    
    # Initialize results list
    results <- list()
    
    # In the enrichment analysis section:
    if ("enrichment" %in% analysis_type) {
        message("Running pathway enrichment analysis...")
        
        significant_results_df <- perform_enrichment_analysis(
            inputMetabolites = inputMetabolites,
            PathwayVsMetabolites = PathwayVsMetabolites,
            top_n = top_n,
            p_value_cutoff = p_value_cutoff
        )
        
        results$pathway_enrichment_results <- significant_results_df
        
        if (run_plots && nrow(significant_results_df) > 0) {
            results$pathway_plot <- create_enrichment_plot(significant_results_df)  # Consistent naming
            results$impact_plot <- create_impact_plot(significant_results_df)       # Consistent naming
        }
    }
    
    # In the GSEA section:
    if ("gsea" %in% analysis_type) {
        if (is.null(example_data)) {
            warning("GSEA analysis requested but example_data not provided")
        } else {
            message("Running GSEA analysis...")
            
            gsea_results <- perform_gsea_analysis(example_data, PathwayVsMetabolites)
            results$gsea_results <- gsea_results
            
            if (run_plots && nrow(gsea_results) > 0) {
                results$gsea_plot <- create_gsea_plot(gsea_results)  # Consistent naming
            }
        }
    }
    # Centrality Analysis
    if ("centrality" %in% analysis_type) {
        message("Running centrality analysis...")
        
        # Calculate centrality for ALL metabolites first
        all_centrality <- calculate_metabolite_centrality(PathwayVsMetabolites)
        
        # Then filter for input metabolites only
        centrality_results <- all_centrality %>%
            dplyr::filter(Metabolite %in% inputMetabolites) %>%
            dplyr::arrange(desc(RBC_Metabolite))
        
        # Add KEGG names if available
        if (!is.null(kegg_lookup)) {
            centrality_results <- centrality_results %>%
                dplyr::left_join(kegg_lookup, by = c("Metabolite" = "kegg_id")) %>%
                dplyr::mutate(Display_Name = ifelse(!is.na(name), Metabolite)) %>%
                dplyr::select(-name)
        } else {
            centrality_results$Display_Name <- centrality_results$Metabolite
        }
        
        results$metabolite_centrality <- centrality_results
        
        if (run_plots && nrow(centrality_results) > 0) {
            results$rbc_plot <- create_centrality_plot(centrality_results, kegg_lookup = kegg_lookup)
        }
    }
    
    # Network Analysis
    if ("network" %in% analysis_type && run_plots) {
        message("Generating metabolite-pathway network visualization...")
        results$network_plot <- create_network_plot(inputMetabolites, PathwayVsMetabolites, kegg_lookup)
    }
    
    # Enrichment Heatmap
    if ("heatmap" %in% analysis_type && run_plots) {
        message("Generating enrichment heatmap...")
        if ("enrichment" %in% analysis_type && !is.null(results$pathway_enrichment_results)) {
            results$heatmap_plot <- create_heatmap_plot(
                results$pathway_enrichment_results, 
                PathwayVsMetabolites, 
                inputMetabolites, 
                kegg_lookup
            )
        }
    }
    
    # Pathway Membership Plot
    if ("membership" %in% analysis_type && run_plots) {
        message("Generating pathway membership plot...")
        results$membership_plot <- create_membership_plot(PathwayVsMetabolites, inputMetabolites, kegg_lookup)
    }
    
    # STITCH Interaction Network Analysis
    if ("interaction" %in% analysis_type && run_plots) {
        message("Generating STITCH interaction network...")
        results$interaction_plot <- create_interaction_plot(inputMetabolites, mapping_df, stitch_df, kegg_lookup)
    }
    
    # Save results to files if requested
    if (save_excel && !is.null(output_dir)) {
        if ("enrichment" %in% analysis_type && !is.null(results$pathway_enrichment_results)) {
            output_file1 <- file.path(output_dir, "pathway_enrichment_results.xlsx")
            openxlsx::write.xlsx(results$pathway_enrichment_results, output_file1)
            message("Pathway enrichment results saved to: ", output_file1)
        }
        
        if ("gsea" %in% analysis_type && !is.null(results$gsea_results)) {
            output_file2 <- file.path(output_dir, "gsea_results.xlsx")
            openxlsx::write.xlsx(results$gsea_results, output_file2)
            message("GSEA results saved to: ", output_file2)
        }
        
        if ("centrality" %in% analysis_type && !is.null(results$metabolite_centrality)) {
            output_file3 <- file.path(output_dir, "metabolite_centrality_results.xlsx")
            openxlsx::write.xlsx(results$metabolite_centrality, output_file3)
            message("Centrality results saved to: ", output_file3)
        }
    }
    
    return(results)
}
