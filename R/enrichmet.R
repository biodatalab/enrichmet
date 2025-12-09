# utils.R or anywhere near the top of your R script
utils::globalVariables(
    c("G1",
      "G2",
      "P.Value",
      "adj.P.Val",
      "logFC",
      "AveExpr",
      "ave_expr",
      "Count",
      "Min_Adjusted_P",
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
        "logp_vec",
        "log2fc",
        "Significant",
        "kegg_id",
        "padj"
    )
)
# Add these import statements
#' @import dplyr
#' @import tidyr
#' @importFrom tibble as_tibble
#' @importFrom stringr str_pad
NULL
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
#'        Can handle complex formats like "C00042|C02170" or "C08356|C00137|C00936".
#'        Alternatively, can be a data frame from run_de()$kegg_ready containing
#'        metabolite statistics. Required if da_results is not provided.
#' @param PathwayVsMetabolites A data frame containing pathways and their associated 
#'        metabolites. Must include columns 'Pathway' and 'Metabolites'.
#' @param example_data A data frame containing metabolite-level data for GSEA analysis. 
#'        Should include columns "met_id", "pval", and "log2fc". Required for GSEA analysis.
#' @param da_results Optional output from run_de() function. If provided, 
#'        inputMetabolites will be extracted from da_results$kegg_ready and used for
#'        enrichment analysis of significant metabolites. Mutually exclusive with inputMetabolites.
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
#' @param network_top_n Number of top pathways to include in network visualization (default = 20)
#' @param heatmap_top_n Number of top pathways to include in heatmap visualization (default = 20)
#' @param membership_top_n Number of top pathways to include in membership matrix (default = 20)
#' @param min_pathway_occurrence Minimum occurrence threshold for pathways in heatmap/membership (default = 1)
#' @param min_metabolite_occurrence Minimum occurrence threshold for metabolites in heatmap/membership (default = 1)
#' @param use_significant_only Logical indicating whether to use only significant metabolites
#'        from da_results for enrichment analysis (default = TRUE). Only applies when 
#'        da_results is provided.
#' @param split_complex_ids Logical indicating whether to split complex KEGG IDs like 
#'        "C00042|C02170" into individual metabolites (default = TRUE).
#' @param significance_threshold Significance threshold for filtering metabolites when
#'        using da_results. Options: "up", "down", "both" (default: "both").
#' @param include_volcano Logical indicating whether to include volcano plot when 
#'        da_results is provided (default = TRUE).
#' @param fc_cutoff_up Fold change cutoff for upregulated metabolites (default = 1, i.e., 2-fold).
#' @param fc_cutoff_down Fold change cutoff for downregulated metabolites (default = -1, i.e., 0.5-fold).
#' @param fdr_cutoff_da FDR cutoff for differential analysis significance (default = 0.05).
#' @param force_custom_filters Logical indicating whether to force using custom FC and FDR 
#'        cutoffs even when "Significant" column exists (default = FALSE).
#'
#' @return A list containing results from the specified analyses. Possible components:
#' \itemize{
#'   \item \code{input_metabolites_used} - Character vector of metabolites used for analysis
#'   \item \code{pathway_enrichment_all} - Data frame of all pathway enrichment results
#'   \item \code{pathway_enrichment_results} - Data frame of filtered pathway enrichment results
#'   \item \code{gsea_results} - Data frame of GSEA analysis results
#'   \item \code{metabolite_centrality} - Data frame of centrality analysis results
#'   \item \code{volcano_plot} - ggplot object for volcano plot (only when da_results provided)
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
#' # Example: Comprehensive enrichment analysis showing all plots
#' 
#' # 1. Create realistic pathway data with real KEGG pathways
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = c("Glycolysis / Gluconeogenesis",
#'               "Citrate cycle (TCA cycle)",
#'               "Pentose phosphate pathway",
#'               "Pyruvate metabolism",
#'               "Amino sugar and nucleotide sugar metabolism",
#'               "Butanoate metabolism",
#'               "Propanoate metabolism",
#'               "Valine, leucine and isoleucine degradation"),
#'   Metabolites = c(
#'     "C00031,C00022,C00197,C00221,C00631,C01172,C00074,C00186",
#'     "C00022,C00036,C00024,C00158,C00149,C00311,C00417,C05125",
#'     "C00117,C00257,C00121,C00085,C00118,C00198,C00231,C00279",
#'     "C00022,C00024,C00036,C00074,C00122,C00186,C00236,C00267",
#'     "C00031,C00095,C00103,C00140,C00159,C00185,C00267,C00645",
#'     "C00122,C00149,C00164,C00232,C00356,C00441,C00533,C01074",
#'     "C00163,C00179,C00233,C00356,C00408,C00533,C00874,C00979",
#'     "C00141,C00183,C00188,C00233,C00322,C00407,C00507,C01048"
#'   )
#' )
#' 
#' # 2. Create input metabolites that will show enrichment
#' inputMetabolites <- c(
#'   "C00031", "C00022", "C00197", "C00221", "C00631", 
#'   "C01172", "C00074", "C00186", "C00036", "C00158"
#' )
#' 
#' # 3. Create KEGG lookup table for better metabolite names
#' kegg_lookup <- data.frame(
#'   kegg_id = c(
#'     "C00031", "C00022", "C00197", "C00221", "C00631", "C01172",
#'     "C00074", "C00186", "C00036", "C00158", "C00024", "C00149",
#'     "C00311", "C00117", "C00257", "C00121", "C00085", "C00118"
#'   ),
#'   name = c(
#'     "D-Glucose 6-phosphate", "Oxaloacetate", 
#'     "D-Glyceraldehyde 3-phosphate", "D-Glucose", 
#'     "Glycerone phosphate", "sn-Glycerol 3-phosphate",
#'     "Phosphoenolpyruvate", "3-Phospho-D-glycerate", 
#'     "2-Oxoglutarate", "Citrate", "Acetyl-CoA", "Acetoacetyl-CoA",
#'     "Isocitrate", "D-Ribulose 5-phosphate", 
#'     "Sedoheptulose 7-phosphate", "D-Xylulose 5-phosphate",
#'     "D-Glucose 6-phosphate", "D-Fructose 6-phosphate"
#'   )
#' )
#' 
#' # 4. Create example metabolite statistics for GSEA
#' # This data simulates differential expression results
#' set.seed(123)
#' n_metab <- 50
#' example_data <- data.frame(
#'   met_id = c(
#'     paste0("C", sprintf("%05d", 1:30)),
#'     "C00031", "C00022", "C00197", "C00221", "C00631", 
#'     "C01172", "C00074", "C00186", "C00036", "C00158",
#'     "C00024", "C00149", "C00311", "C00117", "C00257",
#'     "C00121", "C00085", "C00118"
#'   ),
#'   pval = c(
#'     runif(30, 0.001, 0.1),
#'     c(0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 
#'       0.005, 0.006, 0.007, 0.008, 0.01, 0.02, 0.03, 
#'       0.04, 0.05, 0.06, 0.07, 0.08)
#'   ),
#'   log2fc = c(
#'     rnorm(30, mean = 0, sd = 0.5),
#'     c(2.5, 2.3, 2.1, 1.9, 1.7, 1.5, 1.3, 1.1, 0.9, 0.7,
#'       -0.5, -0.7, -0.9, -1.1, -1.3, -1.5, -1.7, -1.9)
#'   )
#' )
#' 
#' # 5. Create minimal mapping data for STITCH network
#' mapping_df <- data.frame(
#'   KEGG_ID = c(
#'     "C00031", "C00022", "C00197", "C00221", "C00631", 
#'     "C01172", "C00074", "C00186", "C00036", "C00158"
#'   ),
#'   PubChem_CID = c(
#'     "439284", "164619", "729", "79025", "751",
#'     "753", "1005", "724", "51", "311"
#'   ),
#'   STITCH_ID = c(
#'     "CIDs000439284", "CIDs000164619", "CIDs000000729", 
#'     "CIDs000079025", "CIDs000000751", "CIDs000000753",
#'     "CIDs000001005", "CIDs000000724", "CIDs000000051", 
#'     "CIDs000000311"
#'   )
#' )
#' 
#' # 6. Create synthetic STITCH interactions
#' set.seed(42)
#' stitch_combinations <- combn(mapping_df$STITCH_ID, 2)
#' stitch_pairs <- data.frame(
#'   chemical1 = stitch_combinations[1, ],
#'   chemical2 = stitch_combinations[2, ]
#' )
#' 
#' # Randomly select 15 interactions
#' selected_pairs <- stitch_pairs[sample(nrow(stitch_pairs), 15), ]
#' 
#' stitch_df <- data.frame(
#'   chemical1 = selected_pairs$chemical1,
#'   chemical2 = selected_pairs$chemical2,
#'   similarity = runif(15, 0.6, 0.95),
#'   experimental = sample(200:800, 15, replace = TRUE),
#'   database = sample(c(0, 300, 600, 900), 15, replace = TRUE),
#'   textmining = sample(0:900, 15, replace = TRUE)
#' )
#' 
#' stitch_df$combined_score <- with(
#'   stitch_df, 
#'   similarity * 200 + experimental + database + textmining
#' )
#' 
#' # 7. Run comprehensive enrichment analysis with ALL plot types
#' 
#' if (requireNamespace("igraph", quietly = TRUE) &&
#'     requireNamespace("ggraph", quietly = TRUE) &&
#'     requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#'     
#'   results <- enrichmet(
#'     inputMetabolites = inputMetabolites,
#'     PathwayVsMetabolites = PathwayVsMetabolites,
#'     example_data = example_data,
#'     kegg_lookup = kegg_lookup,
#'     mapping_df = mapping_df,
#'     stitch_df = stitch_df,
#'     top_n = 10,
#'     p_value_cutoff = 1,
#'     analysis_type = c(
#'       "enrichment", "gsea", "centrality", "network",
#'       "heatmap", "membership", "interaction"
#'     ),
#'     network_top_n = 8,
#'     heatmap_top_n = 8,
#'     membership_top_n = 8,
#'     min_pathway_occurrence = 2,
#'     min_metabolite_occurrence = 1
#'   )
#'   
#'   # Print summary
#'   cat("=== ENRICHMENT ANALYSIS RESULTS ===\n")
#'   cat("Input metabolites used:", 
#'       length(results$input_metabolites_used), "\n")
#'   cat("Pathways tested:", 
#'       nrow(results$pathway_enrichment_all), "\n")
#'   cat("Significant pathways (p < 0.05):", 
#'       nrow(results$pathway_enrichment_results), "\n")
#'   
#'   if (!is.null(results$pathway_enrichment_results)) {
#'     cat("\n=== TOP 5 ENRICHED PATHWAYS ===\n")
#'     print(head(
#'       results$pathway_enrichment_results[
#'         , c("Pathway", "P_value", "Adjusted_P_value", "Enrichment_Ratio")
#'       ], 
#'       5
#'     ))
#'   }
#'   
#'   if (!is.null(results$gsea_results)) {
#'     cat("\n=== TOP 5 GSEA PATHWAYS ===\n")
#'     print(head(
#'       results$gsea_results[
#'         , c("pathway", "pval", "padj", "NES")
#'       ], 
#'       5
#'     ))
#'   }
#'   
#'   if (!is.null(results$metabolite_centrality)) {
#'     cat("\n=== TOP 5 CENTRAL METABOLITES ===\n")
#'     print(head(
#'       results$metabolite_centrality[
#'         , c("Display_Name", "RBC_Metabolite")
#'       ], 
#'       5
#'     ))
#'   }
#'   
#'   # Display all available plots
#'   cat("\n=== AVAILABLE PLOTS ===\n")
#'   available_plots <- names(results)[sapply(results, function(x) {
#'     any(class(x) %in% c("ggplot", "gg", "ggraph", "Heatmap", "HeatmapList"))
#'   })]
#'   cat("Plots generated:", paste(available_plots, collapse = ", "), "\n")
#'   
#'   # Show key plots (if they exist)
#'   if (!is.null(results$pathway_plot)) {
#'     cat("\nDisplaying pathway enrichment plot...\n")
#'     print(results$pathway_plot)
#'   }
#'   
#'   if (!is.null(results$impact_plot)) {
#'     cat("\nDisplaying impact plot...\n")
#'     print(results$impact_plot)
#'   }
#'   
#'   if (!is.null(results$gsea_plot)) {
#'     cat("\nDisplaying GSEA plot...\n")
#'     print(results$gsea_plot)
#'   }
#'   
#'   if (!is.null(results$rbc_plot)) {
#'     cat("\nDisplaying centrality plot...\n")
#'     print(results$rbc_plot)
#'   }
#' 
#'  if (!is.null(results$heatmap_plot)) {
#'     cat("\nDisplaying heatmap plot...\n")
#'     print(results$heatmap_plot)
#'   }
#'   
#'   if (!is.null(results$membership_plot)) {
#'     cat("\nDisplaying membership plot...\n")
#'     print(results$membership_plot)
#'   }
#'   
#'   if (!is.null(results$interaction_plot)) {
#'     cat("\nDisplaying interaction plot...\n")
#'     print(results$interaction_plot)
#'   }
#'   
#'   if (!is.null(results$network_plot)) {
#'     cat("\nDisplaying network plot...\n")
#'     print(results$network_plot)
#'   }
#' }
#' 
#' @export
enrichmet <- function(inputMetabolites = NULL,
                      PathwayVsMetabolites,
                      example_data = NULL,
                      da_results = NULL,
                      top_n = 100,
                      p_value_cutoff = 1,
                      kegg_lookup = NULL,
                      mapping_df = NULL,
                      stitch_df = NULL,
                      output_dir = NULL,
                      save_excel = FALSE,
                      analysis_type = c("enrichment", "gsea", "centrality", "network", "heatmap", "membership", "interaction"),
                      run_plots = TRUE,
                      network_top_n = 20,
                      heatmap_top_n = 20,
                      membership_top_n = 20,
                      min_pathway_occurrence = 1,    
                      min_metabolite_occurrence = 1,
                      use_significant_only = TRUE,
                      split_complex_ids = TRUE,
                      significance_threshold = "both",
                      include_volcano = TRUE,
                      fc_cutoff_up = 1,
                      fc_cutoff_down = -1,
                      fdr_cutoff_da = 0.05,
                      force_custom_filters = FALSE) {
    
    # ---- Helper function to process complex KEGG IDs ----
    process_kegg_ids <- function(ids) {
        if (is.null(ids)) return(NULL)
        
        # Split by | and remove empty strings
        all_ids <- unlist(strsplit(ids, "\\|"))
        # Remove any empty strings that might result from consecutive ||
        all_ids <- all_ids[all_ids != ""]
        # Remove duplicates and return
        unique(all_ids)
    }
    
    # ---- Input validation and metabolite extraction ----
    metabolites_to_use <- NULL
    
    # Case 1: da_results provided
    if (!is.null(da_results)) {
        if (!is.null(inputMetabolites)) {
            warning("Both da_results and inputMetabolites provided. Using da_results and ignoring inputMetabolites.")
        }
        
        message("Using metabolites from da_results")
        
        if (!"kegg_ready" %in% names(da_results)) {
            stop("da_results must contain 'kegg_ready' element")
        }
        
        kegg_data <- da_results$kegg_ready
        
        if (use_significant_only) {
            # Check if we should use existing Significant column or apply custom filters
            use_existing_significant <- "Significant" %in% colnames(kegg_data) && !force_custom_filters
            
            if (use_existing_significant) {
                # Use existing Significant column if available and not forcing custom filters
                message("Using existing 'Significant' column from DA results")
                
                if (significance_threshold == "up") {
                    metabolites_to_use <- kegg_data$kegg_id[kegg_data$Significant == "Up"]
                    message("Using ", length(metabolites_to_use), " upregulated metabolites from existing Significant column")
                } else if (significance_threshold == "down") {
                    metabolites_to_use <- kegg_data$kegg_id[kegg_data$Significant == "Down"]
                    message("Using ", length(metabolites_to_use), " downregulated metabolites from existing Significant column")
                } else {
                    metabolites_to_use <- kegg_data$kegg_id[kegg_data$Significant != "Not significant"]
                    message("Using ", length(metabolites_to_use), " significant metabolites (both up and down) from existing Significant column")
                }
                
                # Warn if custom filters were provided but not used
                if (fc_cutoff_up != 1 || fc_cutoff_down != -1 || fdr_cutoff_da != 0.05) {
                    message("Note: Custom FC/FDR filters provided but using existing Significant column. ",
                            "Set force_custom_filters = TRUE to apply custom filters.")
                }
                
            } else {
                # Apply custom fold change and FDR cutoffs
                message("Applying custom significance filters: FC up > ", fc_cutoff_up, 
                        ", FC down < ", fc_cutoff_down, ", FDR < ", fdr_cutoff_da)
                
                if (significance_threshold == "up") {
                    metabolites_to_use <- kegg_data$kegg_id[
                        kegg_data$log2fc > fc_cutoff_up & kegg_data$padj < fdr_cutoff_da
                    ]
                    message("Using ", length(metabolites_to_use), " upregulated metabolites (FDR < ", fdr_cutoff_da, ", FC > ", fc_cutoff_up, ")")
                } else if (significance_threshold == "down") {
                    metabolites_to_use <- kegg_data$kegg_id[
                        kegg_data$log2fc < fc_cutoff_down & kegg_data$padj < fdr_cutoff_da
                    ]
                    message("Using ", length(metabolites_to_use), " downregulated metabolites (FDR < ", fdr_cutoff_da, ", FC < ", fc_cutoff_down, ")")
                } else {
                    metabolites_to_use <- kegg_data$kegg_id[
                        (kegg_data$log2fc > fc_cutoff_up | kegg_data$log2fc < fc_cutoff_down) & 
                            kegg_data$padj < fdr_cutoff_da
                    ]
                    message("Using ", length(metabolites_to_use), " significant metabolites (FDR < ", fdr_cutoff_da, 
                            ", |FC| > ", abs(fc_cutoff_up), ")")
                }
            }
        } else {
            # Use all metabolites from kegg_ready
            metabolites_to_use <- kegg_data$kegg_id
            message("Using all ", length(metabolites_to_use), " metabolites from DA results")
        }
        
        # Remove NA values from kegg_id
        metabolites_to_use <- metabolites_to_use[!is.na(metabolites_to_use)]
        metabolites_to_use <- metabolites_to_use[metabolites_to_use != ""]
        
        # Store the DA results for potential use in GSEA
        if (is.null(example_data) && "full_results" %in% names(da_results)) {
            example_data <- da_results$full_results
            message("Using da_results$full_results for GSEA analysis")
        }
        
    } 
    # Case 2: inputMetabolites provided directly (no da_results)
    else if (!is.null(inputMetabolites)) {
        if (is.data.frame(inputMetabolites)) {
            # If inputMetabolites is a data frame (like kegg_ready)
            if ("kegg_id" %in% colnames(inputMetabolites)) {
                metabolites_to_use <- inputMetabolites$kegg_id
                message("Using ", length(metabolites_to_use), " metabolites from inputMetabolites data frame")
                
                # If example_data not provided, use this for GSEA if it has required columns
                if (is.null(example_data) && all(c("pval", "log2fc") %in% colnames(inputMetabolites))) {
                    example_data <- inputMetabolites
                    message("Using inputMetabolites data frame for GSEA analysis")
                }
            } else if ("met_id" %in% colnames(inputMetabolites)) {
                # Fallback to met_id if kegg_id not available
                metabolites_to_use <- inputMetabolites$met_id
                message("Using ", length(metabolites_to_use), " metabolites from inputMetabolites data frame (met_id column)")
            } else {
                stop("If inputMetabolites is a data frame, it must contain 'kegg_id' or 'met_id' column")
            }
        } else if (is.character(inputMetabolites)) {
            # If inputMetabolites is a character vector
            metabolites_to_use <- inputMetabolites
            message("Using ", length(metabolites_to_use), " metabolites from inputMetabolites character vector")
        } else {
            stop("inputMetabolites must be a character vector, data frame with 'kegg_id' or 'met_id' column, or NULL when da_results is provided")
        }
    } 
    # Case 3: Neither provided
    else {
        stop("Either inputMetabolites or da_results must be provided")
    }
    
    # ---- Process complex KEGG IDs ----
    if (split_complex_ids && !is.null(metabolites_to_use)) {
        message("Processing complex KEGG IDs (splitting by |)...")
        original_count <- length(metabolites_to_use)
        
        # Split all complex IDs and flatten
        all_split_ids <- unlist(lapply(metabolites_to_use, process_kegg_ids))
        
        # Remove duplicates
        metabolites_to_use <- unique(all_split_ids)
        
        message("Split ", original_count, " complex IDs into ", length(metabolites_to_use), " unique KEGG IDs")
    }
    
    # Final validation
    if (length(metabolites_to_use) == 0) {
        stop("No metabolites found for analysis. Check your input data.")
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
    
    # ---- Initialize results ----
    results <- list()
    results$input_metabolites_used <- metabolites_to_use
    
    # ---- Include volcano plot from DA results if requested ----
    if (!is.null(da_results) && include_volcano && "volcano_plot" %in% names(da_results)) {
        message("Including volcano plot from differential analysis results")
        results$volcano_plot <- da_results$volcano_plot
    }
    
    # ---- Enrichment Analysis ----
    if ("enrichment" %in% analysis_type) {
        message("Running pathway enrichment analysis...")
        
        # Compute all pathways first
        all_enrichment <- perform_enrichment_analysis(
            inputMetabolites = metabolites_to_use,
            PathwayVsMetabolites = PathwayVsMetabolites,
            top_n = NULL,       # keep all
            p_value_cutoff = 1  # keep all
        )
        
        results$pathway_enrichment_all <- all_enrichment
        
        # Filter for significant pathways
        significant_results <- all_enrichment %>%
            dplyr::filter(Adjusted_P_value < p_value_cutoff) %>%
            dplyr::arrange(desc(Log_P_value))
        
        if (!is.null(top_n) && nrow(significant_results) > top_n) {
            significant_results <- head(significant_results, top_n)
        }
        
        results$pathway_enrichment_results <- significant_results
        
        # Optional plots - ALWAYS CREATE BOTH PATHWAY AND IMPACT PLOTS
        if (run_plots && nrow(significant_results) > 0) {
            results$pathway_plot <- create_enrichment_plot(significant_results)
            results$impact_plot <- create_impact_plot(significant_results)
        }
    }
    
    # ---- GSEA Analysis ----
    if ("gsea" %in% analysis_type) {
        if (is.null(example_data)) {
            warning("GSEA analysis requested but example_data not provided")
        } else {
            message("Running GSEA analysis...")
            gsea_results <- perform_gsea_analysis(example_data, PathwayVsMetabolites)
            results$gsea_results <- gsea_results
            if (run_plots && nrow(gsea_results) > 0) {
                results$gsea_plot <- create_gsea_plot(gsea_results, 
                                                      top_n = 20, 
                                                      kegg_lookup = kegg_lookup)
            }
        }
    }
    
    # ---- Centrality Analysis ----
    if ("centrality" %in% analysis_type) {
        message("Running centrality analysis...")
        all_centrality <- calculate_metabolite_centrality(PathwayVsMetabolites)
        centrality_results <- all_centrality %>%
            dplyr::filter(Metabolite %in% metabolites_to_use) %>%
            dplyr::arrange(desc(RBC_Metabolite))
        
        if (!is.null(kegg_lookup)) {
            centrality_results <- centrality_results %>%
                dplyr::left_join(kegg_lookup, by = c("Metabolite" = "kegg_id")) %>%
                dplyr::mutate(Display_Name = ifelse(!is.na(name), name, Metabolite)) %>%
                dplyr::select(-name)
        } else {
            centrality_results$Display_Name <- centrality_results$Metabolite
        }
        
        results$metabolite_centrality <- centrality_results
        if (run_plots && nrow(centrality_results) > 0) {
            results$rbc_plot <- create_centrality_plot(centrality_results, kegg_lookup = kegg_lookup)
        }
    }
    
    # ---- Network ----
    if ("network" %in% analysis_type && run_plots) {
        message("Generating metabolite-pathway network visualization...")
        results$network_plot <- create_network_plot(
            metabolites_to_use, 
            PathwayVsMetabolites, 
            kegg_lookup,
            top_n = network_top_n
        )
    }
    
    # ---- Heatmap ----
    if ("heatmap" %in% analysis_type && run_plots) {
        message("Generating enrichment heatmap...")
        if ("enrichment" %in% analysis_type && !is.null(results$pathway_enrichment_results)) {
            results$heatmap_plot <- create_heatmap_plot(
                results$pathway_enrichment_results, 
                PathwayVsMetabolites, 
                metabolites_to_use, 
                kegg_lookup,
                top_n = heatmap_top_n,
                min_pathways = min_pathway_occurrence,
                min_metabolites = min_metabolite_occurrence
            )
        }
    }
    
    # ---- Membership ----
    if ("membership" %in% analysis_type && run_plots) {
        message("Generating pathway membership plot...")
        results$membership_plot <- create_membership_plot(
            PathwayVsMetabolites, 
            metabolites_to_use, 
            kegg_lookup,
            top_n = membership_top_n,
            min_pathway_occurrence = min_pathway_occurrence,
            min_metabolite_occurrence = min_metabolite_occurrence
        )
    }
    
    # ---- STITCH interaction ----
    if ("interaction" %in% analysis_type && run_plots) {
        message("Generating STITCH interaction network...")
        results$interaction_plot <- create_interaction_plot(metabolites_to_use, mapping_df, stitch_df, kegg_lookup)
    }
    
    # ---- Save Excel ----
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