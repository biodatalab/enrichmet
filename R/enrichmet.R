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
      "padj",
      "KEGG",
      "Reaction",
      "shared_reactions",
      # used by name / mixed mapping
      "name_lc"
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
#' Input can be KEGG compound IDs, metabolite names, or a mix of both (controlled by
#' \code{input_type}). Name mapping uses \code{kegg_lookup} (e.g. from
#' \code{\link{fetch_kegg_compound_lookup}}) and matches the behaviour of the Shiny app.
#'
#' @param inputMetabolites A character vector of metabolite identifiers and/or names,
#'        or a data frame from \code{run_de()$kegg_ready}. For \code{input_type = "kegg"}
#'        (default), values should be KEGG IDs (e.g. \code{"C00031"}); complex forms like
#'        \code{"C00042|C02170"} are supported when \code{split_complex_ids = TRUE}.
#'        For \code{input_type = "name"} or \code{"mixed"}, common names (e.g. \code{"glucose"})
#'        are mapped to KEGG IDs via \code{kegg_lookup}. Required if \code{da_results} is not provided.
#' @param PathwayVsMetabolites A data frame containing pathways and their associated
#'        metabolites. Must include columns \code{Pathway} and \code{Metabolites}.
#' @param example_data A data frame containing metabolite-level data for GSEA analysis.
#'        Should include columns \code{met_id}, \code{pval}, and \code{log2fc}. Required for GSEA.
#' @param da_results Optional output from \code{run_de()}. If provided, metabolites are
#'        taken from \code{da_results$kegg_ready}. Mutually exclusive with prioritising
#'        \code{inputMetabolites} when both are supplied (da_results wins).
#' @param top_n An integer specifying the number of top pathways to include in the
#'        pathway enrichment results (default is 100). Use \code{NULL} to return all pathways.
#' @param p_value_cutoff A numeric value for the adjusted p-value threshold for
#'        filtering significant pathways (default is 1, no filtering).
#' @param kegg_lookup Optional data frame for KEGG ID ↔ name mapping. Must contain
#'        columns \code{kegg_id} and \code{name}. **Required** when
#'        \code{input_type} is \code{"name"} or \code{"mixed"} (e.g.
#'        \code{fetch_kegg_compound_lookup()}).
#' @param reactome_df Optional data frame containing the KEGG-to-Reactome
#'        reaction mapping. Must contain columns \code{KEGG} and \code{Reaction}.
#' @param output_dir Optional directory path for saving output files. If \code{NULL}
#'        (default), no files are written.
#' @param save_excel Logical; whether to save results as Excel files (default = FALSE).
#' @param analysis_type Character vector specifying which analyses to run. Options:
#'        \code{"enrichment"}, \code{"gsea"}, \code{"centrality"}, \code{"network"},
#'        \code{"interaction"}, \code{"heatmap"}, \code{"membership"}. Default: all.
#' @param run_plots Logical; whether to generate visualization plots (default = TRUE).
#' @param network_top_n Number of top pathways in the network plot (default = 20).
#' @param heatmap_top_n Number of top pathways in the heatmap (default = 20).
#' @param membership_top_n Number of top pathways in the membership matrix (default = 20).
#' @param min_pathway_occurrence Minimum pathway occurrence for heatmap/membership (default = 1).
#' @param min_metabolite_occurrence Minimum metabolite occurrence for heatmap/membership (default = 1).
#' @param use_significant_only Logical; use only significant metabolites from
#'        \code{da_results} (default = TRUE). Only applies when \code{da_results} is provided.
#' @param split_complex_ids Logical; split complex KEGG IDs like \code{"C00042|C02170"}
#'        into individual IDs (default = TRUE).
#' @param significance_threshold When using \code{da_results}: \code{"up"}, \code{"down"},
#'        or \code{"both"} (default: \code{"both"}).
#' @param include_volcano Logical; include volcano plot from \code{da_results} (default = TRUE).
#' @param fc_cutoff_up Fold-change cutoff for upregulated metabolites (default = 1).
#' @param fc_cutoff_down Fold-change cutoff for downregulated metabolites (default = -1).
#' @param fdr_cutoff_da FDR cutoff for DA significance (default = 0.05).
#' @param force_custom_filters Logical; force custom FC/FDR even if a \code{Significant}
#'        column exists (default = FALSE).
#' @param backgroundMetabolites Optional character vector of KEGG IDs defining the
#'        statistical background (universe) for background-corrected enrichment.
#'        If \code{NULL}, the full pathway database is used.
#' @param input_type Character; how to interpret \code{inputMetabolites}:
#'        \itemize{
#'          \item \code{"kegg"} (default) – treat entries as KEGG IDs (\code{C#####}).
#'          \item \code{"name"} – treat entries as metabolite names and map them to
#'                KEGG IDs via \code{kegg_lookup}.
#'          \item \code{"mixed"} – keep valid KEGG IDs and map the rest as names.
#'        }
#'        Same behaviour as the Shiny app input-type control.
#'
#' @return A list containing results from the specified analyses. Possible components:
#' \itemize{
#'   \item \code{input_metabolites_used} – Character vector of KEGG IDs used for analysis
#'   \item \code{pathway_enrichment_all} – Data frame of all pathway enrichment results
#'   \item \code{pathway_enrichment_results} – Data frame of filtered pathway enrichment results
#'   \item \code{gsea_results} – Data frame of GSEA analysis results
#'   \item \code{metabolite_centrality} – Data frame of centrality analysis results
#'   \item \code{volcano_plot} – ggplot object (only when \code{da_results} provided)
#'   \item \code{pathway_plot} – ggplot object for pathway enrichment
#'   \item \code{impact_plot} – ggplot object for impact vs significance
#'   \item \code{gsea_plot} – ggplot object for GSEA results
#'   \item \code{rbc_plot} – ggplot object for relative betweenness centrality
#'   \item \code{network_plot} – ggraph object for metabolite–pathway network
#'   \item \code{heatmap_plot} – ComplexHeatmap object for enrichment significance
#'   \item \code{membership_plot} – ComplexHeatmap object for pathway membership
#'   \item \code{interaction_plot} – ggraph object for Reactome interaction network
#' }
#'
#' @details
#' This function performs comprehensive pathway analysis for metabolomics data.
#' Users can run the complete workflow or call individual modular functions for
#' specific analyses. See the individual function documentation for more details.
#'
#' When \code{input_type} is \code{"name"} or \code{"mixed"}, names are matched
#' case-insensitively against \code{kegg_lookup$name} (exact, then starts-with,
#' then contains). Unmapped entries are dropped with a warning. All downstream
#' analyses (including the Reactome interaction network) receive the resulting
#' clean KEGG ID vector.
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
#' \code{\link{create_centrality_plot}} for centrality visualization,
#' \code{\link{fetch_kegg_compound_lookup}} for building \code{kegg_lookup}
#'
#' @examples
#' 
#' 
#' # Pathway map and compound names (KEGG REST API)
#' PathwayVsMetabolites <- fetch_kegg_pathway_metabolites(organism = "hsa")
#' kegg_lookup <- fetch_kegg_compound_lookup()
#'
#' # Metabolomics matrix (bundled example data)
#' metabolomics_path <- system.file(
#'     "extdata", "example_data.csv",
#'     package = "enrichmet"
#' )
#'
#' if (metabolomics_path == "") {
#'     stop("Example file 'example_data.csv' not found in inst/extdata/")
#' }
#'
#' metabolomics_mat <- read.csv(
#'     metabolomics_path,
#'     row.names = 1,
#'     check.names = FALSE
#' )
#' 
#' da_out <- run_de(
#'     metabolomics_mat, "TK-CMV", "K-CMV",
#'     fc_threshold = 1, pval_threshold = 0.05
#' )
#'
#' # Enrichment from DA results (query set taken from da_out)
#' results <- enrichmet(
#'     inputMetabolites = NULL,
#'     PathwayVsMetabolites = PathwayVsMetabolites,
#'     da_results = da_out,
#'     example_data = example_data,
#'     kegg_lookup = kegg_lookup,
#'     analysis_type = c("enrichment", "gsea", "centrality"),
#'     top_n = 10,
#'     p_value_cutoff = 1
#' )
#'
#' head(results$pathway_enrichment_results, 3)
#' 
#' @export
enrichmet <- function(inputMetabolites = NULL,
                      PathwayVsMetabolites,
                      example_data = NULL,
                      da_results = NULL,
                      top_n = 100,
                      p_value_cutoff = 1,
                      kegg_lookup = NULL,
                      reactome_df = NULL,
                      output_dir = NULL,
                      save_excel = FALSE,
                      analysis_type = c("enrichment", "gsea", "centrality", "network",
                                        "heatmap", "membership", "interaction"),
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
                      force_custom_filters = FALSE,
                      backgroundMetabolites = NULL,
                      input_type = c("kegg", "name", "mixed")) {
    
    input_type <- match.arg(input_type)
    
    # ---- Helper: process complex KEGG IDs ----
    process_kegg_ids <- function(ids) {
        if (is.null(ids)) return(NULL)
        all_ids <- unlist(strsplit(ids, "\\|"))
        all_ids <- all_ids[all_ids != ""]
        unique(all_ids)
    }
    
    # ---- Helper: map names → KEGG (same logic as the Shiny app) ----
    map_names_to_kegg <- function(x, lookup) {
        if (is.null(lookup) || !is.data.frame(lookup) || nrow(lookup) == 0) {
            stop("input_type is '", input_type,
                 "' but kegg_lookup is missing or empty. ",
                 "Pass kegg_lookup = fetch_kegg_compound_lookup().",
                 call. = FALSE)
        }
        if (!all(c("kegg_id", "name") %in% colnames(lookup))) {
            stop("kegg_lookup must contain columns 'kegg_id' and 'name'.", call. = FALSE)
        }
        
        lookup$name_lc <- tolower(trimws(as.character(lookup$name)))
        lookup$kegg_id <- trimws(as.character(lookup$kegg_id))
        
        map_one <- function(val) {
            val <- trimws(as.character(val))
            if (is.na(val) || val == "") return(NA_character_)
            
            # Already a KEGG ID?
            if (grepl("^C[0-9]{5}$", val, ignore.case = TRUE)) {
                return(toupper(val))
            }
            
            val_lc <- tolower(val)
            
            # Exact name match
            hits <- lookup$kegg_id[lookup$name_lc == val_lc]
            if (length(hits) > 0) return(hits[1])
            
            # Starts-with match
            hits <- lookup$kegg_id[startsWith(lookup$name_lc, val_lc)]
            if (length(hits) > 0) return(hits[1])
            
            # Contains match (last resort)
            hits <- lookup$kegg_id[grepl(val_lc, lookup$name_lc, fixed = TRUE)]
            if (length(hits) > 0) return(hits[1])
            
            NA_character_
        }
        
        mapped <- vapply(x, map_one, character(1), USE.NAMES = FALSE)
        n_mapped   <- sum(!is.na(mapped) & mapped != "")
        n_unmapped <- sum(is.na(mapped) | mapped == "")
        
        if (n_unmapped > 0) {
            warning(n_unmapped, " entries could not be mapped to KEGG IDs and were dropped.",
                    call. = FALSE)
        }
        if (n_mapped == 0) {
            stop("None of the entered names could be mapped to KEGG IDs.", call. = FALSE)
        }
        
        out <- unique(mapped[!is.na(mapped) & mapped != ""])
        out <- out[grepl("^C[0-9]{5}$", out)]
        if (length(out) == 0) {
            stop("After mapping, no valid KEGG IDs remained.", call. = FALSE)
        }
        
        message("Mapped ", n_mapped, " entries → ", length(out), " unique KEGG IDs")
        out
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
            use_existing_significant <- "Significant" %in% colnames(kegg_data) && !force_custom_filters
            
            if (use_existing_significant) {
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
                
                if (fc_cutoff_up != 1 || fc_cutoff_down != -1 || fdr_cutoff_da != 0.05) {
                    message("Note: Custom FC/FDR filters provided but using existing Significant column. ",
                            "Set force_custom_filters = TRUE to apply custom filters.")
                }
                
            } else {
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
            metabolites_to_use <- kegg_data$kegg_id
            message("Using all ", length(metabolites_to_use), " metabolites from DA results")
        }
        
        metabolites_to_use <- metabolites_to_use[!is.na(metabolites_to_use)]
        metabolites_to_use <- metabolites_to_use[metabolites_to_use != ""]
        
        if (is.null(example_data) && "full_results" %in% names(da_results)) {
            example_data <- da_results$full_results
            message("Using da_results$full_results for GSEA analysis")
        }
        
        # Case 2: inputMetabolites provided directly
    } else if (!is.null(inputMetabolites)) {
        if (is.data.frame(inputMetabolites)) {
            if ("kegg_id" %in% colnames(inputMetabolites)) {
                metabolites_to_use <- inputMetabolites$kegg_id
                message("Using ", length(metabolites_to_use), " metabolites from inputMetabolites data frame")
                
                if (is.null(example_data) && all(c("pval", "log2fc") %in% colnames(inputMetabolites))) {
                    example_data <- inputMetabolites
                    message("Using inputMetabolites data frame for GSEA analysis")
                }
            } else if ("met_id" %in% colnames(inputMetabolites)) {
                metabolites_to_use <- inputMetabolites$met_id
                message("Using ", length(metabolites_to_use), " metabolites from inputMetabolites data frame (met_id column)")
            } else {
                stop("If inputMetabolites is a data frame, it must contain 'kegg_id' or 'met_id' column")
            }
        } else if (is.character(inputMetabolites)) {
            metabolites_to_use <- inputMetabolites
            message("Using ", length(metabolites_to_use), " metabolites from inputMetabolites character vector")
        } else {
            stop("inputMetabolites must be a character vector, data frame with 'kegg_id' or 'met_id' column, or NULL when da_results is provided")
        }
    } else {
        stop("Either inputMetabolites or da_results must be provided")
    }
    
    # ---- Name / mixed → KEGG (consistent with Shiny app) ----
    if (is.character(metabolites_to_use)) {
        if (input_type %in% c("name", "mixed")) {
            metabolites_to_use <- map_names_to_kegg(metabolites_to_use, kegg_lookup)
        } else {
            # input_type == "kegg": keep only valid C##### tokens
            keep <- grepl("^C[0-9]{5}$", metabolites_to_use, ignore.case = TRUE)
            if (sum(!keep) > 0) {
                warning(sum(!keep), " entries were not valid KEGG IDs (C#####) and were dropped.",
                        call. = FALSE)
            }
            metabolites_to_use <- unique(toupper(metabolites_to_use[keep]))
            if (length(metabolites_to_use) == 0) {
                stop("No valid KEGG IDs (C#####) found in inputMetabolites.", call. = FALSE)
            }
        }
    }
    
    # ---- Process complex KEGG IDs ----
    if (split_complex_ids && !is.null(metabolites_to_use)) {
        message("Processing complex KEGG IDs (splitting by |)...")
        original_count <- length(metabolites_to_use)
        all_split_ids <- unlist(lapply(metabolites_to_use, process_kegg_ids))
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
        message("Running pathway enrichment analysis",
                if (!is.null(backgroundMetabolites)) " (background-corrected)" else "",
                "...")
        
        if (!is.null(backgroundMetabolites)) {
            all_enrichment <- perform_enrichment_analysis_bg(
                inputMetabolites      = metabolites_to_use,
                PathwayVsMetabolites  = PathwayVsMetabolites,
                top_n                 = NULL,
                p_value_cutoff        = 1,
                backgroundMetabolites = backgroundMetabolites
            )
        } else {
            all_enrichment <- perform_enrichment_analysis(
                inputMetabolites     = metabolites_to_use,
                PathwayVsMetabolites = PathwayVsMetabolites,
                top_n                = NULL,
                p_value_cutoff       = 1
            )
        }
        
        results$pathway_enrichment_all <- all_enrichment
        
        significant_results <- all_enrichment %>%
            dplyr::filter(Adjusted_P_value < p_value_cutoff) %>%
            dplyr::arrange(dplyr::desc(Log_P_value))
        
        if (!is.null(top_n) && nrow(significant_results) > top_n) {
            significant_results <- head(significant_results, top_n)
        }
        
        results$pathway_enrichment_results <- significant_results
        
        if (run_plots && nrow(significant_results) > 0) {
            results$pathway_plot <- create_enrichment_plot(significant_results)
            results$impact_plot  <- create_impact_plot(significant_results)
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
    
    # ---- Reactome interaction ----
    if ("interaction" %in% analysis_type && run_plots) {
        message("Generating Reactome interaction network...")
        results$interaction_plot <- create_interaction_plot(
            inputMetabolites = metabolites_to_use,
            reactome_df = reactome_df,
            kegg_lookup = kegg_lookup
        )
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