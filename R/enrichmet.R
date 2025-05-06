#' ENRICHMET: R package for Quick and Easy Pathway Analysis in Metabolomics
#'
#' This function performs pathway enrichment analysis using Fisher’s exact test, impact plot, computes betweenness centrality for metabolites,
#' and performs Metabolite Set Enrichment Analysis (MetSEA).
#' It also generates plots for pathway enrichment, MetSEA,relative betweenness centrality (RBC), network plot, heatmap plot, membership plot, interaction plot.
#'
#' @param inputMetabolites A vector of metabolites for which pathway enrichment and centrality analysis are to be performed.
#' @param PathwayVsMetabolites A data frame containing pathways and their associated metabolites.
#' @param example_data A data frame containing example data for GSEA. This should include columns such as "met_id", "pval", and "log2fc".
#' @param top_n An integer specifying the number of top pathways to include in the pathway enrichment results (default is 100).
#' @param p_value_cutoff A numeric value for adjusting the p-value threshold for filtering significant pathways (default is 1).
#' @param mapping_df A data frame containing the mapping of metabolite IDs (KEGG_ID) to their corresponding STITCH IDs.
#' @param stitch_df A data frame containing STITCH interaction data, including STITCH IDs and chemical–chemical interaction information downloaded from the STITCH database.
#' @return A list containing three ggplot objects: pathway enrichment plot, GSEA plot, and RBC plot.
#' @examples
#' ## ** Examples
#' library(igraph)
#' library(ggraph)
#' library(tidygraph)
#' library(ggplot2)
#' library(openxlsx)
#' library(dplyr)
#' library(tidyr)
#' library(fgsea)
#' library(purrr)
#' library(ggrepel)
#' library(stringr)
#' library(pheatmap)
#' library(ComplexHeatmap)
#' library(ggplotify)
#' library(enrichmet)
#' # Generate example data with at least n=50 metabolites
#' set.seed(1234)
#'
#' # Create 50 unique metabolites
#' inputMetabolites <- paste0("M", 1:20)
#'
#' # Generate 10 pathways with random metabolites assigned
#' pathway_names <- paste0("Pathway", 1:50)
#' PathwayVsMetabolites <- data.frame(
#'   Pathway = rep(pathway_names, each = 1),
#'   Metabolites = sapply(1:50, function(x) paste(sample(inputMetabolites, sample(5:15, 1)), collapse = ","))
#' )
#'
#' # Add new pathway entries (Pathway101 and Pathway102)
#' new_rows <- data.frame(
#'   Pathway = c("Pathway101", "Pathway102", "Pathway103", "Pathway104", "pathway105"),
#'   Metabolites = c(
#'     "M12,M13,M14,M15,M16,M1,M18,M3,M29,M6,M16,M4",
#'     "M6,M7,M8,M9,M10,M11,M9,M29,M6,M6,M16,M4",
#'     "M24,M25,M26,M27,M28,M29,M30,M29,M26,M5",
#'     "M13,M14,M15,M16,M17,M24,M27,M14",
#'     "M15,M16,M17,M18,M19,M20,M21,M4,M8,M10"
#'   )
#' )
#'
#' # Combine with existing PathwayVsMetabolites
#' PathwayVsMetabolites <- rbind(PathwayVsMetabolites, new_rows)
#'
#' # Generate example metabolite-level data
#' example_data <- data.frame(
#'   met_id = inputMetabolites,
#'   pval = runif(20, 0.001, 0.05),  # Random p-values between 0.001 and 0.05
#'   log2fc = rnorm(20, mean = 0, sd = 1)  # Log2 fold changes from normal distribution
#' )
#'
#'# ---- 4. Create mapping_df ----
#' set.seed(42)
#' mapping_df <- data.frame(
#'  KEGG_ID = inputMetabolites,
#'  PubChem_CID = as.character(sample(10000:99999, length(inputMetabolites))),
#' STITCH_ID = paste0("CIDs", str_pad(sample(1000:9999, length(inputMetabolites)), 8, pad = "0"))
#')
#'
#' # ---- 5. Create synthetic STITCH interaction data ----
#' stitch_ids <- mapping_df$STITCH_ID
#'
#' stitch_pairs <- expand.grid(chemical1 = stitch_ids, chemical2 = stitch_ids) %>%
#'  filter(chemical1 != chemical2)
#'
#' set.seed(123)
#' stitch_df <- stitch_pairs %>%
#' slice_sample(n = 200) %>%
#' mutate(
#'   similarity = runif(n(), 0, 1),
#'   experimental = sample(0:500, n(), replace = TRUE),
#'   database = sample(c(0, 300, 600, 900), n(), replace = TRUE),
#'   textmining = sample(0:1000, n(), replace = TRUE),
#'   combined_score = similarity * 200 + experimental + database + textmining
#'  ) %>%
#'  as_tibble()
#'
# ---- 6. Run enrichment analysis ----
#' enrichmet(
#' inputMetabolites = inputMetabolites,
#' PathwayVsMetabolites = PathwayVsMetabolites, example_data,
#' top_n = 20,
#'  mapping_df = mapping_df,
#' stitch_df = stitch_df
#' )
#' @export
enrichmet <- function(inputMetabolites, PathwayVsMetabolites, example_data, top_n = 100, p_value_cutoff = 1, kegg_lookup = NULL, mapping_df = NULL, stitch_df = NULL) {

  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(purrr)
  library(ggrepel)
  library(stringr)
  library(pheatmap)
  library(ComplexHeatmap)
  library(ggplotify)

  matrix_to_list <- function(pws) {
    pws.l <- list()
    for (pw in colnames(pws)) {
      pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
    }
    return(pws.l)
  }

  prepare_gmt <- function(gmt_file, metabolites_in_data, savefile = FALSE) {
    gmt <- gmtPathways(gmt_file)
    hidden <- unique(unlist(gmt))
    mat <- matrix(NA, dimnames = list(hidden, names(gmt)), nrow = length(hidden), ncol = length(gmt))
    for (i in 1:dim(mat)[2]) {
      mat[, i] <- as.numeric(hidden %in% gmt[[i]])
    }
    hidden1 <- intersect(metabolites_in_data, hidden)
    mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1, ]) > 5)]]
    final_list <- matrix_to_list(mat)
    if (savefile) {
      saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
    }
    return(final_list)
  }

  PathwayVsMetabolites$description <- "https://www.genome.jp/kegg/pathway.html#metabolism"

  convert_to_gmt <- function(pathway, description, metabolites) {
    pathway_underscore <- gsub(" ", "_", pathway)
    metabolites_vector <- unlist(strsplit(metabolites, ","))
    gmt_line <- paste(pathway_underscore, description, paste(metabolites_vector, collapse = "\t"), sep = "\t")
    return(gmt_line)
  }

  gmt_data <- mapply(convert_to_gmt, PathwayVsMetabolites$Pathway, PathwayVsMetabolites$description, PathwayVsMetabolites$Metabolites, SIMPLIFY = TRUE)
  gmt_file <- "output.gmt"
  writeLines(gmt_data, gmt_file)

  data <- PathwayVsMetabolites %>%
    mutate(Metabolites = strsplit(Metabolites, ",")) %>%
    unnest(Metabolites)

  allMetabolitesSet <- unique(data$Metabolites)

  # Centrality
  edge_list_metabolites <- data.frame(from = unlist(data$Metabolites), to = rep(data$Pathway, lengths(data$Metabolites)))
  g_metabolites <- graph_from_data_frame(d = edge_list_metabolites, directed = FALSE)
  betweenness_metabolites <- betweenness(g_metabolites, directed = FALSE, normalized = TRUE)
  metabolite_centrality <- data.frame(Metabolite = names(betweenness_metabolites), RBC_Metabolite = betweenness_metabolites)

  input_metabolite_centrality <- metabolite_centrality %>%
    filter(Metabolite %in% inputMetabolites) %>%
    arrange(desc(RBC_Metabolite))

  # Pathway enrichment
  results <- list()
  for (i in 1:nrow(PathwayVsMetabolites)) {
    row <- PathwayVsMetabolites[i, ]
    pathway <- row$Pathway
    pathwayMetabolites <- unlist(strsplit(row$Metabolites, ","))
    matchedMet <- intersect(pathwayMetabolites, inputMetabolites)

    a <- length(matchedMet)
    b <- length(setdiff(inputMetabolites, pathwayMetabolites))
    c <- length(setdiff(pathwayMetabolites, inputMetabolites))
    d <- length(setdiff(allMetabolitesSet, union(inputMetabolites, pathwayMetabolites)))

    contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    fisher_test_result <- fisher.test(contingency_table, alternative = "two.sided")

    matched_centrality <- metabolite_centrality %>%
      filter(Metabolite %in% matchedMet)

    all_centrality <- metabolite_centrality %>%
      filter(Metabolite %in% pathwayMetabolites)

    impact <- ifelse(nrow(all_centrality) > 0, sum(matched_centrality$RBC_Metabolite) / sum(all_centrality$RBC_Metabolite), 0)
    coverage <- length(matchedMet) / length(pathwayMetabolites)

    results[[i]] <- list(
      Pathway = pathway,
      P_value = fisher_test_result$p.value,
      Log_P_value = -log10(fisher_test_result$p.value),
      Impact = impact,
      Coverage = coverage
    )
  }

  results_df <- do.call(rbind, lapply(results, as.data.frame))
  results_df$Adjusted_P_value <- p.adjust(results_df$P_value, method = "BH")

  significant_results_df <- results_df %>%
    filter(Adjusted_P_value < p_value_cutoff) %>%
    arrange(desc(Log_P_value))

  if (!is.null(top_n)) {
    significant_results_df <- head(significant_results_df, top_n)
  }

  if (nrow(significant_results_df) > 0) {
    # Updated Pathway Plot: Bubble plot sorted by -logPvalue
    pathway_plot <- ggplot(significant_results_df, aes(x = reorder(Pathway, Log_P_value), y = Log_P_value, size = Log_P_value, fill = Adjusted_P_value)) +
      geom_point(shape = 21, color = "NA", alpha = 0.7) +  # Bubble plot
      scale_size_continuous(range = c(3, 10), guide = "none") +  # Hide size legend
      scale_fill_gradient(
        low = "red",
        high = "blue",
        breaks = c(0.05, 0.5, 1),  # Specify the breaks
        labels = c("0", "0.05", "1")  # Set the labels
      ) +  # Color based on adjusted p-value
      labs(x = "Pathway", y = "-log10(P-value)", fill = "Adj P-value", size = "-log10(P-value)") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for better readability
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      )


    impact_plot <- ggplot(significant_results_df, aes(x = Impact, y = Log_P_value, size = Impact, color = P_value, label = Pathway)) +
      geom_point(alpha = 0.85) +
      geom_text_repel(size = 3, max.overlaps = 100, box.padding = 0.3, show.legend = FALSE) +
      scale_size_continuous(range = c(3, 10)) +
      scale_color_gradient(low = "red", high = "blue", name = "P-value") +
      labs(
        title = "Pathway Impact vs Significance",
        x = "Pathway Impact",
        y = "-log10(P-value)",
        size = "Impact"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
  } else {
    warning("No pathways passed the p-value cutoff.")
    pathway_plot <- NULL
    impact_plot <- NULL
  }

  # GSEA Plot
  example_filtered_data <- example_data %>%
    arrange(pval) %>%
    distinct(met_id, .keep_all = TRUE) %>%
    filter(met_id != "No Metabolites found")

  meta <- example_filtered_data$met_id
  bg_metabolites <- prepare_gmt(gmt_file, meta, savefile = FALSE)
  rankings <- sign(example_filtered_data$log2fc) * (-log10(as.numeric(example_filtered_data$pval)))
  names(rankings) <- example_filtered_data$met_id
  rankings <- sort(rankings, decreasing = TRUE)

  MSEAres <- fgsea(pathways = bg_metabolites, stats = rankings, scoreType = 'std', minSize = 10, maxSize = 500, nproc = 1)
  MSEAres$input_count <- sapply(MSEAres$leadingEdge, length)

  # Sort the GSEA results by p-value (ascending order)
  MSEAres <- MSEAres %>%
    arrange(pval)  # Sort by p-value in ascending order

  # Reversing the factor levels to make sure the lowest p-values (most significant) are at the top
  MSEAres$pathway <- factor(MSEAres$pathway, levels = rev(MSEAres$pathway))

  # Now create the gsea_plot
  gsea_plot <- ggplot(MSEAres, aes(x = -log10(pval), y = pathway, size = input_count, color = NES)) +
    geom_point() +
    labs(x = "-log10(p-value)", y = "Pathway", size = "Metabolite count", color = "NES") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal()

  write.xlsx(significant_results_df, "pathway_enrichment_results.xlsx")
  write.xlsx(MSEAres, "gsea_results.xlsx")

  # Centrality Plot (RBC Plot)
  if (!is.null(kegg_lookup)) {
    input_metabolite_centrality <- input_metabolite_centrality %>%
      left_join(kegg_lookup, by = c("Metabolite" = "kegg_id")) %>%
      mutate(Metabolite = coalesce(name, Metabolite))
  }

  rbc_plot <- ggplot(input_metabolite_centrality, aes(x = reorder(Metabolite, RBC_Metabolite), y = RBC_Metabolite, fill = RBC_Metabolite)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "orange", high = "red3") +
    labs(x = "Metabolite", y = "Relative Betweenness Centrality", fill = "RBC") +
    theme_minimal()

  # Network Graph (Metabolite-Pathway Network)
  df <- PathwayVsMetabolites %>%
    mutate(Metabolite = strsplit(Metabolites, ",")) %>%
    unnest(Metabolite) %>%
    filter(Metabolite %in% inputMetabolites)

  edge_list_all <- PathwayVsMetabolites %>%
    mutate(Metabolite = strsplit(Metabolites, ",")) %>%
    unnest(Metabolite) %>%
    select(Pathway, Metabolite)

  g_all <- graph_from_data_frame(edge_list_all, directed = FALSE)
  bet_all <- betweenness(g_all, directed = FALSE, normalized = TRUE)
  bet_df <- tibble(name = names(bet_all), centrality = as.numeric(bet_all))

  edges <- df %>% select(Pathway, Metabolite)
  g <- graph_from_data_frame(edges, directed = FALSE)

  V(g)$type <- ifelse(V(g)$name %in% df$Pathway, "Pathway", "Metabolite")
  V(g)$centrality <- NA
  V(g)$centrality[V(g)$type == "Metabolite"] <- bet_df$centrality[match(V(g)$name[V(g)$type == "Metabolite"], bet_df$name)]
  V(g)$group <- ifelse(V(g)$type == "Pathway", V(g)$name, df$Pathway[match(V(g)$name, df$Metabolite)])

  if (!is.null(kegg_lookup)) {
    V(g)$name <- ifelse(V(g)$type == "Metabolite",
                        kegg_lookup$name[match(V(g)$name, kegg_lookup$kegg_id)],
                        V(g)$name)
  }

  network_layout <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(edge_alpha = 0.5), show.legend = FALSE) +
    geom_node_point(aes(color = centrality, size = centrality)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(3, 8)) +
    theme_minimal() +
    labs(title = "Metabolite-Pathway Network", color = "Centrality")

  network_plot <- network_layout + theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

  # Enrichment Heatmap
  enriched_pathways <- significant_results_df$Pathway
  data_filtered <- data %>%
    filter(Pathway %in% enriched_pathways, Metabolites %in% inputMetabolites)

  # Convert KEGG IDs to names for display if kegg_lookup is provided
  if (!is.null(kegg_lookup)) {
    data_filtered <- data_filtered %>%
      left_join(kegg_lookup, by = c("Metabolites" = "kegg_id")) %>%
      mutate(Metabolites = coalesce(name, Metabolites))
  }

  heatmap_matrix <- table(data_filtered$Metabolites, data_filtered$Pathway)
  heatmap_matrix <- as.matrix(heatmap_matrix)

  logp_vec <- significant_results_df$Log_P_value
  names(logp_vec) <- significant_results_df$Pathway

  heatmap_values <- sweep(heatmap_matrix, 2, logp_vec[colnames(heatmap_matrix)], `*`)
  # Render heatmap with modified legend
  # Calculate appropriate breaks
  heatmap_plot <- Heatmap(
    heatmap_values,
    name = "-logP-value",
    col = circlize::colorRamp2(
      c(0, 3, 6),
      c("white", "blue", "red")
    ),
    cluster_columns = TRUE,
    border=TRUE,
    column_names_gp = gpar(fontsize = 9),  # X-axis label size
    row_names_gp = gpar(fontsize = 9)
  )
  # === Pathway Membership Plot ===
  if (!is.null(PathwayVsMetabolites)) {
    long_pathway_df <- PathwayVsMetabolites %>%
      tidyr::separate_rows(Metabolites, sep = ",") %>%
      dplyr::rename(pathway_name = Pathway, metabolite = Metabolites)

    matched_df <- long_pathway_df %>%
      dplyr::filter(metabolite %in% inputMetabolites)

    if (nrow(matched_df) > 0) {
      # Convert KEGG IDs to names if kegg_lookup is provided
      if (!is.null(kegg_lookup)) {
        matched_df <- matched_df %>%
          dplyr::left_join(kegg_lookup, by = c("metabolite" = "kegg_id")) %>%
          dplyr::mutate(name = ifelse(is.na(name), metabolite, name))
      } else {
        matched_df <- matched_df %>% dplyr::mutate(name = metabolite)
      }

      membership_matrix <- table(matched_df$name, matched_df$pathway_name)
      membership_matrix <- as.matrix(membership_matrix)
      membership_matrix[membership_matrix > 1] <- 1

      n_rows <- nrow(membership_matrix)
      n_cols <- ncol(membership_matrix)

      # Dynamically set size based on number of rows/cols
      heatmap_width <- unit(0.3 * n_cols, "cm")
      heatmap_height <- unit(0.4 * n_rows, "cm")

      membership_plot <- ComplexHeatmap::Heatmap(
        membership_matrix,
        name = "Membership",
        col = c("0" = "white", "1" = "steelblue"),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_title = "Metabolites",
        column_title = "Pathways",
        heatmap_legend_param = list(title = "Member"),
        row_names_gp = grid::gpar(fontsize = 8),
        column_names_gp = grid::gpar(fontsize = 8),
        width = heatmap_width,
        height = heatmap_height,
        border = TRUE
      )
    } else {
      warning("No matching input metabolites found in PathwayVsMetabolites.")
      membership_plot <- NULL
    }
  } else {
    warning("PathwayVsMetabolites not provided.")
    membership_plot <- NULL
  }
  # Interaction plot
  interaction_plot <- NULL  # Default if no data provided

  # Optional interaction plot if mapping_df and stitch_df are provided
  if (!is.null(mapping_df) && !is.null(stitch_df)) {

    # Filter mapping_df to only inputMetabolites (KEGG IDs)
    vertex_df <- mapping_df %>%
      filter(KEGG_ID %in% inputMetabolites) %>%
      filter(!is.na(PubChem_CID)) %>%
      distinct(KEGG_ID, .keep_all = TRUE)

    # Merge with KEGG lookup (if available) to get names
    if (!is.null(kegg_lookup)) {
      vertex_df <- vertex_df %>%
        left_join(kegg_lookup, by = c("KEGG_ID" = "kegg_id")) %>%
        mutate(display_name = coalesce(name, KEGG_ID))
    } else {
      vertex_df <- vertex_df %>%
        mutate(display_name = KEGG_ID)
    }

    vertex_df <- vertex_df %>%
      mutate(display_name = stringr::str_trunc(display_name, 25)) %>%
      select(STITCH_ID, everything())

    # Filter STITCH edges to valid vertices
    valid_edges <- stitch_df %>%
      filter(combined_score >= 100) %>%
      semi_join(vertex_df, by = c("chemical1" = "STITCH_ID")) %>%
      semi_join(vertex_df, by = c("chemical2" = "STITCH_ID")) %>%
      distinct(chemical1, chemical2, .keep_all = TRUE)

    # Build graph
    g <- igraph::graph_from_data_frame(
      d = valid_edges %>% mutate(weight = scales::rescale(combined_score, to = c(0.1, 1))),
      directed = FALSE,
      vertices = vertex_df %>%
        filter(STITCH_ID %in% c(valid_edges$chemical1, valid_edges$chemical2))
    )

    # Graph metrics
    igraph::V(g)$degree <- igraph::degree(g)
    igraph::V(g)$betweenness <- igraph::betweenness(g)
    igraph::V(g)$component <- igraph::components(g)$membership

    # Plot
    set.seed(42)
    interaction_plot <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(aes(alpha = weight), color = "grey50") +
      ggraph::geom_node_point(aes(size = degree, color = as.factor(component)), alpha = 0.8) +
      ggraph::geom_node_text(aes(label = stringr::str_wrap(display_name, 12)),
                             size = 3, repel = TRUE, max.overlaps = 20) +
      scale_color_discrete(name = "Network Component") +
      scale_size_continuous(name = "Degree", range = c(3, 10)) +
      labs(title = "Metabolite Interaction Network",
           subtitle = paste(igraph::vcount(g), "compounds with", igraph::ecount(g), "interactions")) +
      ggraph::theme_graph() +
      theme(legend.position = "right",
            plot.title = element_text(hjust = 0.5, face = "bold"))
  }

  # Return everything
  return(list(
    pathway_plot = pathway_plot,
    impact_plot = impact_plot,
    gsea_plot = gsea_plot,
    rbc_plot = rbc_plot,
    network_plot = network_plot,
    heatmap_plot = heatmap_plot,
    membership_plot = membership_plot,
    interaction_plot=interaction_plot
  ))
}


