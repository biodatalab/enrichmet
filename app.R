# =============================================================================
# Metabolite Pathway Enrichment Explorer
# - KEGG IDs / metabolite names / mixed input (interaction plot fixed)
# - Interactive + static plots
# - Matching heatmap colours
# - Exact Reactome replica via reactome_exact.csv
# =============================================================================

options(shiny.maxRequestSize = 720 * 1024^2)

library(shiny)
library(shinyjs)
library(fgsea)
library(igraph)
library(ggraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(plotly)
library(openxlsx)
library(ComplexHeatmap)
library(pheatmap)
library(ggrepel)
library(scales)
library(stringr)
library(grid)
library(circlize)
library(data.table)
library(BiocFileCache)

heatmaply_available <- requireNamespace("heatmaply", quietly = TRUE)

# -----------------------------------------------------------------------------
# Load remaining analysis helpers (required)
# -----------------------------------------------------------------------------
# source("Perform_enrichment_analysis_background_corrected.R")
# source("enrichmet_helpers.R")  # create_* plot functions, etc.
# library(enrichmet)

# =============================================================================
# Inline: fetch_kegg_data helpers
# =============================================================================
.enrichmet_bfc <- function() {
    BiocFileCache::BiocFileCache(
        cache = tools::R_user_dir("enrichmet", which = "cache"),
        ask = FALSE
    )
}

get_cached_file <- function(url) {
    if (!is.character(url) || length(url) != 1 || is.na(url) || url == "") {
        stop("url must be a single non-empty character string.", call. = FALSE)
    }
    bfc <- .enrichmet_bfc()
    cached <- BiocFileCache::bfcquery(bfc, query = url, field = "rname")
    if (nrow(cached) > 0) {
        path <- BiocFileCache::bfcrpath(bfc, rids = cached$rid[1])
        if (length(path) > 0 && file.exists(path[1])) return(path[1])
    }
    rid <- BiocFileCache::bfcadd(bfc, rname = url, fpath = url, download = TRUE)
    path <- BiocFileCache::bfcrpath(bfc, rids = rid)
    if (length(path) == 0 || !file.exists(path[1])) {
        stop("Failed to download and cache resource: ", url, call. = FALSE)
    }
    path[1]
}

fetch_kegg_compound_lookup <- function(clean_names = TRUE) {
    url <- "https://rest.kegg.jp/list/compound"
    file_path <- get_cached_file(url)
    x <- readLines(file_path, warn = FALSE)
    if (length(x) == 0) stop("KEGG compound list is empty.", call. = FALSE)
    lookup <- data.table::fread(
        text = paste(x, collapse = "\n"), sep = "\t", header = FALSE,
        col.names = c("kegg_id", "name"), fill = TRUE
    )
    lookup$kegg_id <- sub("^cpd:", "", lookup$kegg_id)
    if (clean_names) {
        lookup$name <- sub(";.*$", "", lookup$name)
        lookup$name <- trimws(lookup$name)
    }
    lookup <- lookup[!is.na(lookup$kegg_id) & lookup$kegg_id != "", , drop = FALSE]
    unique(lookup)
}

fetch_kegg_pathway_metabolites <- function(organism = "hsa", clean_pathway_names = TRUE) {
    if (!is.character(organism) || length(organism) != 1 || is.na(organism) || organism == "") {
        stop("organism must be a single non-empty character string.", call. = FALSE)
    }
    pathway_url <- paste0("https://rest.kegg.jp/list/pathway/", organism)
    link_url <- "https://rest.kegg.jp/link/pathway/compound"
    pathway_file <- get_cached_file(pathway_url)
    link_file <- get_cached_file(link_url)
    
    pathway_data <- read.delim(pathway_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(pathway_data) < 2) stop("Unexpected KEGG pathway file format.", call. = FALSE)
    pathway_data <- pathway_data[, 1:2, drop = FALSE]
    colnames(pathway_data) <- c("PathwayID", "Pathway")
    pathway_data$PathwayID <- sub("^path:", "", pathway_data$PathwayID)
    if (clean_pathway_names) {
        pathway_data$Pathway <- sub("^.*?:", "", pathway_data$Pathway)
        pathway_data$Pathway <- sub(" - [^-]+ \\([^)]*\\)$", "", pathway_data$Pathway)
        pathway_data$Pathway <- trimws(pathway_data$Pathway)
    }
    pathway_data$PathwayID_map <- sub(paste0("^", organism), "map", pathway_data$PathwayID)
    
    link_data <- read.delim(link_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(link_data) < 2) stop("Unexpected KEGG pathway-compound link file format.", call. = FALSE)
    link_data <- link_data[, 1:2, drop = FALSE]
    colnames(link_data) <- c("Metabolite", "PathwayID_map")
    link_data$Metabolite <- sub("^cpd:", "", link_data$Metabolite)
    link_data$PathwayID_map <- sub("^path:", "", link_data$PathwayID_map)
    
    merged_data <- merge(
        pathway_data[, c("PathwayID", "Pathway", "PathwayID_map"), drop = FALSE],
        link_data, by = "PathwayID_map"
    )
    if (nrow(merged_data) == 0) {
        stop("No KEGG pathway-metabolite relationships found for '", organism, "'.", call. = FALSE)
    }
    result <- aggregate(Metabolite ~ PathwayID + Pathway, data = merged_data,
                        FUN = function(x) paste(unique(x), collapse = ","))
    colnames(result)[colnames(result) == "Metabolite"] <- "Metabolites"
    result <- result[, c("PathwayID", "Pathway", "Metabolites"), drop = FALSE]
    rownames(result) <- NULL
    result
}

fetch_kegg_chebi <- function(kegg_ids) {
    if (!is.character(kegg_ids) || length(kegg_ids) == 0) {
        stop("kegg_ids must be a non-empty character vector.", call. = FALSE)
    }
    kegg_ids <- unique(trimws(as.character(kegg_ids)))
    kegg_ids <- kegg_ids[grepl("^C[0-9]+$", kegg_ids)]
    if (length(kegg_ids) == 0) stop("No valid KEGG compound IDs were provided.", call. = FALSE)
    
    url <- "https://rest.kegg.jp/conv/chebi/compound"
    file_path <- get_cached_file(url)
    x <- read.delim(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(x) < 2) stop("Unexpected KEGG-ChEBI conversion file format.", call. = FALSE)
    
    col1 <- trimws(as.character(x[[1]]))
    col2 <- trimws(as.character(x[[2]]))
    is_kegg_1  <- grepl("^cpd:", col1, ignore.case = TRUE)
    is_kegg_2  <- grepl("^cpd:", col2, ignore.case = TRUE)
    is_chebi_1 <- grepl("^chebi:", col1, ignore.case = TRUE)
    is_chebi_2 <- grepl("^chebi:", col2, ignore.case = TRUE)
    
    if (any(is_kegg_1) && any(is_chebi_2)) {
        kegg_col <- col1; chebi_col <- col2
    } else if (any(is_kegg_2) && any(is_chebi_1)) {
        kegg_col <- col2; chebi_col <- col1
    } else {
        stop("Could not identify KEGG and ChEBI columns.", call. = FALSE)
    }
    
    kegg_chebi <- data.frame(
        KEGG  = sub("^cpd:", "", kegg_col, ignore.case = TRUE),
        ChEBI = sub("^chebi:", "", chebi_col, ignore.case = TRUE),
        stringsAsFactors = FALSE
    )
    kegg_chebi <- unique(kegg_chebi[kegg_chebi$KEGG %in% kegg_ids, , drop = FALSE])
    if (nrow(kegg_chebi) == 0) stop("No KEGG-ChEBI mappings found.", call. = FALSE)
    rownames(kegg_chebi) <- NULL
    kegg_chebi
}

fetch_reactome_reactions <- function(species = "Homo sapiens") {
    if (!is.character(species) || length(species) != 1 || is.na(species) || species == "") {
        stop("species must be a single non-empty character string.", call. = FALSE)
    }
    url <- "https://reactome.org/download/current/ChEBI2ReactomeReactions.txt"
    file_path <- get_cached_file(url)
    reactome <- read.delim(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (ncol(reactome) < 6) stop("Unexpected Reactome file format.", call. = FALSE)
    reactome <- reactome[, 1:6, drop = FALSE]
    colnames(reactome) <- c("ChEBI", "Reactome", "Reactome_URL", "Reaction", "Evidence", "Species")
    reactome$ChEBI   <- trimws(sub("^chebi:", "", reactome$ChEBI, ignore.case = TRUE))
    reactome$Species <- trimws(reactome$Species)
    reactome <- unique(reactome[reactome$Species == species, , drop = FALSE])
    if (nrow(reactome) == 0) warning("No Reactome reactions for species: ", species, call. = FALSE)
    rownames(reactome) <- NULL
    reactome
}

fetch_kegg_reactome <- function(kegg_ids, species = "Homo sapiens") {
    if (!is.character(kegg_ids) || length(kegg_ids) == 0) {
        stop("kegg_ids must be a non-empty character vector.", call. = FALSE)
    }
    kegg_ids <- unique(trimws(as.character(kegg_ids)))
    kegg_chebi <- fetch_kegg_chebi(kegg_ids = kegg_ids)
    reactome_reactions <- fetch_reactome_reactions(species = species)
    reactome_df <- unique(merge(kegg_chebi, reactome_reactions, by = "ChEBI", all = FALSE))
    if (nrow(reactome_df) == 0) {
        warning("No KEGG compounds mapped to Reactome for species: ", species, call. = FALSE)
    }
    rownames(reactome_df) <- NULL
    reactome_df
}

# =============================================================================
# create_interaction_plot (package version)
# =============================================================================
create_interaction_plot <- function(inputMetabolites, reactome_df, kegg_lookup = NULL,
                                    min_shared_reactions = 1, max_metabolites_per_reaction = 25) {
    if (is.null(reactome_df) || nrow(reactome_df) == 0) return(NULL)
    if (!all(c("KEGG", "Reaction") %in% colnames(reactome_df))) return(NULL)
    
    extract_kegg_id <- function(met_id) {
        if (is.na(met_id) || met_id == "") return(NA_character_)
        if (grepl("^C\\d{5}$", met_id)) return(met_id)
        if (grepl("_C\\d", met_id)) {
            parts <- unlist(strsplit(met_id, "_"))
            kegg_parts <- grep("^C\\d", parts, value = TRUE)
            if (length(kegg_parts) > 0) {
                if (grepl("\\|", kegg_parts[1])) {
                    all_kegg <- unlist(strsplit(kegg_parts[1], "\\|"))
                    valid_kegg <- grep("^C\\d", all_kegg, value = TRUE)
                    if (length(valid_kegg) > 0) return(valid_kegg[1])
                } else return(kegg_parts[1])
            }
        }
        return(NA_character_)
    }
    
    extracted <- vapply(inputMetabolites, extract_kegg_id, character(1))
    valid_kegg_ids <- unique(extracted[!is.na(extracted) & extracted != ""])
    if (length(valid_kegg_ids) == 0) return(NULL)
    
    reactome_sub <- reactome_df %>%
        dplyr::filter(KEGG %in% valid_kegg_ids) %>%
        dplyr::distinct(KEGG, Reaction, .keep_all = TRUE)
    if (nrow(reactome_sub) == 0) return(NULL)
    
    reactions_split <- split(reactome_sub$KEGG, reactome_sub$Reaction)
    edge_list <- lapply(reactions_split, function(mets) {
        mets <- unique(mets)
        if (length(mets) < 2 || length(mets) > max_metabolites_per_reaction) return(NULL)
        pairs <- utils::combn(sort(mets), 2, simplify = FALSE)
        do.call(rbind, lapply(pairs, function(p) data.frame(from = p[1], to = p[2], stringsAsFactors = FALSE)))
    })
    edge_list <- edge_list[!vapply(edge_list, is.null, logical(1))]
    if (length(edge_list) == 0) return(NULL)
    
    all_edges <- do.call(rbind, edge_list)
    valid_edges <- all_edges %>%
        dplyr::group_by(from, to) %>%
        dplyr::summarise(shared_reactions = dplyr::n(), .groups = "drop") %>%
        dplyr::filter(shared_reactions >= min_shared_reactions)
    if (nrow(valid_edges) == 0) return(NULL)
    
    vertex_ids <- unique(c(valid_edges$from, valid_edges$to))
    vertex_df  <- data.frame(KEGG_ID = vertex_ids, stringsAsFactors = FALSE)
    
    if (!is.null(kegg_lookup) && all(c("kegg_id", "name") %in% colnames(kegg_lookup))) {
        vertex_df <- vertex_df %>%
            dplyr::left_join(kegg_lookup, by = c("KEGG_ID" = "kegg_id")) %>%
            dplyr::mutate(display_name = ifelse(!is.na(name), name, KEGG_ID)) %>%
            dplyr::select(-name)
    } else {
        vertex_df <- vertex_df %>% dplyr::mutate(display_name = KEGG_ID)
    }
    vertex_df <- vertex_df %>%
        dplyr::mutate(display_name = stringr::str_trunc(display_name, 25)) %>%
        dplyr::select(name = KEGG_ID, display_name, dplyr::everything())
    
    g <- igraph::graph_from_data_frame(
        d = valid_edges %>% dplyr::mutate(weight = scales::rescale(shared_reactions, to = c(0.1, 1))),
        directed = FALSE, vertices = vertex_df
    )
    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) return(NULL)
    
    igraph::V(g)$degree      <- igraph::degree(g)
    igraph::V(g)$betweenness <- igraph::betweenness(g)
    igraph::V(g)$component   <- igraph::components(g)$membership
    
    degree_range  <- range(igraph::V(g)$degree)
    degree_breaks <- seq(degree_range[1], degree_range[2], by = 1)
    if (length(degree_breaks) > 5) degree_breaks <- pretty(degree_range, n = 4)
    
    comps <- igraph::components(g)
    if (comps$no > 1) {
        best_layout <- "gem"
    } else {
        layouts_to_try <- c("fr", "kk", "dh", "gem", "lgl")
        best_layout <- NULL; best_spacing <- 0
        for (ln in layouts_to_try) {
            tryCatch({
                pos <- ggraph::create_layout(g, layout = ln)
                d <- as.matrix(stats::dist(pos[, 1:2])); diag(d) <- NA
                sp <- mean(d, na.rm = TRUE)
                if (sp > best_spacing) { best_spacing <- sp; best_layout <- ln }
            }, error = function(e) NULL)
        }
        if (is.null(best_layout)) best_layout <- "fr"
    }
    
    ggraph::ggraph(g, layout = best_layout) +
        ggraph::geom_edge_link(aes(width = weight), color = "#606060", alpha = 0.4, show.legend = TRUE) +
        ggraph::geom_node_point(aes(size = degree, color = as.factor(component)), alpha = 0.8, stroke = 0.5) +
        ggraph::geom_node_text(aes(label = display_name), size = 3.5, repel = TRUE,
                               box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf,
                               min.segment.length = 0.2, segment.color = "grey30",
                               segment.alpha = 0.6, segment.size = 0.3,
                               family = "sans", fontface = "bold") +
        ggplot2::scale_color_discrete(name = "Network Component") +
        ggplot2::scale_size_continuous(name = "Degree (Connections)", range = c(3, 10),
                                       breaks = degree_breaks,
                                       guide = ggplot2::guide_legend(override.aes = list(color = "grey50"),
                                                                     nrow = min(4, length(degree_breaks)))) +
        ggraph::scale_edge_width(name = "Shared Reactions", range = c(0.5, 2)) +
        ggplot2::labs(title = "Metabolite Interaction Network (Reactome)",
                      subtitle = paste(igraph::vcount(g), "compounds with",
                                       igraph::ecount(g), "shared-reaction edges | Layout:",
                                       toupper(best_layout))) +
        ggraph::theme_graph(base_family = "sans") +
        ggplot2::theme(legend.position = "right",
                       plot.title = element_text(hjust = 0.5, face = "bold"),
                       plot.subtitle = element_text(hjust = 0.5),
                       legend.key.height = unit(0.8, "lines"))
}

build_reactome_graph <- function(inputMetabolites, reactome_df, kegg_lookup = NULL,
                                 min_shared_reactions = 1, max_metabolites_per_reaction = 25) {
    if (is.null(reactome_df) || nrow(reactome_df) == 0) return(list(error = "No reactome_df"))
    if (!all(c("KEGG", "Reaction") %in% colnames(reactome_df))) return(list(error = "Missing columns"))
    
    extract_kegg_id <- function(met_id) {
        if (is.na(met_id) || met_id == "") return(NA_character_)
        if (grepl("^C\\d{5}$", met_id)) return(met_id)
        if (grepl("_C\\d", met_id)) {
            parts <- unlist(strsplit(met_id, "_"))
            kegg_parts <- grep("^C\\d", parts, value = TRUE)
            if (length(kegg_parts) > 0) {
                if (grepl("\\|", kegg_parts[1])) {
                    all_kegg <- unlist(strsplit(kegg_parts[1], "\\|"))
                    valid_kegg <- grep("^C\\d", all_kegg, value = TRUE)
                    if (length(valid_kegg) > 0) return(valid_kegg[1])
                } else return(kegg_parts[1])
            }
        }
        return(NA_character_)
    }
    
    extracted <- vapply(inputMetabolites, extract_kegg_id, character(1))
    valid_kegg_ids <- unique(extracted[!is.na(extracted) & extracted != ""])
    if (length(valid_kegg_ids) == 0) return(list(error = "No valid KEGG IDs"))
    
    reactome_sub <- reactome_df %>%
        dplyr::filter(KEGG %in% valid_kegg_ids) %>%
        dplyr::distinct(KEGG, Reaction, .keep_all = TRUE)
    if (nrow(reactome_sub) == 0) return(list(error = "No reactions"))
    
    reactions_split <- split(reactome_sub$KEGG, reactome_sub$Reaction)
    edge_list <- lapply(reactions_split, function(mets) {
        mets <- unique(mets)
        if (length(mets) < 2 || length(mets) > max_metabolites_per_reaction) return(NULL)
        pairs <- utils::combn(sort(mets), 2, simplify = FALSE)
        do.call(rbind, lapply(pairs, function(p) data.frame(from = p[1], to = p[2], stringsAsFactors = FALSE)))
    })
    edge_list <- edge_list[!vapply(edge_list, is.null, logical(1))]
    if (length(edge_list) == 0) return(list(error = "No edges"))
    
    all_edges <- do.call(rbind, edge_list)
    valid_edges <- all_edges %>%
        dplyr::group_by(from, to) %>%
        dplyr::summarise(shared_reactions = dplyr::n(), .groups = "drop") %>%
        dplyr::filter(shared_reactions >= min_shared_reactions)
    if (nrow(valid_edges) == 0) return(list(error = "No valid edges"))
    
    vertex_ids <- unique(c(valid_edges$from, valid_edges$to))
    vertex_df  <- data.frame(KEGG_ID = vertex_ids, stringsAsFactors = FALSE)
    
    if (!is.null(kegg_lookup) && all(c("kegg_id", "name") %in% colnames(kegg_lookup))) {
        vertex_df <- vertex_df %>%
            dplyr::left_join(kegg_lookup, by = c("KEGG_ID" = "kegg_id")) %>%
            dplyr::mutate(display_name = ifelse(!is.na(name), name, KEGG_ID)) %>%
            dplyr::select(-name)
    } else {
        vertex_df <- vertex_df %>% dplyr::mutate(display_name = KEGG_ID)
    }
    vertex_df <- vertex_df %>%
        dplyr::mutate(display_name = stringr::str_trunc(display_name, 25)) %>%
        dplyr::select(name = KEGG_ID, display_name, dplyr::everything())
    
    g <- igraph::graph_from_data_frame(
        d = valid_edges %>% dplyr::mutate(weight = scales::rescale(shared_reactions, to = c(0.1, 1))),
        directed = FALSE, vertices = vertex_df
    )
    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) return(list(error = "Empty graph"))
    
    igraph::V(g)$degree      <- igraph::degree(g)
    igraph::V(g)$betweenness <- igraph::betweenness(g)
    igraph::V(g)$component   <- igraph::components(g)$membership
    
    comps <- igraph::components(g)
    best_layout <- if (comps$no > 1) "gem" else {
        layouts_to_try <- c("fr", "kk", "dh", "gem", "lgl")
        best <- "fr"; best_sp <- 0
        for (ln in layouts_to_try) {
            tryCatch({
                pos <- ggraph::create_layout(g, layout = ln)
                d <- as.matrix(stats::dist(pos[, 1:2])); diag(d) <- NA
                sp <- mean(d, na.rm = TRUE)
                if (sp > best_sp) { best_sp <- sp; best <- ln }
            }, error = function(e) NULL)
        }
        best
    }
    layout_obj <- ggraph::create_layout(g, layout = best_layout)
    list(error = NULL, g = g, layout_obj = layout_obj, valid_edges = valid_edges,
         best_layout = best_layout, n_nodes = igraph::vcount(g), n_edges = igraph::ecount(g))
}

render_reactome_plotly <- function(build) {
    nodes <- build$layout_obj
    edges <- build$valid_edges
    edge_coords <- edges %>%
        dplyr::left_join(nodes %>% dplyr::select(name, x, y), by = c("from" = "name")) %>%
        dplyr::rename(x0 = x, y0 = y) %>%
        dplyr::left_join(nodes %>% dplyr::select(name, x, y), by = c("to" = "name")) %>%
        dplyr::rename(x1 = x, y1 = y)
    edge_x <- as.vector(rbind(edge_coords$x0, edge_coords$x1, NA))
    edge_y <- as.vector(rbind(edge_coords$y0, edge_coords$y1, NA))
    
    p <- plotly::plot_ly()
    p <- plotly::add_trace(p, x = edge_x, y = edge_y, type = "scatter", mode = "lines",
                           line = list(color = "rgba(96,96,96,0.45)", width = 1.2),
                           hoverinfo = "none", showlegend = FALSE)
    p <- plotly::add_trace(p, data = nodes, x = ~x, y = ~y, type = "scatter", mode = "markers+text",
                           text = ~display_name, textposition = "top center",
                           marker = list(size = ~scales::rescale(degree, to = c(12, 28)),
                                         color = ~as.numeric(as.factor(component)),
                                         colorscale = "Set2", line = list(width = 1, color = "#333")),
                           hovertext = ~paste0("<b>", display_name, "</b><br>KEGG: ", name,
                                               "<br>Degree: ", degree, "<br>Component: ", component),
                           hoverinfo = "text", showlegend = FALSE)
    p <- plotly::layout(p,
                        title = list(text = paste0("Metabolite Interaction Network (Reactome)<br>",
                                                   "<sup>", build$n_nodes, " compounds • ", build$n_edges,
                                                   " edges | Layout: ", toupper(build$best_layout), "</sup>")),
                        xaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE),
                        yaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE),
                        hovermode = "closest")
    plotly::config(p, displaylogo = FALSE, modeBarButtonsToRemove = c("select2d", "lasso2d"))
}

heatmap_colors <- colorRampPalette(c("#F7F7F7", "#EF8A62", "#B2182B"))(100)

extract_heatmap_matrix <- function(ht) {
    tryCatch({
        if (inherits(ht, "HeatmapList")) ht <- ht@ht_list[[1]]
        if (inherits(ht, "Heatmap")) return(ht@matrix)
        NULL
    }, error = function(e) NULL)
}

# =============================================================================
# enrichmet() – requires helper functions to exist
# =============================================================================
enrichmet <- function(inputMetabolites = NULL, PathwayVsMetabolites, example_data = NULL,
                      da_results = NULL, top_n = 100, p_value_cutoff = 1, kegg_lookup = NULL,
                      reactome_df = NULL, output_dir = NULL, save_excel = FALSE,
                      analysis_type = c("enrichment", "gsea", "centrality", "network",
                                        "heatmap", "membership", "interaction"),
                      run_plots = TRUE, network_top_n = 20, heatmap_top_n = 20,
                      membership_top_n = 20, min_pathway_occurrence = 1,
                      min_metabolite_occurrence = 1, use_significant_only = TRUE,
                      split_complex_ids = TRUE, significance_threshold = "both",
                      include_volcano = TRUE, fc_cutoff_up = 1, fc_cutoff_down = -1,
                      fdr_cutoff_da = 0.05, force_custom_filters = FALSE,
                      backgroundMetabolites = NULL) {
    
    process_kegg_ids <- function(ids) {
        if (is.null(ids)) return(NULL)
        unique(unlist(strsplit(ids, "\\|"))[nzchar(unlist(strsplit(ids, "\\|")))])
    }
    
    metabolites_to_use <- NULL
    if (!is.null(da_results)) {
        if (!"kegg_ready" %in% names(da_results)) stop("da_results must contain 'kegg_ready'")
        kegg_data <- da_results$kegg_ready
        if (use_significant_only) {
            if ("Significant" %in% colnames(kegg_data) && !force_custom_filters) {
                if (significance_threshold == "up") {
                    metabolites_to_use <- kegg_data$kegg_id[kegg_data$Significant == "Up"]
                } else if (significance_threshold == "down") {
                    metabolites_to_use <- kegg_data$kegg_id[kegg_data$Significant == "Down"]
                } else {
                    metabolites_to_use <- kegg_data$kegg_id[kegg_data$Significant != "Not significant"]
                }
            } else {
                if (significance_threshold == "up") {
                    metabolites_to_use <- kegg_data$kegg_id[kegg_data$log2fc > fc_cutoff_up & kegg_data$padj < fdr_cutoff_da]
                } else if (significance_threshold == "down") {
                    metabolites_to_use <- kegg_data$kegg_id[kegg_data$log2fc < fc_cutoff_down & kegg_data$padj < fdr_cutoff_da]
                } else {
                    metabolites_to_use <- kegg_data$kegg_id[
                        (kegg_data$log2fc > fc_cutoff_up | kegg_data$log2fc < fc_cutoff_down) &
                            kegg_data$padj < fdr_cutoff_da]
                }
            }
        } else {
            metabolites_to_use <- kegg_data$kegg_id
        }
        metabolites_to_use <- metabolites_to_use[!is.na(metabolites_to_use) & metabolites_to_use != ""]
        if (is.null(example_data) && "full_results" %in% names(da_results)) {
            example_data <- da_results$full_results
        }
    } else if (!is.null(inputMetabolites)) {
        if (is.data.frame(inputMetabolites)) {
            if ("kegg_id" %in% colnames(inputMetabolites)) {
                metabolites_to_use <- inputMetabolites$kegg_id
            } else if ("met_id" %in% colnames(inputMetabolites)) {
                metabolites_to_use <- inputMetabolites$met_id
            } else stop("data frame must contain kegg_id or met_id")
        } else if (is.character(inputMetabolites)) {
            metabolites_to_use <- inputMetabolites
        } else stop("inputMetabolites must be character or data frame")
    } else stop("Either inputMetabolites or da_results must be provided")
    
    if (split_complex_ids && !is.null(metabolites_to_use)) {
        metabolites_to_use <- unique(unlist(lapply(metabolites_to_use, process_kegg_ids)))
    }
    if (length(metabolites_to_use) == 0) stop("No metabolites found for analysis.")
    if (!is.data.frame(PathwayVsMetabolites) ||
        !all(c("Pathway", "Metabolites") %in% colnames(PathwayVsMetabolites))) {
        stop("PathwayVsMetabolites must have Pathway and Metabolites columns")
    }
    
    results <- list(input_metabolites_used = metabolites_to_use)
    
    if ("enrichment" %in% analysis_type) {
        if (!is.null(backgroundMetabolites) && exists("perform_enrichment_analysis_bg")) {
            all_enrichment <- perform_enrichment_analysis_bg(
                inputMetabolites = metabolites_to_use,
                PathwayVsMetabolites = PathwayVsMetabolites,
                top_n = NULL, p_value_cutoff = 1,
                backgroundMetabolites = backgroundMetabolites)
        } else if (exists("perform_enrichment_analysis")) {
            all_enrichment <- perform_enrichment_analysis(
                inputMetabolites = metabolites_to_use,
                PathwayVsMetabolites = PathwayVsMetabolites,
                top_n = NULL, p_value_cutoff = 1)
        } else {
            stop("perform_enrichment_analysis() not found. Source helpers or library(enrichmet).")
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
            if (exists("create_enrichment_plot")) results$pathway_plot <- create_enrichment_plot(significant_results)
            if (exists("create_impact_plot"))     results$impact_plot  <- create_impact_plot(significant_results)
        }
    }
    
    if ("gsea" %in% analysis_type && !is.null(example_data) && exists("perform_gsea_analysis")) {
        gsea_results <- perform_gsea_analysis(example_data, PathwayVsMetabolites)
        results$gsea_results <- gsea_results
        if (run_plots && nrow(gsea_results) > 0 && exists("create_gsea_plot")) {
            results$gsea_plot <- create_gsea_plot(gsea_results, top_n = 20, kegg_lookup = kegg_lookup)
        }
    }
    
    if ("centrality" %in% analysis_type && exists("calculate_metabolite_centrality")) {
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
        if (run_plots && nrow(centrality_results) > 0 && exists("create_centrality_plot")) {
            results$rbc_plot <- create_centrality_plot(centrality_results, kegg_lookup = kegg_lookup)
        }
    }
    
    if ("network" %in% analysis_type && run_plots && exists("create_network_plot")) {
        results$network_plot <- create_network_plot(
            metabolites_to_use, PathwayVsMetabolites, kegg_lookup, top_n = network_top_n)
    }
    
    if ("heatmap" %in% analysis_type && run_plots &&
        !is.null(results$pathway_enrichment_results) && exists("create_heatmap_plot")) {
        results$heatmap_plot <- create_heatmap_plot(
            results$pathway_enrichment_results, PathwayVsMetabolites,
            metabolites_to_use, kegg_lookup,
            top_n = heatmap_top_n,
            min_pathways = min_pathway_occurrence,
            min_metabolites = min_metabolite_occurrence)
    }
    
    if ("membership" %in% analysis_type && run_plots && exists("create_membership_plot")) {
        results$membership_plot <- create_membership_plot(
            PathwayVsMetabolites, metabolites_to_use, kegg_lookup,
            top_n = membership_top_n,
            min_pathway_occurrence = min_pathway_occurrence,
            min_metabolite_occurrence = min_metabolite_occurrence)
    }
    
    if ("interaction" %in% analysis_type && run_plots) {
        results$interaction_plot <- create_interaction_plot(
            inputMetabolites = metabolites_to_use,
            reactome_df = reactome_df,
            kegg_lookup = kegg_lookup)
    }
    
    return(results)
}

# =============================================================================
# Species mapping
# =============================================================================
species_choices <- c(
    "Human (Homo sapiens)" = "hsa", "Mouse (Mus musculus)" = "mmu",
    "Rat (Rattus norvegicus)" = "rno", "Zebrafish (Danio rerio)" = "dre",
    "Fruit fly (Drosophila melanogaster)" = "dme", "Worm (Caenorhabditis elegans)" = "cel",
    "Yeast (Saccharomyces cerevisiae)" = "sce", "Arabidopsis (Arabidopsis thaliana)" = "ath",
    "Chicken (Gallus gallus)" = "gga", "Dog (Canis familiaris)" = "cfa",
    "Pig (Sus scrofa)" = "ssc", "Cow (Bos taurus)" = "bta"
)

reactome_species_map <- c(
    hsa = "Homo sapiens", mmu = "Mus musculus", rno = "Rattus norvegicus",
    dre = "Danio rerio", dme = "Drosophila melanogaster", cel = "Caenorhabditis elegans",
    sce = "Saccharomyces cerevisiae", ath = "Arabidopsis thaliana", gga = "Gallus gallus",
    cfa = "Canis familiaris", ssc = "Sus scrofa", bta = "Bos taurus"
)

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
    useShinyjs(),
    titlePanel("Metabolite Pathway Enrichment Explorer"),
    actionButton("toggleSidebar", "Toggle Sidebar", icon = icon("arrows-alt-h")),
    tags$hr(),
    fluidRow(
        column(width = 3, id = "sidebarPanel",
               wellPanel(
                   h4("1. Species"),
                   selectInput("species", "Model organism", choices = species_choices, selected = "hsa"),
                   
                   tags$hr(),
                   h4("2. Summary statistics"),
                   fileInput("summary_file", "Summary statistics data (CSV) – required for GSEA", accept = ".csv"),
                   
                   tags$hr(),
                   h4("3. Background correction (optional)"),
                   checkboxInput("use_background", "Restrict universe to measured metabolites", FALSE),
                   
                   tags$hr(),
                   h4("4. Input metabolites"),
                   radioButtons("input_type", "Input type",
                                choices = c("KEGG IDs" = "kegg",
                                            "Metabolite names" = "name",
                                            "Mixed (auto-detect)" = "mixed"),
                                selected = "kegg", inline = TRUE),
                   textAreaInput("metabolites",
                                 "Enter one ID or name per line (or comma-separated)",
                                 value = paste(c(
                                     "C00209","C00249","C06424","C01727","C04025","C00366","C00385",
                                     "C08261","C00079","C00506","C00093","C02979","C00106","C01384",
                                     "C00219","C00242","C19806","C00105","C00606","C00144","C01216",
                                     "C06429","C00475","C00386","C01046","C02180","C00082","C00299",
                                     "C00051","C16513","C00387","C00022","C00026","C01062","C00418",
                                     "C00295","C01762","C00525","C00074","C00491","C02067","C00015",
                                     "C14829","C00519","C06428","C05472","C00328","C16358","C16357",
                                     "C16353","C02632","C00246","C05842","C00785","C03739","C00212",
                                     "C01181","C00140","C05021","C02862","C02721","C01104","C10172",
                                     "C00380","C01262","C05122","C03736","C00599","C00446","C00275",
                                     "C00103","C00881","C15587","C00361","C00429","C01042","C11439",
                                     "C05637","C00127","C00035","C01909","C00836","C00847","C04294",
                                     "C00043","C01481","C00301","C05635","C02140","C05488","C00570",
                                     "C01551","C00092","C03139"
                                 ), collapse = "\n"), rows = 8),
                   helpText("Wait until pathways are ready, then use names or mixed. Matching is case-insensitive."),
                   
                   tags$hr(),
                   h4("5. Analyses"),
                   checkboxGroupInput("analysis_type", NULL,
                                      choices = c("enrichment","gsea","centrality","network",
                                                  "heatmap","membership","interaction"),
                                      selected = c("enrichment","gsea","centrality","network",
                                                   "heatmap","membership","interaction")),
                   numericInput("top_n", "Top N pathways", value = 20, min = 1),
                   numericInput("pval_cutoff", "Adjusted p-value cutoff", value = 0.05, step = 0.01),
                   
                   tags$hr(),
                   h4("6. Display options"),
                   radioButtons("plot_mode", "Plot mode",
                                choices = c("Static (matches package)" = "static",
                                            "Interactive (Plotly)" = "interactive"),
                                selected = "static", inline = TRUE),
                   sliderInput("plot_height", "Plot height (px)", min = 500, max = 1200, value = 800, step = 50),
                   
                   tags$hr(),
                   actionButton("run", "Run Analysis", class = "btn-primary btn-lg", width = "100%")
               )
        ),
        column(width = 9,
               tabsetPanel(
                   tabPanel("Pathway Plot",                uiOutput("pathway_ui")),
                   tabPanel("Impact Plot",                 uiOutput("impact_ui")),
                   tabPanel("Metabolite Set Enrichment",   uiOutput("gsea_ui")),
                   tabPanel("Metabolite Centrality",       uiOutput("rbc_ui")),
                   tabPanel("Membership Plot",             uiOutput("membership_ui")),
                   tabPanel("Network Plot",                uiOutput("network_ui")),
                   tabPanel("Interaction Plot (Reactome)", uiOutput("interaction_ui")),
                   tabPanel("Heatmap",                     uiOutput("heatmap_ui"))
               )
        )
    )
)

# =============================================================================
# Server
# =============================================================================
server <- function(input, output, session) {
    
    observeEvent(input$toggleSidebar, {
        runjs("$('#sidebarPanel').toggle(); if ($('#sidebarPanel').is(':visible')) { $('.col-sm-9').removeClass('col-sm-12').addClass('col-sm-9'); } else { $('.col-sm-9').removeClass('col-sm-9').addClass('col-sm-12'); }")
    })
    
    bg <- reactiveValues(PathwayVsMetabolites = NULL, kegg_lookup = NULL, current_species = NULL)
    reactome_build <- reactiveVal(NULL)
    
    observeEvent(input$species, {
        org <- input$species
        if (identical(org, bg$current_species) && !is.null(bg$PathwayVsMetabolites)) return()
        
        showNotification(paste("Loading data for", names(species_choices)[species_choices == org], "…"),
                         type = "message", duration = 5)
        
        withProgress(message = "Fetching KEGG pathways…", value = 0.3, {
            bg$PathwayVsMetabolites <- tryCatch(
                fetch_kegg_pathway_metabolites(organism = org, clean_pathway_names = TRUE),
                error = function(e) {
                    showNotification(paste("Pathway fetch failed:", e$message), type = "error")
                    NULL
                }
            )
        })
        
        if (is.null(bg$kegg_lookup)) {
            withProgress(message = "Fetching KEGG compound names…", value = 0.6, {
                bg$kegg_lookup <- tryCatch(
                    fetch_kegg_compound_lookup(clean_names = TRUE),
                    error = function(e) {
                        showNotification(paste("Compound lookup failed:", e$message), type = "warning")
                        NULL
                    }
                )
            })
        }
        
        bg$current_species <- org
        if (!is.null(bg$PathwayVsMetabolites)) {
            showNotification(paste(nrow(bg$PathwayVsMetabolites), "pathways ready"), type = "message", duration = 3)
        }
    }, ignoreNULL = FALSE)
    
    results <- eventReactive(input$run, {
        req(bg$PathwayVsMetabolites)
        req(input$summary_file)
        
        example_data <- read.csv(input$summary_file$datapath, stringsAsFactors = FALSE)
        if (!"pval" %in% names(example_data) && "padj" %in% names(example_data)) {
            example_data$pval <- example_data$padj
        }
        
        # ---- Parse input → always clean KEGG IDs ----
        raw_input <- strsplit(input$metabolites, "[\n,]")[[1]] %>%
            trimws() %>% gsub('["\']', "", .) %>% .[. != ""] %>% unique()
        
        if (length(raw_input) == 0) {
            showNotification("Please enter at least one metabolite", type = "warning")
            return(NULL)
        }
        
        inputMetabolites <- raw_input
        
        if (input$input_type %in% c("name", "mixed")) {
            
            if (is.null(bg$kegg_lookup) || nrow(bg$kegg_lookup) == 0) {
                showNotification(
                    "KEGG compound names are still loading (or failed). Wait for 'pathways ready' then try again.",
                    type = "error", duration = 8
                )
                return(NULL)
            }
            
            lookup <- bg$kegg_lookup
            lookup$name_lc <- tolower(trimws(as.character(lookup$name)))
            lookup$kegg_id <- trimws(as.character(lookup$kegg_id))
            
            map_one <- function(x) {
                x <- trimws(x)
                if (x == "") return(NA_character_)
                if (grepl("^C[0-9]{5}$", x, ignore.case = TRUE)) return(toupper(x))
                x_lc <- tolower(x)
                hits <- lookup$kegg_id[lookup$name_lc == x_lc]
                if (length(hits) > 0) return(hits[1])
                hits <- lookup$kegg_id[startsWith(lookup$name_lc, x_lc)]
                if (length(hits) > 0) return(hits[1])
                hits <- lookup$kegg_id[grepl(x_lc, lookup$name_lc, fixed = TRUE)]
                if (length(hits) > 0) return(hits[1])
                NA_character_
            }
            
            mapped <- vapply(raw_input, map_one, character(1), USE.NAMES = FALSE)
            n_mapped   <- sum(!is.na(mapped) & mapped != "")
            n_unmapped <- sum(is.na(mapped) | mapped == "")
            
            if (n_unmapped > 0) {
                showNotification(
                    paste(n_unmapped, "entries could not be mapped to KEGG IDs and were dropped."),
                    type = "warning", duration = 6
                )
            }
            if (n_mapped == 0) {
                showNotification("None of the entered names could be mapped to KEGG IDs.", type = "error")
                return(NULL)
            }
            
            inputMetabolites <- unique(mapped[!is.na(mapped) & mapped != ""])
            inputMetabolites <- inputMetabolites[grepl("^C[0-9]{5}$", inputMetabolites)]
            
            if (length(inputMetabolites) == 0) {
                showNotification("After mapping, no valid KEGG IDs remained.", type = "error")
                return(NULL)
            }
            
            showNotification(
                paste("Mapped", n_mapped, "→", length(inputMetabolites), "unique KEGG IDs"),
                type = "message", duration = 4
            )
            
        } else {
            inputMetabolites <- unique(raw_input[grepl("^C[0-9]{5}$", raw_input, ignore.case = TRUE)])
            inputMetabolites <- toupper(inputMetabolites)
            if (length(inputMetabolites) == 0) {
                showNotification("No valid KEGG IDs (C#####) found in the input.", type = "error")
                return(NULL)
            }
        }
        
        # ---- Reactome (always uses cleaned KEGG IDs) ----
        reactome_df <- NULL
        if ("interaction" %in% input$analysis_type) {
            exact_file <- "reactome_exact.csv"
            if (file.exists(exact_file)) {
                reactome_df <- read.csv(exact_file, stringsAsFactors = FALSE)
                showNotification(
                    paste("Using exact Reactome file →", length(unique(reactome_df$KEGG)), "metabolites"),
                    type = "message"
                )
            } else {
                reactome_sp <- reactome_species_map[[input$species]]
                if (is.null(reactome_sp)) reactome_sp <- "Homo sapiens"
                withProgress(message = paste("Fetching Reactome (", reactome_sp, ")…"), value = 0.5, {
                    reactome_df <- tryCatch(
                        fetch_kegg_reactome(kegg_ids = inputMetabolites, species = reactome_sp),
                        error = function(e) {
                            showNotification(paste("Reactome fetch failed:", e$message), type = "warning")
                            NULL
                        }
                    )
                })
                if (!is.null(reactome_df)) {
                    showNotification(
                        paste("Live Reactome →", length(unique(reactome_df$KEGG)), "metabolites"),
                        type = "message"
                    )
                }
            }
            
            if (!is.null(reactome_df) && nrow(reactome_df) > 0) {
                build <- build_reactome_graph(
                    inputMetabolites = inputMetabolites,
                    reactome_df      = reactome_df,
                    kegg_lookup      = bg$kegg_lookup
                )
                reactome_build(build)
                if (!is.null(build$error)) {
                    showNotification(paste("Interaction graph:", build$error), type = "warning")
                }
            } else {
                reactome_build(NULL)
                showNotification("No Reactome data available for interaction plot.", type = "warning")
            }
        }
        
        background_ids <- NULL
        if (isTRUE(input$use_background) && exists("extract_measured_kegg_ids")) {
            background_ids <- tryCatch(extract_measured_kegg_ids(example_data), error = function(e) NULL)
        }
        
        tryCatch({
            enrichmet(
                inputMetabolites          = inputMetabolites,
                PathwayVsMetabolites      = bg$PathwayVsMetabolites,
                example_data              = example_data,
                kegg_lookup               = bg$kegg_lookup,
                reactome_df               = reactome_df,
                top_n                     = input$top_n,
                p_value_cutoff            = input$pval_cutoff,
                analysis_type             = input$analysis_type,
                network_top_n             = input$top_n,
                heatmap_top_n             = input$top_n,
                membership_top_n          = input$top_n,
                min_pathway_occurrence    = 1,
                min_metabolite_occurrence = 1,
                backgroundMetabolites     = background_ids
            )
        }, error = function(e) {
            showNotification(paste("Analysis error:", e$message), type = "error")
            NULL
        })
    })
    
    ph <- reactive(paste0(input$plot_height, "px"))
    
    make_gg <- function(get_plot, placeholder) {
        list(
            static = renderPlot({
                req(results())
                p <- get_plot()
                if (!is.null(p)) p else ggplot() + annotate("text", 0.5, 0.5, label = placeholder, size = 5) + theme_void()
            }),
            plotly = plotly::renderPlotly({
                req(results())
                p <- get_plot()
                if (is.null(p)) return(plotly::plotly_empty())
                plotly::ggplotly(p) %>% plotly::config(displaylogo = FALSE)
            })
        )
    }
    
    for (nm in c("pathway", "impact", "gsea", "rbc")) {
        local({
            id <- nm
            getter <- switch(id,
                             pathway = function() results()$pathway_plot,
                             impact  = function() results()$impact_plot,
                             gsea    = function() results()$gsea_plot,
                             rbc     = function() results()$rbc_plot
            )
            outs <- make_gg(getter, paste("No", id, "plot"))
            output[[paste0(id, "_static")]] <<- outs$static
            output[[paste0(id, "_plotly")]] <<- outs$plotly
            output[[paste0(id, "_ui")]] <<- renderUI({
                if (identical(input$plot_mode, "interactive")) {
                    plotly::plotlyOutput(paste0(id, "_plotly"), height = ph())
                } else {
                    plotOutput(paste0(id, "_static"), height = ph())
                }
            })
        })
    }
    
    output$interaction_static <- renderPlot({
        req(results())
        if (!is.null(results()$interaction_plot)) results()$interaction_plot
        else ggplot() + annotate("text", 0.5, 0.5, label = "No interaction plot", size = 5) + theme_void()
    })
    output$interaction_plotly <- plotly::renderPlotly({
        req(results())
        build <- reactome_build()
        if (is.null(build) || !is.null(build$error)) return(plotly::plotly_empty())
        render_reactome_plotly(build)
    })
    output$interaction_ui <- renderUI({
        if (identical(input$plot_mode, "interactive")) {
            plotly::plotlyOutput("interaction_plotly", height = ph())
        } else {
            plotOutput("interaction_static", height = ph())
        }
    })
    
    output$network_static <- renderPlot({
        req(results())
        if (!is.null(results()$network_plot)) results()$network_plot
        else ggplot() + annotate("text", 0.5, 0.5, label = "No network plot", size = 5) + theme_void()
    })
    output$network_ui <- renderUI(plotOutput("network_static", height = ph()))
    
    make_ht <- function(get_ht, placeholder) {
        list(
            static = renderPlot({
                req(results())
                ht <- get_ht()
                if (!is.null(ht)) { grid::grid.newpage(); ComplexHeatmap::draw(ht) }
                else ggplot() + annotate("text", 0.5, 0.5, label = placeholder, size = 5) + theme_void()
            }),
            plotly = plotly::renderPlotly({
                req(results())
                ht  <- get_ht()
                mat <- if (!is.null(ht)) extract_heatmap_matrix(ht) else NULL
                if (is.null(mat) || !heatmaply_available) return(plotly::plotly_empty())
                heatmaply::heatmaply(mat, colors = heatmap_colors, scale = "none",
                                     showticklabels = c(TRUE, TRUE),
                                     fontsize_row = 10, fontsize_col = 10) %>%
                    plotly::config(displaylogo = FALSE)
            })
        )
    }
    
    for (nm in c("membership", "heatmap")) {
        local({
            id <- nm
            getter <- if (id == "membership") function() results()$membership_plot else function() results()$heatmap_plot
            outs <- make_ht(getter, paste("No", id, "plot"))
            output[[paste0(id, "_static")]] <<- outs$static
            output[[paste0(id, "_plotly")]] <<- outs$plotly
            output[[paste0(id, "_ui")]] <<- renderUI({
                ht <- getter()
                mat_ok <- !is.null(ht) && !is.null(extract_heatmap_matrix(ht)) && heatmaply_available
                if (identical(input$plot_mode, "interactive") && mat_ok) {
                    plotly::plotlyOutput(paste0(id, "_plotly"), height = ph())
                } else {
                    plotOutput(paste0(id, "_static"), height = ph())
                }
            })
        })
    }
}

shinyApp(ui = ui, server = server)