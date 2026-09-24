# Internal: package-specific cache directory, with a tempdir fallback
.enrichmet_cache_dir <- function() {
    dir <- tools::R_user_dir("enrichmet", which = "cache")
    ok <- tryCatch({
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        file.access(dir, 2) == 0
    }, error = function(e) FALSE)
    if (!isTRUE(ok)) {
        dir <- file.path(tempdir(), "enrichmet_cache")
    }
    dir
}

#' Download and cache a remote file
#'
#' Downloads a remote resource and stores it in the EnrichMet cache.
#' Existing cached resources are reused automatically. Failed downloads
#' are retried, and a plain download is used if the cache fails.
#'
#' @param url Character. URL of the remote file.
#' @param retries Integer. Number of download attempts. Default is 3.
#'
#' @return Character path to the cached file.
#'
#' @examples
#' \donttest{
#' if (curl::has_internet()) {
#'     f <- get_cached_file("https://rest.kegg.jp/list/pathway/hsa")
#'     head(readLines(f))
#' }
#' }
#'
#' @export
get_cached_file <- function(url, retries = 3L) {
    if (!is.character(url) || length(url) != 1 ||
        is.na(url) || url == "") {
        stop("url must be a single non-empty character string.",
             call. = FALSE)
    }
    
    old <- options(timeout = max(300, getOption("timeout")))
    on.exit(options(old), add = TRUE)
    
    bfc <- tryCatch(
        BiocFileCache::BiocFileCache(
            cache = .enrichmet_cache_dir(), ask = FALSE
        ),
        error = function(e) NULL
    )
    
    # 1. Reuse a valid cache entry; drop stale ones (file missing)
    if (!is.null(bfc)) {
        cached <- BiocFileCache::bfcquery(
            bfc, query = url, field = "rname", exact = TRUE
        )
        for (rid in cached$rid) {
            path <- tryCatch(
                BiocFileCache::bfcrpath(bfc, rids = rid),
                error = function(e) NULL
            )
            if (length(path) > 0 && file.exists(path[1])) {
                return(unname(path[1]))
            }
        }
        if (nrow(cached) > 0) {
            tryCatch(
                BiocFileCache::bfcremove(bfc, cached$rid),
                error = function(e) NULL
            )
        }
    }
    
    # 2. Download, with retries
    for (i in seq_len(retries)) {
        if (!is.null(bfc)) {
            # bfcadd() returns the cached file path, named by the rid
            path <- tryCatch(
                BiocFileCache::bfcadd(
                    bfc, rname = url, fpath = url,
                    download = TRUE, rtype = "web"
                ),
                error = function(e) NULL
            )
            if (length(path) > 0 && file.exists(path[1])) {
                return(unname(path[1]))
            }
        }
        
        tmp <- tempfile()
        status <- tryCatch(
            utils::download.file(url, tmp, mode = "wb", quiet = TRUE),
            error = function(e) 1L,
            warning = function(w) 1L
        )
        if (identical(status, 0L) && file.exists(tmp) &&
            file.size(tmp) > 0) {
            return(tmp)
        }
        
        if (i < retries) {
            Sys.sleep(2^(i - 1))
        }
    }
    
    stop("Failed to download resource after ", retries,
         " attempts: ", url, call. = FALSE)
}
#' Fetch KEGG compound names
#'
#' Downloads the KEGG compound list and returns KEGG compound IDs
#' with their corresponding names.
#'
#' @param clean_names Logical. If TRUE, only the first compound name
#'   is retained when multiple names are provided.
#'
#' @return A data frame with columns \code{kegg_id} and \code{name}.
#'
#' 
#' kegg_lookup <- fetch_kegg_compound_lookup()
#' head(kegg_lookup)
#' 
#'
#' @export
fetch_kegg_compound_lookup <- function(
        clean_names = TRUE
) {
    
    url <- "https://rest.kegg.jp/list/compound"
    
    file_path <- get_cached_file(url)
    
    x <- readLines(
        file_path,
        warn = FALSE
    )
    
    if (length(x) == 0) {
        
        stop(
            "KEGG compound list is empty.",
            call. = FALSE
        )
    }
    
    lookup <- data.table::fread(
        text = paste(x, collapse = "\n"),
        sep = "\t",
        header = FALSE,
        col.names = c(
            "kegg_id",
            "name"
        ),
        fill = TRUE
    )
    
    lookup$kegg_id <- sub(
        "^cpd:",
        "",
        lookup$kegg_id
    )
    
    if (clean_names) {
        
        lookup$name <- sub(
            ";.*$",
            "",
            lookup$name
        )
        
        lookup$name <- trimws(
            lookup$name
        )
    }
    
    lookup <- lookup[
        !is.na(lookup$kegg_id) &
            lookup$kegg_id != "",
        ,
        drop = FALSE
    ]
    
    lookup <- unique(lookup)
    
    if (nrow(lookup) == 0) {
        
        stop(
            "No KEGG compound records were found.",
            call. = FALSE
        )
    }
    
    lookup
}

#' Fetch KEGG pathway-metabolite relationships
#'
#' Downloads KEGG pathway and compound relationship data and returns
#' pathway-level metabolite mappings via the public KEGG REST API.
#'
#' @param organism Character. KEGG organism code. Default is \code{"hsa"}.
#' @param clean_pathway_names Logical. If TRUE, KEGG organism information
#'   is removed from pathway names.
#'
#' @return A data frame with columns \code{PathwayID}, \code{Pathway},
#'   and \code{Metabolites}.
#'
#' @examples
#' 
#' PathwayVsMetabolites <- fetch_kegg_pathway_metabolites(organism = "hsa")
#' head(PathwayVsMetabolites)
#' 
#'
#' @export
fetch_kegg_pathway_metabolites <- function(
        organism = "hsa",
        clean_pathway_names = TRUE,
        metaboanalyst_style = TRUE
) {
    
    if (!is.character(organism) || length(organism) != 1 ||
        is.na(organism) || organism == "") {
        stop("organism must be a single non-empty character string.", call. = FALSE)
    }
    
    pathway_url <- paste0("https://rest.kegg.jp/list/pathway/", organism)
    link_url    <- "https://rest.kegg.jp/link/pathway/compound"
    
    pathway_file <- get_cached_file(pathway_url)
    link_file    <- get_cached_file(link_url)
    
    pathway_data <- read.delim(pathway_file, header = FALSE, sep = "\t",
                               stringsAsFactors = FALSE)
    
    if (ncol(pathway_data) < 2) {
        stop("Unexpected KEGG pathway file format.", call. = FALSE)
    }
    
    pathway_data <- pathway_data[, 1:2, drop = FALSE]
    colnames(pathway_data) <- c("PathwayID", "Pathway")
    pathway_data$PathwayID <- sub("^path:", "", pathway_data$PathwayID)
    
    # -------------------------------------------------------
    # Curation
    # -------------------------------------------------------
    if (isTRUE(metaboanalyst_style)) {
        
        pathway_num <- as.integer(
            sub(paste0("^", organism), "", pathway_data$PathwayID)
        )
        
        # Keep only metabolic pathways (numeric ID < 2000)
        # and exclude Global (011xx) + Overview (012xx) maps
        keep <- pathway_num < 2000 &
            !(pathway_num >= 1100 & pathway_num < 1300)
        
        n_excluded <- sum(!keep)
        pathway_data <- pathway_data[keep, , drop = FALSE]
        
        message(
            "Filter applied: ",
            "excluded ", n_excluded, " pathways. ",
            "Remaining: ", nrow(pathway_data)
        )
    }
    
    if (clean_pathway_names) {
        pathway_data$Pathway <- sub("^.*?:", "", pathway_data$Pathway)
        pathway_data$Pathway <- sub(" - [^-]+ \\([^)]*\\)$", "", pathway_data$Pathway)
        pathway_data$Pathway <- trimws(pathway_data$Pathway)
    }
    
    pathway_data$PathwayID_map <- sub(
        paste0("^", organism), "map", pathway_data$PathwayID
    )
    
    link_data <- read.delim(link_file, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE)
    
    if (ncol(link_data) < 2) {
        stop("Unexpected KEGG pathway-compound link file format.", call. = FALSE)
    }
    
    link_data <- link_data[, 1:2, drop = FALSE]
    colnames(link_data) <- c("Metabolite", "PathwayID_map")
    link_data$Metabolite    <- sub("^cpd:", "", link_data$Metabolite)
    link_data$PathwayID_map <- sub("^path:", "", link_data$PathwayID_map)
    
    merged_data <- merge(
        pathway_data[, c("PathwayID", "Pathway", "PathwayID_map"), drop = FALSE],
        link_data,
        by = "PathwayID_map"
    )
    
    if (nrow(merged_data) == 0) {
        stop("No KEGG pathway-metabolite relationships found for '", organism, "'.",
             call. = FALSE)
    }
    
    result <- aggregate(
        Metabolite ~ PathwayID + Pathway,
        data = merged_data,
        FUN = function(x) paste(unique(x), collapse = ",")
    )
    
    colnames(result)[colnames(result) == "Metabolite"] <- "Metabolites"
    result <- result[, c("PathwayID", "Pathway", "Metabolites"), drop = FALSE]
    rownames(result) <- NULL
    
    result
}
#' Extract KEGG compound IDs from pathway data
#'
#' Extracts unique KEGG compound IDs from the \code{Metabolites}
#' column of a KEGG pathway-metabolite table.
#'
#' @param PathwayVsMetabolites Data frame containing a
#'   \code{Metabolites} column with comma-separated KEGG compound IDs.
#'
#' @return A character vector containing unique KEGG compound IDs.
#'
#' @examples
#' 
#' PathwayVsMetabolites <- fetch_kegg_pathway_metabolites(organism = "hsa")
#' kegg_ids <- extract_kegg_ids(PathwayVsMetabolites)
#' head(kegg_ids)
#' 
#'
#' @export
extract_kegg_ids <- function(
        PathwayVsMetabolites
) {
    
    if (!is.data.frame(PathwayVsMetabolites)) {
        
        stop(
            "PathwayVsMetabolites must be a data frame.",
            call. = FALSE
        )
    }
    
    if (!"Metabolites" %in% colnames(PathwayVsMetabolites)) {
        
        stop(
            "PathwayVsMetabolites must contain a 'Metabolites' column.",
            call. = FALSE
        )
    }
    
    kegg_ids <- PathwayVsMetabolites %>%
        tidyr::separate_rows(
            Metabolites,
            sep = ","
        ) %>%
        dplyr::distinct(Metabolites) %>%
        dplyr::filter(
            !is.na(Metabolites),
            Metabolites != ""
        ) %>%
        dplyr::pull(Metabolites)
    
    if (length(kegg_ids) == 0) {
        
        stop(
            "No KEGG compound IDs were found in PathwayVsMetabolites.",
            call. = FALSE
        )
    }
    
    message(
        "Pathway universe contains ",
        length(kegg_ids),
        " unique KEGG compound IDs"
    )
    
    kegg_ids
}

#' Fetch KEGG-ChEBI mappings
#'
#' Downloads the KEGG-to-ChEBI conversion table and filters
#' it to the supplied KEGG compound IDs.
#'
#' @param kegg_ids Character vector of KEGG compound IDs.
#'
#' @return A data frame with columns \code{KEGG} and \code{ChEBI}.
#'
#' @examples
#' 
#' kegg_chebi <- fetch_kegg_chebi(c("C00022", "C00031"))
#' head(kegg_chebi)
#' 
#'
#' @export
fetch_kegg_chebi <- function(
        kegg_ids
) {
    
    if (!is.character(kegg_ids) ||
        length(kegg_ids) == 0) {
        
        stop(
            "kegg_ids must be a non-empty character vector.",
            call. = FALSE
        )
    }
    
    kegg_ids <- unique(
        trimws(
            as.character(kegg_ids)
        )
    )
    
    kegg_ids <- kegg_ids[
        grepl(
            "^C[0-9]+$",
            kegg_ids
        )
    ]
    
    if (length(kegg_ids) == 0) {
        
        stop(
            "No valid KEGG compound IDs were provided.",
            call. = FALSE
        )
    }
    
    url <- "https://rest.kegg.jp/conv/chebi/compound"
    
    file_path <- get_cached_file(url)
    
    x <- read.delim(
        file_path,
        header = FALSE,
        sep = "\t",
        stringsAsFactors = FALSE
    )
    
    if (ncol(x) < 2) {
        
        stop(
            "Unexpected KEGG-ChEBI conversion file format.",
            call. = FALSE
        )
    }
    
    col1 <- trimws(
        as.character(x[[1]])
    )
    
    col2 <- trimws(
        as.character(x[[2]])
    )
    
    is_kegg_1 <- grepl(
        "^cpd:",
        col1,
        ignore.case = TRUE
    )
    
    is_kegg_2 <- grepl(
        "^cpd:",
        col2,
        ignore.case = TRUE
    )
    
    is_chebi_1 <- grepl(
        "^chebi:",
        col1,
        ignore.case = TRUE
    )
    
    is_chebi_2 <- grepl(
        "^chebi:",
        col2,
        ignore.case = TRUE
    )
    
    if (any(is_kegg_1) &&
        any(is_chebi_2)) {
        
        kegg_col <- col1
        chebi_col <- col2
        
    } else if (
        any(is_kegg_2) &&
        any(is_chebi_1)
    ) {
        
        kegg_col <- col2
        chebi_col <- col1
        
    } else {
        
        stop(
            "Could not identify KEGG and ChEBI columns in the ",
            "KEGG conversion file.",
            call. = FALSE
        )
    }
    
    kegg_chebi <- data.frame(
        KEGG = sub(
            "^cpd:",
            "",
            kegg_col,
            ignore.case = TRUE
        ),
        ChEBI = sub(
            "^chebi:",
            "",
            chebi_col,
            ignore.case = TRUE
        ),
        stringsAsFactors = FALSE
    )
    
    kegg_chebi <- kegg_chebi[
        kegg_chebi$KEGG %in% kegg_ids,
        ,
        drop = FALSE
    ]
    
    kegg_chebi <- unique(
        kegg_chebi
    )
    
    if (nrow(kegg_chebi) == 0) {
        
        stop(
            "No KEGG-ChEBI mappings were found for the supplied ",
            "KEGG compound IDs.",
            call. = FALSE
        )
    }
    
    rownames(kegg_chebi) <- NULL
    
    kegg_chebi
}

#' Fetch Reactome reactions
#'
#' Downloads the Reactome ChEBI-to-reaction mapping and filters
#' it to the requested species. Data are retrieved at runtime and
#' cached; no Reactome dump is redistributed with the package.
#'
#' @param species Character. Species name used to filter Reactome data.
#'
#' @return A data frame containing ChEBI, Reactome, reaction,
#'   evidence, and species information.
#'
#' @examples
#' 
#' reactome <- fetch_reactome_reactions(species = "Homo sapiens")
#' head(reactome)
#' nrow(reactome)
#' 
#'
#' @export
fetch_reactome_reactions <- function(
        species = "Homo sapiens"
) {
    
    if (!is.character(species) ||
        length(species) != 1 ||
        is.na(species) ||
        species == "") {
        
        stop(
            "species must be a single non-empty character string.",
            call. = FALSE
        )
    }
    
    url <- paste0(
        "https://reactome.org/download/current/",
        "ChEBI2ReactomeReactions.txt"
    )
    
    file_path <- get_cached_file(
        url
    )
    
    reactome <- read.delim(
        file_path,
        header = FALSE,
        sep = "\t",
        stringsAsFactors = FALSE
    )
    
    if (ncol(reactome) < 6) {
        
        stop(
            "Unexpected Reactome ChEBI-to-reaction file format.",
            call. = FALSE
        )
    }
    
    reactome <- reactome[
        ,
        1:6,
        drop = FALSE
    ]
    
    colnames(reactome) <- c(
        "ChEBI",
        "Reactome",
        "Reactome_URL",
        "Reaction",
        "Evidence",
        "Species"
    )
    
    reactome$ChEBI <- sub(
        "^chebi:",
        "",
        reactome$ChEBI,
        ignore.case = TRUE
    )
    
    reactome$ChEBI <- trimws(
        reactome$ChEBI
    )
    
    reactome$Species <- trimws(
        reactome$Species
    )
    
    reactome <- reactome[
        reactome$Species == species,
        ,
        drop = FALSE
    ]
    
    reactome <- unique(
        reactome
    )
    
    if (nrow(reactome) == 0) {
        
        warning(
            "No Reactome reactions were found for species: ",
            species,
            call. = FALSE
        )
    }
    
    rownames(reactome) <- NULL
    
    reactome
}

#' Build KEGG-Reactome reaction mappings
#'
#' Links KEGG compounds to Reactome reactions through ChEBI identifiers.
#'
#' @param kegg_ids Character vector of KEGG compound IDs.
#' @param species Character. Reactome species name.
#'
#' @return A data frame containing KEGG-to-ChEBI-to-Reactome mappings.
#'
#' @examples
#' 
#' reactome_df <- fetch_kegg_reactome(
#'     kegg_ids = c("C00022", "C00031"),
#'     species = "Homo sapiens"
#' )
#' head(reactome_df)
#' 
#'
#' @export
fetch_kegg_reactome <- function(
        kegg_ids,
        species = "Homo sapiens"
) {
    
    if (!is.character(kegg_ids) ||
        length(kegg_ids) == 0) {
        
        stop(
            "kegg_ids must be a non-empty character vector.",
            call. = FALSE
        )
    }
    
    kegg_ids <- unique(
        trimws(
            as.character(kegg_ids)
        )
    )
    
    kegg_chebi <- fetch_kegg_chebi(
        kegg_ids = kegg_ids
    )
    
    reactome_reactions <- fetch_reactome_reactions(
        species = species
    )
    
    reactome_df <- merge(
        kegg_chebi,
        reactome_reactions,
        by = "ChEBI",
        all = FALSE
    )
    
    reactome_df <- unique(
        reactome_df
    )
    
    if (nrow(reactome_df) == 0) {
        
        warning(
            "No KEGG compounds could be mapped to Reactome reactions ",
            "for species: ",
            species,
            call. = FALSE
        )
    }
    
    rownames(reactome_df) <- NULL
    
    reactome_df
}

#' Fetch LION lipid ontology
#'
#' Downloads the LION ontology from BioPortal and converts the
#' OBO file into a pathway-metabolite table suitable for
#' enrichment analysis. Requires a BioPortal API key; no LION
#' data are redistributed with the package.
#'
#' @param api_key BioPortal API key. If NULL, the function uses
#'   the BIOPORTAL_API_KEY environment variable.
#' @param submission LION submission number. Default is 1.
#'
#' @return A data.frame with columns Pathway, Metabolites, and Count.
#'
#' @examples
#' # Expected output structure (no API call)
#' pathway_lipids <- data.frame(
#'     Pathway = c("Glycerophospholipids", "Sphingolipids"),
#'     Metabolites = c("PC(16:0/18:1),PE(18:0/18:1)", "SM(d18:1/16:0)"),
#'     Count = c(2L, 1L),
#'     stringsAsFactors = FALSE
#' )
#' head(pathway_lipids)
#'
#' 
#' # Live download when a BioPortal key is available
#' if (nzchar(Sys.getenv("BIOPORTAL_API_KEY"))) {
#'     pathway_lipids <- fetch_lion_lipid_ontology()
#'     head(pathway_lipids)
#' } else {
#'     message(
#'         "Set Sys.setenv(BIOPORTAL_API_KEY = \"your_key\") ",
#'         "to run the live LION download example."
#'     )
#' }
#' 
#'
#' @export
fetch_lion_lipid_ontology <- function(
        api_key = NULL,
        submission = 1
) {
    
    if (is.null(api_key) || api_key == "") {
        api_key <- Sys.getenv("BIOPORTAL_API_KEY")
    }
    
    if (is.null(api_key) || api_key == "") {
        stop(
            paste0(
                "A BioPortal API key is required to download the ",
                "LION ontology.\n",
                "Set it with:\n",
                "Sys.setenv(BIOPORTAL_API_KEY = \"your_api_key\")"
            ),
            call. = FALSE
        )
    }
    
    old_timeout <- getOption("timeout")
    options(timeout = max(300, old_timeout))
    
    on.exit(
        options(timeout = old_timeout),
        add = TRUE
    )
    
    url_obo <- paste0(
        "https://data.bioontology.org/ontologies/",
        "LION/submissions/",
        submission,
        "/download?apikey=",
        api_key
    )
    
    tmp_obo <- tempfile(
        fileext = ".obo"
    )
    
    on.exit(
        if (file.exists(tmp_obo)) unlink(tmp_obo),
        add = TRUE
    )
    
    tryCatch(
        {
            curl::curl_download(
                url_obo,
                destfile = tmp_obo,
                quiet = FALSE
            )
        },
        error = function(e) {
            stop(
                paste0(
                    "Failed to download the LION ontology from BioPortal.\n",
                    "Error: ",
                    conditionMessage(e)
                ),
                call. = FALSE
            )
        }
    )
    
    if (!file.exists(tmp_obo)) {
        stop(
            "LION ontology file was not created.",
            call. = FALSE
        )
    }
    
    file_size <- file.info(tmp_obo)$size
    
    if (is.na(file_size) || file_size == 0) {
        stop(
            "The downloaded LION ontology file is empty.",
            call. = FALSE
        )
    }
    
    lines <- readLines(
        tmp_obo,
        warn = FALSE
    )
    
    term_blocks <- split(
        lines,
        cumsum(
            grepl(
                "^\\[Term\\]",
                lines
            )
        )
    )
    
    parse_term <- function(block) {
        
        id_line <- grep(
            "^id:",
            block,
            value = TRUE
        )
        
        name_line <- grep(
            "^name:",
            block,
            value = TRUE
        )
        
        isa_line <- grep(
            "^is_a:",
            block,
            value = TRUE
        )
        
        if (
            length(id_line) == 0 ||
            length(name_line) == 0 ||
            length(isa_line) == 0
        ) {
            return(NULL)
        }
        
        id <- sub(
            "^id:\\s*",
            "",
            id_line
        )
        
        name <- sub(
            "^name:\\s*",
            "",
            name_line
        )
        
        ontology <- sub(
            ".*!\\s*",
            "",
            isa_line
        )
        
        ontology <- trimws(
            ontology
        )
        
        data.frame(
            Lipid_ID = id,
            Lipid_Name = name,
            Ontology_Classification = ontology,
            stringsAsFactors = FALSE
        )
    }
    
    parsed_list <- lapply(
        term_blocks,
        parse_term
    )
    
    parsed_list <- parsed_list[
        !vapply(
            parsed_list,
            is.null,
            logical(1)
        )
    ]
    
    if (length(parsed_list) == 0) {
        stop(
            "No valid LION terms could be parsed.",
            call. = FALSE
        )
    }
    
    parsed_terms <- do.call(
        rbind,
        parsed_list
    )
    
    parsed_terms <- unique(
        parsed_terms
    )
    
    parsed_terms <- parsed_terms[
        !is.na(parsed_terms$Lipid_Name),
        ,
        drop = FALSE
    ]
    
    parsed_terms_filtered <- parsed_terms[
        !grepl(
            "^CAT",
            parsed_terms$Lipid_ID
        ),
        ,
        drop = FALSE
    ]
    
    pathway_metabolites <- aggregate(
        Lipid_Name ~ Ontology_Classification,
        data = parsed_terms_filtered,
        FUN = function(x) {
            paste(
                unique(x),
                collapse = ","
            )
        }
    )
    
    pathway_counts <- aggregate(
        Lipid_Name ~ Ontology_Classification,
        data = parsed_terms_filtered,
        FUN = length
    )
    
    names(pathway_metabolites) <- c(
        "Pathway",
        "Metabolites"
    )
    
    names(pathway_counts) <- c(
        "Pathway",
        "Count"
    )
    
    pathway_metabolites <- merge(
        pathway_metabolites,
        pathway_counts,
        by = "Pathway",
        all.x = TRUE
    )
    
    pathway_metabolites <- pathway_metabolites[
        order(
            -pathway_metabolites$Count
        ),
        ,
        drop = FALSE
    ]
    
    rownames(pathway_metabolites) <- NULL
    
    pathway_metabolites
}