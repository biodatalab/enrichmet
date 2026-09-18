# R/helpers.R

#' Load example datasets for enrichmet
#'
#' This function loads minimal example datasets for demonstrating
#' the functionality of enrichmet package functions.
#'
#' @param dataset Character string specifying which dataset to load.
#'   Options are: "pathway", "summary", "kegg", "mapping", "stitch", or "all".
#' @return The requested dataset or a list of datasets if "all" is specified.
#' @examples
#' # Load pathway data
#' pathway_data <- load_example_data("pathway")
#' 
#' # Load summary statistics
#' summary_data <- load_example_data("summary")
#' 
#' # Load all example data
#' all_data <- load_example_data("all")
#' @export
load_example_data <- function(dataset = c("pathway", "summary", "kegg", "mapping", "stitch", "all")) {
    dataset <- match.arg(dataset)
    
    # Path to the example data file
    example_file <- system.file("extdata", "enrichmet_example.rda", 
                                package = "enrichmet")
    
    if (!file.exists(example_file)) {
        # Try alternative location during development
        dev_file <- "inst/extdata/enrichmet_example.rda"
        if (file.exists(dev_file)) {
            example_file <- dev_file
        } else {
            stop("Example data file not found. Expected at: ", 
                 system.file("extdata", package = "enrichmet"), "/enrichmet_example.rda\n",
                 "Please reinstall the package or run create_example_data.R")
        }
    }
    
    # Load the data
    env <- new.env()
    suppressWarnings({
        loaded <- tryCatch({
            load(example_file, envir = env)
            TRUE
        }, error = function(e) {
            FALSE
        })
    })
    
    if (!loaded || length(ls(env)) == 0) {
        stop("Failed to load example data. File may be corrupted.")
    }
    
    # Define expected objects
    expected_objects <- c("example_pathway", "example_summary", "example_kegg", 
                          "example_mapping", "example_stitch")
    
    # Check which objects are available
    available <- expected_objects[expected_objects %in% ls(env)]
    
    if (dataset == "all") {
        result <- list()
        for (obj in available) {
            result[[gsub("example_", "", obj)]] <- env[[obj]]
        }
        return(result)
    } else {
        # Map dataset name to object name
        obj_name <- switch(dataset,
                           pathway = "example_pathway",
                           summary = "example_summary",
                           kegg = "example_kegg",
                           mapping = "example_mapping",
                           stitch = "example_stitch"
        )
        
        if (!obj_name %in% ls(env)) {
            stop("Requested dataset '", dataset, "' not found in the data file.")
        }
        
        return(env[[obj_name]])
    }
}


#' Get a cached file path using BiocFileCache
#'
#' Downloads and caches a remote file if it is not already present.
#' Designed to work reliably on clean Bioconductor build machines.
#'
#' @param url Character. Single URL to download/cache.
#' @return Character. Local path to the cached file.
#' @keywords internal
get_cached_file <- function(url) {
    
    if (!is.character(url) || length(url) != 1 || is.na(url) || url == "") {
        stop("url must be a single non-empty character string.", call. = FALSE)
    }
    
    bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
    
    # 1. Try to find existing entry
    cached <- BiocFileCache::bfcquery(bfc, query = url, field = "rname")
    
    if (nrow(cached) > 0) {
        rid <- cached$rid[1]
        path <- tryCatch(
            BiocFileCache::bfcrpath(bfc, rids = rid),
            error = function(e) NULL
        )
        if (!is.null(path) && length(path) > 0 && file.exists(path[1])) {
            return(path[1])
        }
    }
    
    # 2. Download + add
    rid <- tryCatch(
        BiocFileCache::bfcadd(
            bfc,
            rname = url,
            fpath = url,
            download = TRUE,
            rtype = "web"
        ),
        error = function(e) {
            stop("Failed to download resource: ", url, "\n", conditionMessage(e), call. = FALSE)
        }
    )
    
    # 3. Safely retrieve the path (retry once if needed)
    path <- tryCatch(
        BiocFileCache::bfcrpath(bfc, rids = rid),
        error = function(e) {
            # Rare race condition on some build machines – re-query
            cached2 <- BiocFileCache::bfcquery(bfc, query = url, field = "rname")
            if (nrow(cached2) > 0) {
                BiocFileCache::bfcrpath(bfc, rids = cached2$rid[1])
            } else {
                stop("Failed to retrieve cached path for: ", url, call. = FALSE)
            }
        }
    )
    
    if (length(path) == 0 || !file.exists(path[1])) {
        stop("Cached file does not exist after download: ", url, call. = FALSE)
    }
    
    path[1]
}
# Optional: Create a test function
#' @keywords internal
.test_example_data <- function() {
    cat("Testing example data...\n")
    tryCatch({
        data <- load_example_data("all")
        cat("✓ Successfully loaded all data\n")
        cat("  Available datasets:", paste(names(data), collapse = ", "), "\n")
        return(TRUE)
    }, error = function(e) {
        cat("✗ Error:", e$message, "\n")
        return(FALSE)
    })
}