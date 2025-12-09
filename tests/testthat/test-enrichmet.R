library(readr)
library(readxl)
library(httr)
library(BiocFileCache)

# InputMetabolites: A vector or list containing the metabolites of interest.
inputMetabolites <- c(
    "C00209", "C00249", "C06424", "C01727", "C04025", "C00366", "C00385", "C08261", "C00079", "C00506", "C00093", "C02979", "C00106", "C01384", "C00219", "C00242", "C19806", "C00105", "C00606", "C00144", "C01216", "C06429",
    "C00475", "C00386", "C01046", "C02180", "C00082", "C00299", "C00051", "C16513", "C00387", "C00022", "C00026",
    "C01062", "C00418", "C00295", "C01762", "C00525", "C00074", "C00491", "C02067", "C00015", "C14829", "C00519",
    "C06428", "C05472", "C00328", "C16358", "C16357", "C16353", "C02632", "C00246", "C05842", "C00785", "C03739",
    "C00212", "C01181", "C00140", "C05021", "C02862", "C02721", "C01104", "C10172", "C00380", "C01262", "C05122",
    "C03736", "C00599", "C00446", "C00275", "C00103", "C00881", "C15587", "C00361", "C00429", "C01042", "C11439",
    "C05637", "C00127", "C00035", "C01909", "C00836", "C00847", "C04294", "C00043", "C01481", "C00301", "C05635",
    "C02140", "C05488", "C00570", "C01551", "C00092", "C03139"
)
# Standard cache location
cache_dir <- tools::R_user_dir("enrichmet", "cache")
bfc <- BiocFileCache(cache_dir, ask = FALSE)

# Helper function to download & cache any file from Zenodo
get_cached_file <- function(url) {
    bfcrpath(bfc, url)
}

# -------------------------
# Load files with caching
# -------------------------

# PathwayVsMetabolites
pathway_url <- "https://zenodo.org/api/records/17819145/files/human_pathway.csv/content"
pathway_path <- get_cached_file(pathway_url)
PathwayVsMetabolites <- read.csv(pathway_path)
str(PathwayVsMetabolites)

# SummaryStatistics
example_url <- "https://zenodo.org/api/records/17819145/files/summary_stat.csv/content"
example_path <- get_cached_file(example_url)
example_data <- read.csv(example_path)
head(example_data)

# KEGGLookup
kegg_url <- "https://zenodo.org/api/records/17819145/files/kegg_lookup.xlsx/content"
kegg_path <- get_cached_file(kegg_url)
kegg_lookup <- read_excel(kegg_path)
head(kegg_lookup)


# Mapping_df
mapping_url <- "https://zenodo.org/api/records/17819145/files/mapping_df.xlsx/content"
mapping_path <- get_cached_file(mapping_url)
mapping_df <- read_excel(mapping_path)
head(mapping_df)

# STITCHInteractions
stitch_url <- "https://zenodo.org/records/17819145/files/stitch.tsv"
stitch_path <- get_cached_file(stitch_url)
stitch_df <- read_tsv(stitch_path)
head(stitch_df)

test_that("enrichmet returns expected output format", {
result <- enrichmet(
    inputMetabolites = inputMetabolites,
    PathwayVsMetabolites = PathwayVsMetabolites,
    example_data = example_data,  # Named parameter
    kegg_lookup = kegg_lookup,    # Named parameter
    top_n = 20,
    p_value_cutoff = 0.05,
    mapping_df = mapping_df,
    stitch_df = stitch_df,
    network_top_n = 10,
    heatmap_top_n = 20,
    membership_top_n = 20,
    min_pathway_occurrence = 2,
    min_metabolite_occurrence = 1
)
  expect_type(result, "list")
  expect_named(result, c("input_metabolites_used",
                         "pathway_enrichment_all",
                         "pathway_enrichment_results",
                         "pathway_plot",
                         "impact_plot",
                         "gsea_results",
                         "gsea_plot",
                         "metabolite_centrality",
                         "rbc_plot",
                         "network_plot",
                         "heatmap_plot",
                         "membership_plot",
                         "interaction_plot"))
})
