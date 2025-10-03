library(readr)
library(readxl)
library(httr)
library(BiocFileCache)

# InputMetabolites: A vector or list containing the metabolites of interest.
inputMetabolites <- c(
    "C02721", "C03736", "C03139", "C00074", "C02862", "C00092", "C00275",
    "C00386", "C10172", "C00361", "C00242", "C19806", "C16358", "C16357",
    "C16353", "C02140", "C05488", "C01181", "C03017", "C02632", "C00246",
    "C00079", "C08317", "C00588", "C16527", "C00065", "C00249", "C00836",
    "C02990", "C03690", "C04025", "C06424", "C00123", "C00606", "C03067",
    "C00539", "C00633", "C00519", "C00366", "C00148", "C00219", "C08261",
    "C00385", "C06429", "C05637", "C00295", "C06428", "C16513", "C01909",
    "C03299", "C01530", "C01571", "C08278", "C00599", "C15587", "C00144",
    "C03739"
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
pathway_url <- "https://zenodo.org/api/records/15498097/files/Human-PathwaysVsMetabolites.csv/content"
pathway_path <- get_cached_file(pathway_url)
PathwayVsMetabolites <- read.csv(pathway_path)
str(PathwayVsMetabolites)

# SummaryStatistics
example_url <- "https://zenodo.org/api/records/15498097/files/example_data.xlsx/content"
example_path <- get_cached_file(example_url)
example_data <- read_excel(example_path)
head(example_data)

# KEGGLookup
kegg_url <- "https://zenodo.org/api/records/15498097/files/kegg_lookup.xlsx/content"
kegg_path <- get_cached_file(kegg_url)
kegg_lookup <- read_excel(kegg_path)
head(kegg_lookup)

# Mapping_df
mapping_url <- "https://zenodo.org/api/records/15498097/files/mapping_df.xlsx/content"
mapping_path <- get_cached_file(mapping_url)
mapping_df <- read_excel(mapping_path)
head(mapping_df)

# STITCHInteractions
stitch_url <- "https://zenodo.org/records/15498097/files/chemical_chemical.tsv"
stitch_path <- get_cached_file(stitch_url)
stitch_df <- read_tsv(stitch_path)
head(stitch_df)

test_that("enrichmet returns expected output format", {
  result <- enrichmet(inputMetabolites, PathwayVsMetabolites=PathwayVsMetabolites, example_data, top_n = 100, p_value_cutoff = 1, kegg_lookup = kegg_lookup, mapping_df = mapping_df, stitch_df = stitch_df)
  expect_type(result, "list")
  expect_named(result, c("pathway_enrichment_results",
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
