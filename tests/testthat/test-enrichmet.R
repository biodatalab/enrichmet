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
# PathwayVsMetabolites: A file mapping pathways to their corresponding metabolites.
PathwayVsMetabolites=read.csv("https://zenodo.org/api/records/15498097/files/Human-PathwaysVsMetabolites.csv/content")
str(PathwayVsMetabolites)
# SummaryStatistics: A file containing differential analysis results or other relevant statistical outputs.
example_data=read_excel("https://zenodo.org/api/records/15498097/files/example_data.xlsx/content")
head(example_data)
# KEGGLookup: A file containing annotated KEGG compound IDs with their corresponding metabolite names.
kegg_lookup=read_excel("https://zenodo.org/api/records/15498097/files/kegg_lookup.xlsx/content")
head(kegg_lookup)
# Mapping_df: A data frame that maps metabolites to their corresponding STITCH IDs.
mapping_df = read_excel("https://zenodo.org/api/records/15498097/files/mapping_df.xlsx/content")
head(mapping_df)
#STITCHInteractions: A file obtained from the STITCH database containing chemical–chemical interaction information.
# Use a standard cache location (BiocFileCache will create it if missing)
cache_dir <- tools::R_user_dir("enrichmet", "cache")  

# Initialize BiocFileCache without prompting
bfc <- BiocFileCache(cache_dir, ask = FALSE)

# Zenodo-hosted STITCH interactions file
url <- "https://zenodo.org/records/15498097/files/chemical_chemical.tsv"

# Retrieve or download + cache the file
path <- bfcrpath(bfc, url)

# Read the file
stitch_df <- read_tsv(path)

head(stitch_df)


test_that("enrichmet returns expected output format", {
  result <- enrichmet(inputMetabolites, PathwayVsMetabolites=PathwayVsMetabolites, example_data, top_n = 100, p_value_cutoff = 1, kegg_lookup = kegg_lookup, mapping_df = mapping_df, stitch_df = stitch_df)
  expect_type(result, "list")
  expect_named(result, c("pathway_plot", "impact_plot", "MetSEA_plot", "rbc_plot", "network_plot", "heatmap_plot", "membership_plot", "interaction_plot"))
})
