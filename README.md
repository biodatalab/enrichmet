
<!-- README.md is generated from README.Rmd. Please edit that file -->

# enrichmet

<!-- badges: start -->
<!-- badges: end -->

**enrichmet** simplifies pathway enrichment analysis by allowing the
complete workflow to be executed through a single R function call. This
design eliminates repetitive steps such as data reformatting and
parameter configuration, improving efficiency, reducing the risk of
errors, and supporting reproducible analysis. For users who wish to run
only a specific analysis or generate selected plots rather than
executing the full workflow, enrichmet also provides dedicated modules
for each individual analysis and visualization.

**enrichmet performs pathway enrichment analysis using Fisher’s exact
test, computes betweenness centrality for metabolites, and performs
Metabolite Set Enrichment Analysis (MetSEA). The **enrichmet\*\*()
function produces three tables (S3 data.frame objects), which may
include the MetSEA table, metabolite centrality, and pathway enrichment
results. In addition, it generates eight plots (S3/S4 plot objects):

- **Pathway enrichment plot**
- **Pathway impact plot**
- **Metabolite Set Enrichment Analysis (MetSEA plot)**
- **Relative Betweenness Centrality (RBC) plot**
- **Network graph**
- **Pathway heatmap**
- **Pathway membership plot**
- **Interaction network plot**

## Installation

You can install enrichmet as:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("enrichmet")
```

## Example

This is a basic example

``` r
## ** Example showing ALL 8 plots
# Load required libraries
library(enrichmet)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(stringr)
library(tidyr)

## Generate example data with proper KEGG IDs
set.seed(1234)

# Create 20 unique metabolites WITH KEGG IDs
inputMetabolites <- paste0("C", sprintf("%05d", 1:20))

# ---- 1. Generate 50 pathways with random metabolites assigned ----
pathway_names <- paste0("Pathway", 1:50)
PathwayVsMetabolites <- data.frame(
    Pathway = rep(pathway_names, each = 1),
    Metabolites = sapply(1:50, function(x) paste(sample(inputMetabolites, sample(5:15, 1)), collapse = ","))
)

# ---- 2. Add new pathway entries ----
new_rows <- data.frame(
    Pathway = c("Pathway101", "Pathway102", "Pathway103", "Pathway104", "Pathway105"),
    Metabolites = c(
        "C00012,C00013,C00014,C00015,C00016,C00001,C00018,C00003,C00006,C00004",
        "C00006,C00007,C00008,C00009,C00010,C00011,C00009,C00006,C00016,C00004",
        "C00024,C00025,C00026,C00027,C00028,C00029,C00030,C00026,C00005",
        "C00013,C00014,C00015,C00016,C00017,C00024,C00027,C00014",
        "C00015,C00016,C00017,C00018,C00019,C00020,C00021,C00004,C00008,C00010"
    )
)

# Combine with existing PathwayVsMetabolites
PathwayVsMetabolites <- rbind(PathwayVsMetabolites, new_rows)

# ---- 3. Generate example metabolite-level data (with KEGG IDs) ----
example_data <- data.frame(
    met_id = c(inputMetabolites, "C00021", "C00022", "C00023"),  # Add extra metabolites
    pval = c(runif(20, 0.001, 0.05), 0.0001, 0.0002, 0.0003),
    log2fc = c(rnorm(20, mean = 0, sd = 1), 2.5, 2.3, -2.1),
    mz = rnorm(23, mean = 200, sd = 50),
    rt = rnorm(23, mean = 10, sd = 3)
)

# ---- 4. Create KEGG lookup table (ESSENTIAL) ----
kegg_lookup <- data.frame(
    kegg_id = c(inputMetabolites, "C00021", "C00022", "C00023"),
    name = c(
        "Glucose", "Lactate", "Pyruvate", "Alanine", "Valine",
        "Leucine", "Isoleucine", "Proline", "Serine", "Threonine",
        "Cysteine", "Methionine", "Aspartate", "Glutamate", "Asparagine",
        "Glutamine", "Lysine", "Arginine", "Histidine", "Phenylalanine",
        "Tyrosine", "Tryptophan", "Creatine"
    )
)

# ---- 5. Create mapping_df with additional IDs ----
set.seed(42)
mapping_df <- data.frame(
    KEGG_ID = inputMetabolites,
    HMDB_ID = paste0("HMDB", stringr::str_pad(sample(10000:99999, length(inputMetabolites)), 6, pad = "0")),
    PubChem_CID = as.character(sample(10000:99999, length(inputMetabolites))),
    CHEBI_ID = paste0("CHEBI:", sample(10000:99999, length(inputMetabolites))),
    STITCH_ID = paste0("CIDs", stringr::str_pad(sample(1000:9999, length(inputMetabolites)), 8, pad = "0")),
    Compound_Name = paste("Compound", sample(LETTERS, length(inputMetabolites), replace = TRUE), 
                          sample(100:999, length(inputMetabolites), replace = TRUE))
)

# ---- 6. Create synthetic STITCH interaction data ----
stitch_ids <- mapping_df$STITCH_ID

stitch_pairs <- expand.grid(chemical1 = stitch_ids, chemical2 = stitch_ids) %>%
    dplyr::filter(chemical1 != chemical2)

set.seed(123)
stitch_df <- stitch_pairs %>%
    dplyr::slice_sample(n = 50) %>%  # Reduced to 50 for faster processing
    dplyr::mutate(
        similarity = runif(dplyr::n(), 0.6, 0.95),  # Higher similarity for better networks
        experimental = sample(200:800, dplyr::n(), replace = TRUE),
        database = sample(c(0, 300, 600, 900), dplyr::n(), replace = TRUE),
        textmining = sample(0:900, dplyr::n(), replace = TRUE),
        combined_score = similarity * 200 + experimental + database + textmining
    )

# ---- 7. Run comprehensive enrichment analysis with ALL features ----
results <- enrichmet(
    inputMetabolites = inputMetabolites,
    PathwayVsMetabolites = PathwayVsMetabolites,
    example_data = example_data,
    kegg_lookup = kegg_lookup,
    top_n = 15,
    p_value_cutoff = 0.1,  # Use a reasonable cutoff
    mapping_df = mapping_df,
    stitch_df = stitch_df,
    analysis_type = c("enrichment", "gsea", "centrality", "network",
                      "heatmap", "membership", "interaction"),
    network_top_n = 10,
    heatmap_top_n = 10,
    membership_top_n = 10,
    min_pathway_occurrence = 2,
    min_metabolite_occurrence = 2
)
#> Using 20 metabolites from inputMetabolites character vector
#> Processing complex KEGG IDs (splitting by |)...
#> Split 20 complex IDs into 20 unique KEGG IDs
#> Running pathway enrichment analysis...
#> Enrichment analysis completed: 55 pathways tested, 55 pathways passed filtering (p <= 1.000000)
#> Running GSEA analysis...
#> Extracting KEGG IDs from 'met_id' column...
#> Successfully extracted KEGG IDs for 23 metabolites
#> Sample extracted KEGG IDs: C00001, C00002, C00003, C00004, C00005, C00006
#> Prepared 23 metabolites with KEGG IDs for GSEA
#> Created rankings for 23 KEGG metabolites using 'signed_pval' method
#> Ranking range: -3.523 to 4
#> Preparing GMT data from 55 pathways
#> GMT preparation summary:
#>   Total pathways processed: 55
#>   Pathways kept (size >= 5): 54
#>   Pathways removed (size < 5): 1
#>   Pathway size distribution:
#>     Min: 5
#>     Median: 10
#>     Max: 15
#>     Mean: 10.2
#> Testing 54 pathways with GSEA
#> GSEA completed: 54 pathways tested, 0 significant at FDR < 0.05
#> Applied KEGG pathway name mapping using 'kegg_id' and 'name' columns
#> Running centrality analysis...
#> Generating metabolite-pathway network visualization...
#> Using top 10 metabolites by centrality for network plot
#> Generating enrichment heatmap...
#> Selected top 10 metabolites by enrichment significance for heatmap
#> Generating pathway membership plot...
#> === MEMBERSHIP PLOT DEBUG ===
#> Input metabolites: 20
#> Top_n parameter: 10
#> Initial matched metabolites: 20
#> Filtering to top 10 metabolites by centrality...
#> Total metabolites with centrality: 28
#> Top metabolites selected: 10
#> Sample top metabolites: C00005, C00015, C00013, C00004, C00017, C00014
#> After centrality filtering - unique metabolites: 10
#> Matrix before occurrence filtering: 10 x 55
#> Pathway counts range: 1 to 10
#> Metabolite counts range: 27 to 31
#> Pathways to keep: 54/55
#> Metabolites to keep: 10/10
#> Final matrix dimensions: 10 x 54
#> === END DEBUG ===
#> Generating STITCH interaction network...
#> Extracting KEGG IDs from input metabolites...
#> Successfully extracted 20 KEGG IDs from 20 input metabolites
#> Sample extracted KEGG IDs: C00001, C00002, C00003, C00004, C00005, C00006
#> Found 20 metabolites in mapping_df with valid PubChem CIDs
#> Applied KEGG pathway name mapping
#> Display names sample: Glucose, Lactate, Pyruvate, Alanine, Valine, Leucine
#> Found 50 valid interactions between 20 metabolites
#> Creating graph with 20 vertices
#> Vertex attributes: name, display_name, KEGG_ID, PubChem_CID
#> Graph vertex attributes: name, display_name, KEGG_ID, PubChem_CID
#> Sample vertex display_names: Glucose, Lactate, Pyruvate, Alanine, Valine, Leucine
#> Using layout: gem (spacing score: 388.19)
results
#> $input_metabolites_used
#>  [1] "C00001" "C00002" "C00003" "C00004" "C00005" "C00006" "C00007" "C00008"
#>  [9] "C00009" "C00010" "C00011" "C00012" "C00013" "C00014" "C00015" "C00016"
#> [17] "C00017" "C00018" "C00019" "C00020"
#> 
#> $pathway_enrichment_all
#>       Pathway      P_value  Log_P_value    Impact  Coverage Count Pathway_Size
#> 1    Pathway1 0.0009661836 3.014940e+00 1.0000000 1.0000000    14           14
#> 2    Pathway2 0.0405295188 1.392229e+00 1.0000000 1.0000000     8            8
#> 3    Pathway3 0.0009661836 3.014940e+00 1.0000000 1.0000000    14           14
#> 4    Pathway4 0.0004140787 3.382917e+00 1.0000000 1.0000000    15           15
#> 5    Pathway5 0.0140786749 1.851438e+00 1.0000000 1.0000000    10           10
#> 6    Pathway6 0.0020703934 2.683947e+00 1.0000000 1.0000000    13           13
#> 7    Pathway7 0.1577533578 8.020214e-01 1.0000000 1.0000000     5            5
#> 8    Pathway8 0.0140786749 1.851438e+00 1.0000000 1.0000000    10           10
#> 9    Pathway9 0.0405295188 1.392229e+00 1.0000000 1.0000000     8            8
#> 10  Pathway10 0.0020703934 2.683947e+00 1.0000000 1.0000000    13           13
#> 11  Pathway11 0.0243177113 1.614077e+00 1.0000000 1.0000000     9            9
#> 12  Pathway12 0.0009661836 3.014940e+00 1.0000000 1.0000000    14           14
#> 13  Pathway13 0.0654707611 1.183953e+00 1.0000000 1.0000000     7            7
#> 14  Pathway14 0.0009661836 3.014940e+00 1.0000000 1.0000000    14           14
#> 15  Pathway15 0.0654707611 1.183953e+00 1.0000000 1.0000000     7            7
#> 16  Pathway16 0.0004140787 3.382917e+00 1.0000000 1.0000000    15           15
#> 17  Pathway17 0.0243177113 1.614077e+00 1.0000000 1.0000000     9            9
#> 18  Pathway18 0.0405295188 1.392229e+00 1.0000000 1.0000000     8            8
#> 19  Pathway19 0.0654707611 1.183953e+00 1.0000000 1.0000000     7            7
#> 20  Pathway20 0.0140786749 1.851438e+00 1.0000000 1.0000000    10           10
#> 21  Pathway21 0.0004140787 3.382917e+00 1.0000000 1.0000000    15           15
#> 22  Pathway22 0.0243177113 1.614077e+00 1.0000000 1.0000000     9            9
#> 23  Pathway23 0.0243177113 1.614077e+00 1.0000000 1.0000000     9            9
#> 24  Pathway24 0.0243177113 1.614077e+00 1.0000000 1.0000000     9            9
#> 25  Pathway25 0.1028826246 9.876580e-01 1.0000000 1.0000000     6            6
#> 26  Pathway26 0.0405295188 1.392229e+00 1.0000000 1.0000000     8            8
#> 27  Pathway27 0.0004140787 3.382917e+00 1.0000000 1.0000000    15           15
#> 28  Pathway28 0.0405295188 1.392229e+00 1.0000000 1.0000000     8            8
#> 29  Pathway29 0.0009661836 3.014940e+00 1.0000000 1.0000000    14           14
#> 30  Pathway30 0.1028826246 9.876580e-01 1.0000000 1.0000000     6            6
#> 31  Pathway31 0.0020703934 2.683947e+00 1.0000000 1.0000000    13           13
#> 32  Pathway32 0.0004140787 3.382917e+00 1.0000000 1.0000000    15           15
#> 33  Pathway33 0.0009661836 3.014940e+00 1.0000000 1.0000000    14           14
#> 34  Pathway34 0.0654707611 1.183953e+00 1.0000000 1.0000000     7            7
#> 35  Pathway35 0.0405295188 1.392229e+00 1.0000000 1.0000000     8            8
#> 36  Pathway36 0.0243177113 1.614077e+00 1.0000000 1.0000000     9            9
#> 37  Pathway37 0.0020703934 2.683947e+00 1.0000000 1.0000000    13           13
#> 38  Pathway38 0.0078214861 2.106711e+00 1.0000000 1.0000000    11           11
#> 39  Pathway39 0.0405295188 1.392229e+00 1.0000000 1.0000000     8            8
#> 40  Pathway40 0.0140786749 1.851438e+00 1.0000000 1.0000000    10           10
#> 41  Pathway41 0.0140786749 1.851438e+00 1.0000000 1.0000000    10           10
#> 42  Pathway42 0.0078214861 2.106711e+00 1.0000000 1.0000000    11           11
#> 43  Pathway43 0.1577533578 8.020214e-01 1.0000000 1.0000000     5            5
#> 44  Pathway44 0.0041407867 2.382917e+00 1.0000000 1.0000000    12           12
#> 45  Pathway45 0.0004140787 3.382917e+00 1.0000000 1.0000000    15           15
#> 46  Pathway46 0.0140786749 1.851438e+00 1.0000000 1.0000000    10           10
#> 47  Pathway47 0.1028826246 9.876580e-01 1.0000000 1.0000000     6            6
#> 48  Pathway48 0.0405295188 1.392229e+00 1.0000000 1.0000000     8            8
#> 49  Pathway49 0.0009661836 3.014940e+00 1.0000000 1.0000000    14           14
#> 50  Pathway50 0.0140786749 1.851438e+00 1.0000000 1.0000000    10           10
#> 51 Pathway101 0.0140786749 1.851438e+00 1.0000000 1.0000000    10           10
#> 52 Pathway102 0.0405295188 1.392229e+00 1.0000000 1.0000000     8            8
#> 53 Pathway103 0.9999996783 1.397297e-07 0.9676646 0.1250000     1            8
#> 54 Pathway104 0.6939900679 1.586467e-01 0.9775252 0.7142857     5            7
#> 55 Pathway105 0.1164690382 9.337895e-01 1.0000000 0.9000000     9           10
#>    Input_Size Adjusted_P_value      Q_value Enrichment_Ratio
#> 1          20      0.004087700 0.0008861659            1.400
#> 2          20      0.051840082 0.0112383292            1.400
#> 3          20      0.004087700 0.0008861659            1.400
#> 4          20      0.003795721 0.0008228684            1.400
#> 5          20      0.027654540 0.0059951839            1.400
#> 6          20      0.006698332 0.0014521207            1.400
#> 7          20      0.163706315 0.0354896324            1.400
#> 8          20      0.027654540 0.0059951839            1.400
#> 9          20      0.051840082 0.0112383292            1.400
#> 10         20      0.006698332 0.0014521207            1.400
#> 11         20      0.039337474 0.0085279087            1.400
#> 12         20      0.004087700 0.0008861659            1.400
#> 13         20      0.076614720 0.0166091838            1.400
#> 14         20      0.004087700 0.0008861659            1.400
#> 15         20      0.076614720 0.0166091838            1.400
#> 16         20      0.003795721 0.0008228684            1.400
#> 17         20      0.039337474 0.0085279087            1.400
#> 18         20      0.051840082 0.0112383292            1.400
#> 19         20      0.076614720 0.0166091838            1.400
#> 20         20      0.027654540 0.0059951839            1.400
#> 21         20      0.003795721 0.0008228684            1.400
#> 22         20      0.039337474 0.0085279087            1.400
#> 23         20      0.039337474 0.0085279087            1.400
#> 24         20      0.039337474 0.0085279087            1.400
#> 25         20      0.113170887 0.0245341372            1.400
#> 26         20      0.051840082 0.0112383292            1.400
#> 27         20      0.003795721 0.0008228684            1.400
#> 28         20      0.051840082 0.0112383292            1.400
#> 29         20      0.004087700 0.0008861659            1.400
#> 30         20      0.113170887 0.0245341372            1.400
#> 31         20      0.006698332 0.0014521207            1.400
#> 32         20      0.003795721 0.0008228684            1.400
#> 33         20      0.004087700 0.0008861659            1.400
#> 34         20      0.076614720 0.0166091838            1.400
#> 35         20      0.051840082 0.0112383292            1.400
#> 36         20      0.039337474 0.0085279087            1.400
#> 37         20      0.006698332 0.0014521207            1.400
#> 38         20      0.021509087 0.0046629208            1.400
#> 39         20      0.051840082 0.0112383292            1.400
#> 40         20      0.027654540 0.0059951839            1.400
#> 41         20      0.027654540 0.0059951839            1.400
#> 42         20      0.021509087 0.0046629208            1.400
#> 43         20      0.163706315 0.0354896324            1.400
#> 44         20      0.012652404 0.0027428946            1.400
#> 45         20      0.003795721 0.0008228684            1.400
#> 46         20      0.027654540 0.0059951839            1.400
#> 47         20      0.113170887 0.0245341372            1.400
#> 48         20      0.051840082 0.0112383292            1.400
#> 49         20      0.004087700 0.0008861659            1.400
#> 50         20      0.027654540 0.0059951839            1.400
#> 51         20      0.027654540 0.0059951839            1.400
#> 52         20      0.051840082 0.0112383292            1.400
#> 53         20      0.999999678 0.2167883450            0.175
#> 54         20      0.706841736 0.1532350994            1.000
#> 55         20      0.125603865 0.0272294627            1.260
#>                                                                                             Metabolite_List
#> 1         C00005,C00012,C00015,C00009,C00020,C00006,C00004,C00002,C00007,C00018,C00010,C00011,C00019,C00017
#> 2                                                   C00014,C00004,C00019,C00008,C00018,C00017,C00003,C00016
#> 3         C00005,C00002,C00015,C00008,C00011,C00004,C00012,C00003,C00007,C00009,C00013,C00006,C00018,C00019
#> 4  C00002,C00015,C00017,C00006,C00008,C00003,C00018,C00001,C00013,C00009,C00016,C00012,C00007,C00019,C00020
#> 5                                     C00003,C00009,C00016,C00020,C00006,C00019,C00010,C00007,C00018,C00008
#> 6                C00003,C00019,C00018,C00006,C00012,C00004,C00009,C00007,C00020,C00017,C00005,C00008,C00014
#> 7                                                                        C00004,C00019,C00009,C00017,C00006
#> 8                                     C00006,C00013,C00017,C00002,C00014,C00020,C00010,C00018,C00005,C00011
#> 9                                                   C00003,C00004,C00010,C00006,C00017,C00009,C00011,C00020
#> 10               C00014,C00008,C00017,C00013,C00009,C00002,C00006,C00019,C00011,C00020,C00003,C00007,C00016
#> 11                                           C00014,C00019,C00007,C00009,C00011,C00012,C00016,C00008,C00003
#> 12        C00016,C00002,C00005,C00020,C00018,C00009,C00003,C00008,C00014,C00013,C00012,C00015,C00004,C00011
#> 13                                                         C00007,C00017,C00003,C00002,C00005,C00015,C00014
#> 14        C00017,C00003,C00016,C00009,C00004,C00012,C00005,C00001,C00019,C00006,C00014,C00011,C00002,C00020
#> 15                                                         C00018,C00011,C00016,C00006,C00007,C00008,C00004
#> 16 C00015,C00017,C00011,C00007,C00010,C00008,C00018,C00003,C00016,C00002,C00020,C00001,C00012,C00006,C00014
#> 17                                           C00020,C00012,C00010,C00001,C00013,C00006,C00015,C00011,C00003
#> 18                                                  C00004,C00011,C00003,C00001,C00018,C00013,C00017,C00009
#> 19                                                         C00010,C00013,C00012,C00001,C00005,C00019,C00020
#> 20                                    C00013,C00008,C00018,C00002,C00019,C00006,C00009,C00005,C00015,C00020
#> 21 C00014,C00009,C00001,C00018,C00016,C00015,C00006,C00005,C00010,C00002,C00019,C00011,C00004,C00003,C00007
#> 22                                           C00002,C00006,C00020,C00005,C00014,C00007,C00013,C00010,C00004
#> 23                                           C00009,C00006,C00015,C00020,C00002,C00011,C00008,C00012,C00014
#> 24                                           C00019,C00007,C00006,C00005,C00015,C00016,C00004,C00020,C00012
#> 25                                                                C00010,C00004,C00018,C00020,C00002,C00001
#> 26                                                  C00010,C00019,C00013,C00014,C00009,C00018,C00011,C00007
#> 27 C00019,C00010,C00013,C00001,C00006,C00005,C00012,C00016,C00003,C00014,C00017,C00020,C00002,C00007,C00015
#> 28                                                  C00014,C00020,C00017,C00010,C00016,C00007,C00011,C00006
#> 29        C00017,C00014,C00006,C00015,C00007,C00005,C00020,C00012,C00011,C00019,C00018,C00010,C00004,C00003
#> 30                                                                C00016,C00012,C00015,C00008,C00004,C00013
#> 31               C00019,C00001,C00004,C00014,C00015,C00010,C00008,C00013,C00011,C00016,C00003,C00005,C00007
#> 32 C00004,C00017,C00014,C00018,C00015,C00002,C00006,C00008,C00020,C00011,C00003,C00019,C00007,C00010,C00016
#> 33        C00018,C00003,C00004,C00015,C00005,C00006,C00002,C00008,C00014,C00016,C00017,C00013,C00020,C00001
#> 34                                                         C00008,C00012,C00003,C00015,C00001,C00017,C00018
#> 35                                                  C00011,C00007,C00004,C00013,C00016,C00015,C00001,C00020
#> 36                                           C00004,C00015,C00018,C00007,C00002,C00005,C00006,C00003,C00010
#> 37               C00007,C00016,C00012,C00011,C00010,C00009,C00013,C00015,C00020,C00008,C00019,C00018,C00017
#> 38                             C00005,C00006,C00001,C00019,C00013,C00010,C00014,C00008,C00004,C00017,C00002
#> 39                                                  C00010,C00015,C00001,C00018,C00017,C00004,C00003,C00016
#> 40                                    C00013,C00011,C00017,C00015,C00012,C00006,C00020,C00009,C00001,C00003
#> 41                                    C00020,C00018,C00017,C00002,C00001,C00005,C00011,C00010,C00013,C00007
#> 42                             C00005,C00011,C00014,C00007,C00010,C00018,C00015,C00017,C00013,C00012,C00019
#> 43                                                                       C00003,C00015,C00010,C00001,C00013
#> 44                      C00017,C00008,C00015,C00016,C00010,C00002,C00001,C00012,C00003,C00018,C00007,C00014
#> 45 C00008,C00019,C00015,C00007,C00013,C00014,C00020,C00005,C00003,C00006,C00011,C00001,C00009,C00017,C00012
#> 46                                    C00003,C00020,C00007,C00011,C00012,C00005,C00016,C00009,C00004,C00013
#> 47                                                                C00014,C00008,C00005,C00004,C00013,C00007
#> 48                                                  C00013,C00014,C00009,C00005,C00011,C00006,C00018,C00002
#> 49        C00005,C00001,C00007,C00009,C00004,C00018,C00010,C00015,C00012,C00019,C00017,C00020,C00002,C00008
#> 50                                    C00018,C00004,C00010,C00014,C00012,C00013,C00005,C00006,C00016,C00019
#> 51                                    C00012,C00013,C00014,C00015,C00016,C00001,C00018,C00003,C00006,C00004
#> 52                                                  C00006,C00007,C00008,C00009,C00010,C00011,C00016,C00004
#> 53                                                                                                   C00005
#> 54                                                                       C00013,C00014,C00015,C00016,C00017
#> 55                                           C00015,C00016,C00017,C00018,C00019,C00020,C00004,C00008,C00010
#> 
#> $pathway_enrichment_results
#>      Pathway      P_value Log_P_value Impact Coverage Count Pathway_Size
#> 1   Pathway4 0.0004140787    3.382917      1        1    15           15
#> 2  Pathway16 0.0004140787    3.382917      1        1    15           15
#> 3  Pathway21 0.0004140787    3.382917      1        1    15           15
#> 4  Pathway27 0.0004140787    3.382917      1        1    15           15
#> 5  Pathway32 0.0004140787    3.382917      1        1    15           15
#> 6  Pathway45 0.0004140787    3.382917      1        1    15           15
#> 7   Pathway1 0.0009661836    3.014940      1        1    14           14
#> 8   Pathway3 0.0009661836    3.014940      1        1    14           14
#> 9  Pathway12 0.0009661836    3.014940      1        1    14           14
#> 10 Pathway14 0.0009661836    3.014940      1        1    14           14
#> 11 Pathway29 0.0009661836    3.014940      1        1    14           14
#> 12 Pathway33 0.0009661836    3.014940      1        1    14           14
#> 13 Pathway49 0.0009661836    3.014940      1        1    14           14
#> 14  Pathway6 0.0020703934    2.683947      1        1    13           13
#> 15 Pathway10 0.0020703934    2.683947      1        1    13           13
#>    Input_Size Adjusted_P_value      Q_value Enrichment_Ratio
#> 1          20      0.003795721 0.0008228684              1.4
#> 2          20      0.003795721 0.0008228684              1.4
#> 3          20      0.003795721 0.0008228684              1.4
#> 4          20      0.003795721 0.0008228684              1.4
#> 5          20      0.003795721 0.0008228684              1.4
#> 6          20      0.003795721 0.0008228684              1.4
#> 7          20      0.004087700 0.0008861659              1.4
#> 8          20      0.004087700 0.0008861659              1.4
#> 9          20      0.004087700 0.0008861659              1.4
#> 10         20      0.004087700 0.0008861659              1.4
#> 11         20      0.004087700 0.0008861659              1.4
#> 12         20      0.004087700 0.0008861659              1.4
#> 13         20      0.004087700 0.0008861659              1.4
#> 14         20      0.006698332 0.0014521207              1.4
#> 15         20      0.006698332 0.0014521207              1.4
#>                                                                                             Metabolite_List
#> 1  C00002,C00015,C00017,C00006,C00008,C00003,C00018,C00001,C00013,C00009,C00016,C00012,C00007,C00019,C00020
#> 2  C00015,C00017,C00011,C00007,C00010,C00008,C00018,C00003,C00016,C00002,C00020,C00001,C00012,C00006,C00014
#> 3  C00014,C00009,C00001,C00018,C00016,C00015,C00006,C00005,C00010,C00002,C00019,C00011,C00004,C00003,C00007
#> 4  C00019,C00010,C00013,C00001,C00006,C00005,C00012,C00016,C00003,C00014,C00017,C00020,C00002,C00007,C00015
#> 5  C00004,C00017,C00014,C00018,C00015,C00002,C00006,C00008,C00020,C00011,C00003,C00019,C00007,C00010,C00016
#> 6  C00008,C00019,C00015,C00007,C00013,C00014,C00020,C00005,C00003,C00006,C00011,C00001,C00009,C00017,C00012
#> 7         C00005,C00012,C00015,C00009,C00020,C00006,C00004,C00002,C00007,C00018,C00010,C00011,C00019,C00017
#> 8         C00005,C00002,C00015,C00008,C00011,C00004,C00012,C00003,C00007,C00009,C00013,C00006,C00018,C00019
#> 9         C00016,C00002,C00005,C00020,C00018,C00009,C00003,C00008,C00014,C00013,C00012,C00015,C00004,C00011
#> 10        C00017,C00003,C00016,C00009,C00004,C00012,C00005,C00001,C00019,C00006,C00014,C00011,C00002,C00020
#> 11        C00017,C00014,C00006,C00015,C00007,C00005,C00020,C00012,C00011,C00019,C00018,C00010,C00004,C00003
#> 12        C00018,C00003,C00004,C00015,C00005,C00006,C00002,C00008,C00014,C00016,C00017,C00013,C00020,C00001
#> 13        C00005,C00001,C00007,C00009,C00004,C00018,C00010,C00015,C00012,C00019,C00017,C00020,C00002,C00008
#> 14               C00003,C00019,C00018,C00006,C00012,C00004,C00009,C00007,C00020,C00017,C00005,C00008,C00014
#> 15               C00014,C00008,C00017,C00013,C00009,C00002,C00006,C00019,C00011,C00020,C00003,C00007,C00016
#> 
#> $pathway_plot
```

<img src="man/figures/README-example-1.png" width="100%" />

    #> 
    #> $impact_plot

<img src="man/figures/README-example-2.png" width="100%" />

    #> 
    #> $gsea_results
    #>        pathway       pval      padj    log2err         ES        NES  size
    #>         <char>      <num>     <num>      <num>      <num>      <num> <int>
    #>  1:  Pathway14 0.03438555 0.7948800 0.32177592 -0.5555556 -1.6626293    14
    #>  2:  Pathway36 0.05400982 0.7948800 0.24891111  0.5757292  1.5332087     9
    #>  3:  Pathway35 0.05505109 0.7948800 0.32177592 -0.5353320 -1.5107447     8
    #>  4: Pathway104 0.06997743 0.7948800 0.25720647 -0.6111111 -1.5165348     5
    #>  5:  Pathway44 0.07360000 0.7948800 0.20895503  0.5478224  1.4848700    12
    #>  6:  Pathway40 0.11413043 0.8915094 0.21925035 -0.4615385 -1.3763556    10
    #>  7:  Pathway18 0.13131313 0.8915094 0.19578900 -0.4698021 -1.3258147     8
    #>  8:  Pathway10 0.13207547 0.8915094 0.20207171 -0.4262116 -1.2924156    13
    #>  9:   Pathway4 0.16825397 0.9145161 0.13284630  0.5000000  1.3163859    15
    #> 10:  Pathway45 0.16935484 0.9145161 0.17669427 -0.4207094 -1.2303829    15
    #> 11:  Pathway34 0.20000000 0.9782609 0.12443417  0.5000000  1.2562309     7
    #> 12: Pathway105 0.24921136 0.9782609 0.10552094  0.4530527  1.2168795    10
    #> 13:  Pathway46 0.26630435 0.9782609 0.13880511 -0.3846154 -1.1469630    10
    #> 14:  Pathway16 0.28888889 0.9782609 0.09688777  0.4293675  1.1304266    15
    #> 15:   Pathway7 0.29119639 0.9782609 0.11881504 -0.4595606 -1.1404467     5
    #> 16:  Pathway27 0.29301075 0.9782609 0.13077714 -0.3884136 -1.1359325    15
    #> 17:   Pathway5 0.40063091 0.9782609 0.07829552  0.3846154  1.0330599    10
    #> 18:  Pathway26 0.40404040 0.9782609 0.10473282 -0.3597180 -1.0151497     8
    #> 19:  Pathway24 0.40759494 0.9782609 0.10434395 -0.3571429 -1.0584396     9
    #> 20:  Pathway47 0.41362530 0.9782609 0.10099059 -0.4023854 -1.0347550     6
    #> 21:  Pathway12 0.41666667 0.9782609 0.10672988 -0.3333333 -0.9975776    14
    #> 22:   Pathway3 0.41666667 0.9782609 0.10672988 -0.3333333 -0.9975776    14
    #> 23:  Pathway33 0.41666667 0.9782609 0.10672988 -0.3333333 -0.9975776    14
    #> 24:   Pathway6 0.56334232 0.9974747 0.08889453 -0.3000000 -0.9097001    13
    #> 25: Pathway101 0.56521739 0.9974747 0.08916471 -0.3076923 -0.9175704    10
    #> 26:  Pathway49 0.57301587 0.9974747 0.06077195  0.3480110  0.9327974    14
    #> 27:  Pathway11 0.61012658 0.9974747 0.08108021 -0.2992943 -0.8869979     9
    #> 28:  Pathway13 0.62521008 0.9974747 0.05934877  0.3385562  0.8506094     7
    #> 29:  Pathway31 0.63072776 0.9974747 0.08266464 -0.2913365 -0.8834296    13
    #> 30:  Pathway37 0.68720379 0.9974747 0.05205700  0.3000000  0.8215613    13
    #> 31:  Pathway20 0.69400631 0.9974747 0.05153091  0.3076923  0.8264479    10
    #> 32:  Pathway21 0.71505376 0.9974747 0.07588869 -0.2683539 -0.7848126    15
    #> 33:  Pathway50 0.71766562 0.9974747 0.04999139  0.3001043  0.8060668    10
    #> 34:  Pathway39 0.73146623 0.9974747 0.05111480  0.3019127  0.7805432     8
    #> 35:  Pathway48 0.80808081 0.9974747 0.06658921 -0.2666667 -0.7525522     8
    #> 36:  Pathway32 0.84408602 0.9974747 0.06751890 -0.2500000 -0.7311359    15
    #> 37:  Pathway23 0.86075949 0.9974747 0.06364241 -0.2277685 -0.6750217     9
    #> 38:   Pathway1 0.88440860 0.9974747 0.06523531 -0.2222222 -0.6650517    14
    #> 39:  Pathway29 0.88440860 0.9974747 0.06523531 -0.2222222 -0.6650517    14
    #> 40:  Pathway22 0.88860759 0.9974747 0.06211242 -0.2172001 -0.6437008     9
    #> 41:  Pathway15 0.89915966 0.9974747 0.04258778  0.2585972  0.6497156     7
    #> 42:   Pathway9 0.90151515 0.9974747 0.06130261 -0.2394046 -0.6756165     8
    #> 43:  Pathway41 0.90378549 0.9974747 0.03943665  0.2307692  0.6198359    10
    #> 44:   Pathway2 0.90909091 0.9974747 0.06090393 -0.2349751 -0.6631162     8
    #> 45:  Pathway19 0.92268908 0.9974747 0.04140443  0.2500000  0.6281155     7
    #> 46:   Pathway8 0.92271293 0.9974747 0.03847869  0.2261323  0.6073812    10
    #> 47:  Pathway38 0.93103448 0.9974747 0.06211242 -0.2012729 -0.6215045    11
    #> 48:  Pathway17 0.95090016 0.9974747 0.03879622  0.2142857  0.5706585     9
    #> 49:  Pathway42 0.95755968 0.9974747 0.06077195 -0.1800970 -0.5561160    11
    #> 50:  Pathway25 0.96134454 0.9974747 0.03951722  0.2278119  0.5478896     6
    #> 51: Pathway102 0.97364086 0.9974747 0.03800562  0.2099911  0.5428958     8
    #> 52:  Pathway43 0.98223801 0.9974747 0.04107133  0.2222222  0.4979655     5
    #> 53:  Pathway30 0.98296837 0.9974747 0.05547933 -0.2065311 -0.5311055     6
    #> 54:  Pathway28 0.99747475 0.9974747 0.05652995 -0.1537743 -0.4339619     8
    #>        pathway       pval      padj    log2err         ES        NES  size
    #>      leadingEdge input_count significance
    #>           <list>       <int>       <char>
    #>  1: C00004, ....          14           NS
    #>  2: C00010, ....           8           NS
    #>  3: C00004, ....           7           NS
    #>  4: C00014, ....           5           NS
    #>  5: C00010, ....          11           NS
    #>  6: C00011, ....          10           NS
    #>  7: C00004, ....           6           NS
    #>  8: C00014, ....          11           NS
    #>  9: C00008, ....          15           NS
    #> 10: C00014, ....          14           NS
    #> 11: C00008, ....           7           NS
    #> 12: C00021, ....           4           NS
    #> 13: C00004, ....          10           NS
    #> 14: C00010, ....          12           NS
    #> 15: C00004, ....           4           NS
    #> 16: C00014, ....          14           NS
    #> 17: C00010, ....          10           NS
    #> 18: C00014, ....           5           NS
    #> 19: C00004, ....           9           NS
    #> 20: C00004, ....           3           NS
    #> 21: C00004, ....          14           NS
    #> 22: C00004, ....          14           NS
    #> 23: C00004, ....          14           NS
    #> 24: C00004, ....          13           NS
    #> 25: C00004, ....          10           NS
    #> 26: C00010, ....          13           NS
    #> 27: C00014, ....           8           NS
    #> 28: C00005, ....           6           NS
    #> 29: C00004, ....           5           NS
    #> 30: C00010, ....          13           NS
    #> 31: C00008, ....          10           NS
    #> 32: C00004, ....          14           NS
    #> 33: C00010, ....           5           NS
    #> 34: C00010, ....           7           NS
    #> 35: C00014, ....           8           NS
    #> 36: C00004, ....          15           NS
    #> 37: C00014, ....           8           NS
    #> 38: C00004, ....          14           NS
    #> 39: C00004, ....          14           NS
    #> 40: C00004, ....           4           NS
    #> 41: C00008, ....           4           NS
    #> 42: C00004, ....           5           NS
    #> 43: C00010, ....          10           NS
    #> 44: C00004, ....           2           NS
    #> 45: C00010, ....           7           NS
    #> 46: C00010, ....           4           NS
    #> 47: C00004, ....           4           NS
    #> 48: C00010, ....           9           NS
    #> 49: C00014, ....           4           NS
    #> 50: C00010, ....           2           NS
    #> 51: C00010, ....           4           NS
    #> 52: C00010, ....           5           NS
    #> 53: C00004, ....           2           NS
    #> 54: C00014, ....           5           NS
    #>      leadingEdge input_count significance
    #> 
    #> $gsea_plot

<img src="man/figures/README-example-3.png" width="100%" />

    #> 
    #> $metabolite_centrality
    #>    Metabolite RBC_Metabolite  Display_Name
    #> 1      C00005     0.15717641        Valine
    #> 2      C00015     0.05152146    Asparagine
    #> 3      C00013     0.04726846     Aspartate
    #> 4      C00004     0.04691644       Alanine
    #> 5      C00017     0.04620232        Lysine
    #> 6      C00014     0.04410009     Glutamate
    #> 7      C00006     0.04113730       Leucine
    #> 8      C00018     0.04088794      Arginine
    #> 9      C00016     0.03934820     Glutamine
    #> 10     C00020     0.03808712 Phenylalanine
    #> 11     C00010     0.03763882     Threonine
    #> 12     C00007     0.03556122    Isoleucine
    #> 13     C00003     0.03494654      Pyruvate
    #> 14     C00011     0.03094167      Cysteine
    #> 15     C00019     0.02768143     Histidine
    #> 16     C00008     0.02736923       Proline
    #> 17     C00012     0.02371414    Methionine
    #> 18     C00009     0.02220371        Serine
    #> 19     C00001     0.02077305       Glucose
    #> 20     C00002     0.01867062       Lactate
    #> 
    #> $rbc_plot

<img src="man/figures/README-example-4.png" width="100%" />

    #> 
    #> $network_plot

<img src="man/figures/README-example-5.png" width="100%" />

    #> 
    #> $heatmap_plot

<img src="man/figures/README-example-6.png" width="100%" />

    #> 
    #> $membership_plot

<img src="man/figures/README-example-7.png" width="100%" />

    #> 
    #> $interaction_plot

<img src="man/figures/README-example-8.png" width="100%" />
