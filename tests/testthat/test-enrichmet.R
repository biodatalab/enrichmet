library(testthat)

inputMetabolites <- c(
    "C00209", "C00249", "C06424", "C01727", "C04025", "C00366",
    "C00385", "C08261", "C00079", "C00506", "C00093", "C02979",
    "C00106", "C01384", "C00219", "C00242", "C19806", "C00105",
    "C00606", "C00144", "C01216", "C06429",
    "C00475", "C00386", "C01046", "C02180", "C00082", "C00299",
    "C00051", "C16513", "C00387", "C00022", "C00026",
    "C01062", "C00418", "C00295", "C01762", "C00525", "C00074",
    "C00491", "C02067", "C00015", "C14829", "C00519",
    "C06428", "C05472", "C00328", "C16358", "C16357", "C16353",
    "C02632", "C00246", "C05842", "C00785", "C03739",
    "C00212", "C01181", "C00140", "C05021", "C02862", "C02721",
    "C01104", "C10172", "C00380", "C01262", "C05122",
    "C03736", "C00599", "C00446", "C00275", "C00103", "C00881",
    "C15587", "C00361", "C00429", "C01042", "C11439",
    "C05637", "C00127", "C00035", "C01909", "C00836", "C00847",
    "C04294", "C00043", "C01481", "C00301", "C05635",
    "C02140", "C05488", "C00570", "C01551", "C00092", "C03139"
)

PathwayVsMetabolites <- fetch_kegg_pathway_metabolites()

kegg_lookup <- fetch_kegg_compound_lookup()

example_path <- system.file(
"extdata", "summary_stat.csv",
package = "enrichmet"
)

if (example_path == "") {
    stop("Example file 'summary_stat.csv' not found in inst/extdata/")
}

example_data <- read.csv(example_path, stringsAsFactors = FALSE)
# Minimal Reactome test data.
# Reactome interaction requires KEGG and Reaction columns.
reactome_df <- data.frame(
    KEGG = c(
        "C00209", "C00249",
        "C00209", "C06424",
        "C00249", "C06424"
    ),
    Reaction = c(
        "R_TEST_1", "R_TEST_1",
        "R_TEST_2", "R_TEST_2",
        "R_TEST_3", "R_TEST_3"
    ),
    stringsAsFactors = FALSE
)


test_that("enrichmet returns expected output format", {
    
    result <- enrichmet(
        inputMetabolites = inputMetabolites,
        PathwayVsMetabolites = PathwayVsMetabolites,
        example_data = example_data,
        kegg_lookup = kegg_lookup,
        reactome_df = reactome_df,
        top_n = 20,
        p_value_cutoff = 0.05,
        network_top_n = 10,
        heatmap_top_n = 20,
        membership_top_n = 20,
        min_pathway_occurrence = 2,
        min_metabolite_occurrence = 1
    )
    
    expect_type(result, "list")
    
    expect_named(
        result,
        c(
            "input_metabolites_used",
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
            "interaction_plot"
        )
    )
    
    expect_true(
        length(result$input_metabolites_used) > 0
    )
    
    expect_true(
        is.data.frame(result$pathway_enrichment_all)
    )
    
    expect_true(
        is.data.frame(result$pathway_enrichment_results)
    )
    
    expect_true(
        is.data.frame(result$gsea_results)
    )
    
    expect_true(
        is.data.frame(result$metabolite_centrality)
    )
    
    expect_true(
        !is.null(result$pathway_plot)
    )
    
    expect_true(
        !is.null(result$impact_plot)
    )
    
    expect_true(
        !is.null(result$gsea_plot)
    )
    
    expect_true(
        !is.null(result$network_plot)
    )
    
    expect_true(
        !is.null(result$heatmap_plot)
    )
    
    expect_true(
        !is.null(result$membership_plot)
    )
    
    expect_true(
        !is.null(result$interaction_plot)
    )
})