# do GO:BP enrichment, write to a file, and import that file into cytoscape


library(gprofiler2)

gmt_file <- "C:/Users/mding/Desktop/BCB430/Mouse_GOBP_AllPathways_noPFOCR_no_GO_iea_September_01_2025_symbol.gmt"
output_dir <- "C:/Users/mding/Desktop/BCB430/"
# 1. Upload your custom GMT file to get a unique token
my_gmt_token <- upload_GMT_file(gmtfile = gmt_file)

# 2. Run gost using the token
# Ensure your query genes match the IDs in your GMT file
mdd_enrich_up <- readRDS("C:/Users/mding/Desktop/BCB430/cluster/PFC_male_aged.h5ad_PFC_male_adult.h5ad/mdd_enrich_up.rds")
mdd_enrich_down <- readRDS("C:/Users/mding/Desktop/BCB430/cluster/PFC_male_aged.h5ad_PFC_male_adult.h5ad/mdd_enrich_down.rds")

# only use our cell types of interest
results <- gost(query = list(increase_oligo = mdd_enrich_up[["age_cat"]][["Oligo NN"]][["Ding2015"]]@gene, 
                             decrease_oligo = mdd_enrich_down[["age_cat"]][["Oligo NN"]][["Ding2015"]]@gene,
                             increase_opc = mdd_enrich_up[["age_cat"]][["OPC NN"]][["Ding2015"]]@gene, 
                             decrease_opc = mdd_enrich_down[["age_cat"]][["OPC NN"]][["Ding2015"]]@gene,
                             increase_cop = mdd_enrich_up[["age_cat"]][["COP NN_1"]][["Ding2015"]]@gene, 
                             decrease_cop = mdd_enrich_down[["age_cat"]][["COP NN_1"]][["Ding2015"]]@gene,
                             increase_ast = mdd_enrich_up[["age_cat"]][["Astro-TE NN"]][["Ding2015"]]@gene, 
                             decrease_ast = mdd_enrich_down[["age_cat"]][["Astro-TE NN"]][["Ding2015"]]@gene,
                             increase_mic = mdd_enrich_up[["age_cat"]][["Microglia NN"]][["Ding2015"]]@gene, 
                             decrease_mic = mdd_enrich_down[["age_cat"]][["Microglia NN"]][["Ding2015"]]@gene,
                             increase_peri = mdd_enrich_up[["age_cat"]][["Peri NN"]][["Ding2015"]]@gene, 
                             decrease_peri = mdd_enrich_down[["age_cat"]][["Peri NN"]][["Ding2015"]]@gene),
                organism = my_gmt_token, correction_method = "fdr", evcodes = TRUE)

# 3. View results
head(results$result)
a <- results$result
a_clean <- data.frame(lapply(a, as.character), stringsAsFactors = FALSE)
a_clean$term_id <- gsub(" ", "_", a_clean$term_id)
# Extract the result data frame
gem_data <- results$result

# Format for GEM

gem_formatted <- gem_data[, c("term_id", "term_name", "p_value", "p_value", "query", "intersection")]

# Rename columns to meet GEM specifications
colnames(gem_formatted) <- c("name", "descr", "pvalue", "FDR", "Cell_Type", "Genes")

gem_formatted$ES <- ifelse(grepl("decrease", gem_formatted$Cell_Type), -1, 1)
gem_formatted$NES <- gem_formatted$ES
gem_formatted$Cell_Type <- gsub(".*_", "", gem_formatted$Cell_Type)
# Reorder columns: Name, Description, p.value, FDR, Phenotype, Genes
gem_formatted <- gem_formatted[, c("name", "descr", "pvalue", "FDR", "ES", "NES", "Genes", "Cell_Type")]
data.table::fwrite(gem_formatted, 
            file = paste0(output_dir, "male_age_gost_cts.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

library(RCy3)
# connect to the running cytoscape session
cytoscapePing()
em_results_filename <- paste0(output_dir, "male_age_gost_cts.txt")
#em_results_filename <- file.path("C:/Users/mding/Desktop/BCB420/Megan_Ding/A2/enr_results.txt")

# min q val = 0.43... eek
# choose q val threshold to be 0.6 
em_command <-  paste(
  'enrichmentmap build analysisType="generic" ',
   'gmtFile=', gmt_file,
   'pvalue=',"0.05",
  'similaritycutoff=',"0.25",
  'coefficients=',"JACCARD",
  'enrichmentsDataset1=',em_results_filename,
  sep="\t")

# get suid of newly created network.
em_network_suid <- commandsRun(em_command)

renameNetwork("Male_MDD_gost_enrichmentmap", network=as.numeric(em_network_suid))


# female =============================

mdd_enrich_up <- readRDS("C:/Users/mding/Desktop/BCB430/cluster/PFC_female_aged.h5ad_PFC_female_adult.h5ad/mdd_enrich_up.rds")
mdd_enrich_down <- readRDS("C:/Users/mding/Desktop/BCB430/cluster/PFC_female_aged.h5ad_PFC_female_adult.h5ad/mdd_enrich_down.rds")

results <- gost(query = list(increase_oligo = mdd_enrich_up[["age_cat"]][["Oligo NN"]][["Ding2015"]]@gene, 
                             decrease_oligo = mdd_enrich_down[["age_cat"]][["Oligo NN"]][["Ding2015"]]@gene,
                             increase_opc = mdd_enrich_up[["age_cat"]][["OPC NN"]][["Ding2015"]]@gene, 
                             decrease_opc = mdd_enrich_down[["age_cat"]][["OPC NN"]][["Ding2015"]]@gene,
                             increase_cop = mdd_enrich_up[["age_cat"]][["COP NN_1"]][["Ding2015"]]@gene, 
                             # decrease_cop = mdd_enrich_down[["age_cat"]][["COP NN_1"]][["Ding2015"]]@gene, <- no COPs here
                             increase_ast = mdd_enrich_up[["age_cat"]][["Astro-TE NN"]][["Ding2015"]]@gene, 
                             decrease_ast = mdd_enrich_down[["age_cat"]][["Astro-TE NN"]][["Ding2015"]]@gene,
                             increase_mic = mdd_enrich_up[["age_cat"]][["Microglia NN"]][["Ding2015"]]@gene, 
                             decrease_mic = mdd_enrich_down[["age_cat"]][["Microglia NN"]][["Ding2015"]]@gene,
                             increase_peri = mdd_enrich_up[["age_cat"]][["Peri NN"]][["Ding2015"]]@gene, 
                             decrease_peri = mdd_enrich_down[["age_cat"]][["Peri NN"]][["Ding2015"]]@gene),
                organism = my_gmt_token, correction_method = "fdr", evcodes = TRUE)

# 3. View results
head(results$result)

gem_data <- results$result

# Format for GEM

gem_formatted <- gem_data[, c("term_id", "term_name", "p_value", "p_value", "query", "intersection")]

# Rename columns to meet GEM specifications
colnames(gem_formatted) <- c("name", "descr", "pvalue", "FDR", "Cell_Type", "Genes")

gem_formatted$ES <- ifelse(grepl("decrease", gem_formatted$Cell_Type), -1, 1)
gem_formatted$NES <- gem_formatted$ES
gem_formatted$Cell_Type <- gsub(".*_", "", gem_formatted$Cell_Type)
# Reorder columns: Name, Description, p.value, FDR, Phenotype, Genes
gem_formatted <- gem_formatted[, c("name", "descr", "pvalue", "FDR", "ES", "NES", "Genes", "Cell_Type")]
data.table::fwrite(gem_formatted, 
                   file = paste0(output_dir, "female_age_gost_cts.txt"), 
                   sep = "\t", quote = FALSE, row.names = FALSE)

library(RCy3)
# connect to the running cytoscape session
cytoscapePing()
em_results_filename <- paste0(output_dir, "female_age_gost_cts.txt")
#em_results_filename <- file.path("C:/Users/mding/Desktop/BCB420/Megan_Ding/A2/enr_results.txt")

# min q val = 0.43... eek
# choose q val threshold to be 0.6 
em_command <-  paste(
  'enrichmentmap build analysisType="generic" ',
  'gmtFile=', gmt_file,
  'pvalue=',"0.05",
  'similaritycutoff=',"0.25",
  'coefficients=',"JACCARD",
  'enrichmentsDataset1=',em_results_filename,
  sep="\t")

# get suid of newly created network.
em_network_suid <- commandsRun(em_command)

renameNetwork("Female_MDD_gost_enrichmentmap", network=as.numeric(em_network_suid))
