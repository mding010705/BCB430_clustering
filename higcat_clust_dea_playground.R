# for my eyes only...

library(Seurat)
library(anndataR)
library(dplyr)
library(bigstatsr)
library(Matrix)
library(RcppParallel)
library(matrixStats)
library(ggplot2)
library(patchwork)

iso_expression_h5_dir <- commandArgs(trailingOnly = TRUE)[1]
# for each data strat

# iso_expression_h5_dir <- "E:/BCB430/split_isocortex_data/nonPFC_female_adult.h5ad"
SeuratObj <- read_h5ad(iso_expression_h5_dir, as = "Seurat")

# Ensure default assay is RNA
DefaultAssay(SeuratObj) <- "RNA"

# Assign layer X to both layer data and layer counts, since layer x is already normalized
LayerData(SeuratObj, assay = "RNA", layer = "data") <- LayerData(SeuratObj, assay = "RNA", layer = "X")
LayerData(SeuratObj, assay = "RNA", layer = "counts") <- LayerData(SeuratObj, assay = "RNA", layer = "X")

# expression matrix
`expression_matrix` <- LayerData(SeuratObj, assay = "RNA", layer = "data")



# compare PFC_female_adult.h5ad to nonPFC_female_adult.h5ad clusters

# for each cluster in each strata, make a feature vector where each entry is the density of each 
# external cell subclass label

# do pair-wise comparisons (euclidean distance) between these cell type proportions to get 
# pairs of clusters to compare (match nonPFC to PFC, this may result in some nonPFC clusters 
# not being compared to a PFC cluster or one nonPFC cluster being compared to multiple
# PFC clusters)

# for each pairing: note the cell composition and use that to title plots
# merge these 2 clusters and their expression data into one Seurat object with metadata labels of
# PFC and nonPFC
# do DEA
# convert mouse genes to human homologs
# do rrho on significant DE genes compared to cell type matched Maitra et al. results 
# Maitra et al. has 3 astrocytes (Ast), 1 endo (End), 3 oli (Oli), 3 OPC for males
# 1 Ast, 2 Mic, 2 broad mixture (Mix), 4 Oli, 2 OPC for females
# 3 Ast, 1 End, 2 Mic, 1 Mix, 4 Oli, 3 OPC for mixed
# mapping
abc2maitra <- c("Astro-TE NN" = "Ast", "VLMC NN" = "Mix", "Endo NN" = "End",
               "Microglia NN" = "Mic", "Oligo NN" = "Oli", "OPC NN" = "OPC",
               "BAM NN" = "Mix", "SMC NN" = "Mix", "ABC NN" = "Mix", 
               "Peri NN" = "Mix", "Hypendymal NN" = "Mix", "DC NN" = "Mix",
               "Lymphoid NN" = "Mix", "Astro-NT NN" = "Ast")

library(Seurat)
library(anndataR)
library(dplyr)
library(bigstatsr)
library(Matrix)
library(RcppParallel)
library(matrixStats)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(readxl)
library(clusterProfiler)
library(RRHO)
library(GSA)
library(fgsea)
library(org.Mm.eg.db)


# read in ding meta analysis genes, nordic genes

ding_age <- readxl::read_excel("E:/BCB430/ding2015genes.xlsx", 
                               sheet = "Sheet1", range = cell_rows(c(2, NA)))
ding_age <- ding_age[-(1:4), -9]
colnames(ding_age)[5:14] <- c("roP_OC_p",	"roP_OC_q",	"REM_P",	"REM_q", 
                              "mixed_p", "mixed_q",	"female_p", "female_q",
                              "male_p", "male_q")
ding_age[, 5:14] <- as.data.frame(apply(ding_age[, 5:14], MARGIN = 2, 
                                        FUN = function(x){as.numeric(x)}))

nordic <- readxl::read_excel("E:/BCB430/nordic_eo_lo_mdd_genes.xlsx", 
                             sheet = "S2 Genome wide sig genes", range = cell_rows(c(3, NA)))
nordic <- as.data.frame(apply(nordic, MARGIN = 2, FUN = function(x){sub(",", ".", x)}))
nordic$ZSTAT <- as.numeric(nordic$ZSTAT)
nordic$P <- as.numeric(nordic$P)
nordic$p.fdr <- as.numeric(nordic$p.fdr)

nordic_lo <- nordic[grep("Late-onset", nordic$Trait), ]
nordic_eo <- nordic[grep("Early-onset", nordic$Trait), ]

mouse_human <- read.delim(file = "E:/BCB430/human_mouse_homologs.txt")
mouse_human_m <- mouse_human[mouse_human$Common.Organism.Name == "mouse, laboratory", 
                             c("DB.Class.Key", "Symbol")]
mouse_human_h <- mouse_human[mouse_human$Common.Organism.Name == "human", 
                             c("DB.Class.Key", "Symbol")]
mouse2human <- merge(mouse_human_m, mouse_human_h, by = "DB.Class.Key")
colnames(mouse2human) <- c("key", "mouse_gene", "human_gene")

all_genes <- (read.delim(file = "E:/BCB430/all_genes.txt",
                         header = FALSE))[, 1]
all_genes_human <- mouse2human[which(mouse2human$mouse_gene %in% all_genes), ]

ding_mouse <- mouse2human[which(mouse2human$human_gene %in% ding_age$SYMBOL), ]
ding_mouse <- merge(ding_mouse, ding_age, by.x = "human_gene", by.y = "SYMBOL", all.y = FALSE)
nordic_mouse <- mouse2human[which(mouse2human$human_gene %in% nordic$GENE), ]
nordic_eo_mouse <- mouse2human[which(mouse2human$human_gene %in% nordic_eo$GENE), ]
nordic_lo_mouse <- mouse2human[which(mouse2human$human_gene %in% nordic_lo$GENE), ]
nordic_mouse <- merge(nordic_mouse, nordic, by.x = "human_gene", by.y = "GENE", all.y = FALSE)
nordic_eo_mouse <- merge(nordic_eo_mouse, nordic, by.x = "human_gene", by.y = "GENE", all.y = FALSE)
nordic_lo_mouse <- merge(nordic_lo_mouse, nordic, by.x = "human_gene", by.y = "GENE", all.y = FALSE)

maitra_male <- readxl::read_excel("E:/BCB430/maitra_MDD_genes_567.xlsx", 
                                  sheet = "SupplementaryData5", range = cell_rows(c(3, NA)))
maitra_female <- readxl::read_excel("E:/BCB430/maitra_MDD_genes_567.xlsx", 
                                  sheet = "SupplementaryData6", range = cell_rows(c(3, NA)))
maitra_mixed <- readxl::read_excel("E:/BCB430/maitra_MDD_genes_567.xlsx", 
                                  sheet = "SupplementaryData7", range = cell_rows(c(3, NA)))

maitra_female$gene <- gsub("\\..*", "", maitra_female$gene)
maitra_male$gene <- gsub("\\..*", "", maitra_male$gene)
maitra_mixed$gene <- gsub("\\..*", "", maitra_mixed$gene)

maitra <- list(female = maitra_female, male = maitra_male)
gwas <- list(Ding2015 = ding_mouse, Shorter2025_mixed = nordic_mouse, 
             Shorter2025_early = nordic_eo_mouse,
             Shorter2025_late = nordic_lo_mouse)
intersect_maitra_gwas <- list()
for (m in names(maitra)){
  for (g in names(gwas)){
    for (ct in unique(maitra[[m]][[2]])) {
      print(paste(m, g, ct))
      print(intersect(maitra[[m]][maitra[[m]][,2] == ct, ]$gene, gwas[[g]]$human_gene))
      intersect_maitra_gwas[[m]][[g]][[ct]] <- intersect(maitra[[m]][maitra[[m]][,2] == ct, ]$gene, gwas[[g]]$human_gene)
    }
  }
}

male_de <- readRDS(file = "E:/BCB430/cluster/PFC_male_aged.h5ad_PFC_male_adult.h5ad/DE_results.rds")
female_de <- readRDS(file = "E:/BCB430/cluster/PFC_female_aged.h5ad_PFC_female_adult.h5ad/DE_results.rds")
for (i in names(female_de)){
  if (length(female_de[[i]]) == 0){
    next
  }
  female_de[[i]][[1]] <- female_de[[i]][[1]][female_de[[i]][[1]]$p_val_adj < 0.05, ]
}

for (i in names(male_de)){
  if (length(male_de[[i]]) == 0){
    next
  }
  male_de[[i]][[1]] <- male_de[[i]][[1]][male_de[[i]][[1]]$p_val_adj < 0.05, ]
}
hicat <- list(female = female_de, male = male_de)
female_matrix <- matrix(ncol = length(female_de), nrow = length(unique(maitra[["female"]]$cluster_id)))
male_matrix <- matrix(ncol = length(male_de), nrow = length(unique(maitra[["male"]]$cluster_id)))
rownames(female_matrix) <- unique(maitra[["female"]]$cluster_id)
colnames(female_matrix) <- names(female_de)
rownames(male_matrix) <- unique(maitra[["male"]]$cluster_id)
colnames(male_matrix) <- names(male_de)

fm_matrix <- matrix(ncol = length(male_de), nrow = length(unique(maitra[["female"]]$cluster_id)))
mf_matrix <- matrix(ncol = length(female_de), nrow = length(unique(maitra[["male"]]$cluster_id)))
rownames(fm_matrix) <- unique(maitra[["female"]]$cluster_id)
colnames(fm_matrix) <- names(male_de)
rownames(mf_matrix) <- unique(maitra[["male"]]$cluster_id)
colnames(mf_matrix) <- names(female_de)

intersect_maitra_hicat <- list(female_female = as.data.frame(female_matrix), 
                               male_male = as.data.frame(male_matrix),
                               female_male = as.data.frame(fm_matrix), 
                               male_female = as.data.frame(mf_matrix))
maitra_mouse <- lapply(maitra, function(x){
  x_mouse <- mouse2human[which(mouse2human$human_gene %in% x$gene), ]
  x <- merge(x_mouse, x, by.x = "human_gene", by.y = "gene", all.y = FALSE)
  
})
for (m in names(maitra_mouse)){
  for (g in names(hicat)){
    for (ct in unique(maitra_mouse[[m]][[4]])) {
      for (ct2 in names(hicat[[g]])){
        if (length(hicat[[g]][[ct2]]) == 0){
          intersect_maitra_hicat[[paste(m, g, sep = "_")]][ct, ct2] <- ""
          next
        }
        intersect_maitra_hicat[[paste(m, g, sep = "_")]][ct, ct2] <- paste(intersect(maitra_mouse[[m]][maitra_mouse[[m]][,4] == ct, ]$mouse_gene, rownames(hicat[[g]][[ct2]][[1]])), sep = "/", collapse = ",") 
      }
    }
  }
}

intersect_maitra_hicat_count <- list(female_female = as.data.frame(female_matrix), 
                               male_male = as.data.frame(male_matrix),
                               female_male = as.data.frame(fm_matrix), 
                               male_female = as.data.frame(mf_matrix))

for (m in names(intersect_maitra_hicat)){
  for (ct in rownames(intersect_maitra_hicat[[m]])) {
    for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
      intersect_maitra_hicat_count[[m]][ct, ct2] <- length(strsplit(intersect_maitra_hicat[[m]][ct, ct2], ",")[[1]])
      
    }
  }
}

sig_maitra_hicat_count <- list(female_female = as.data.frame(female_matrix), 
                                     male_male = as.data.frame(male_matrix),
                                     female_male = as.data.frame(fm_matrix), 
                                     male_female = as.data.frame(mf_matrix))

for (m in names(intersect_maitra_hicat)){
  for (ct in rownames(intersect_maitra_hicat[[m]])) {
    for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
      if (length(hicat[[gsub(".*_", "", m)]][[ct2]]) == 0){
        sig_maitra_hicat_count[[m]][ct, ct2] <- 1
        next
      }
      sig_maitra_hicat_count[[m]][ct, ct2] <- phyper(q = intersect_maitra_hicat_count[[m]][ct, ct2] - 1,
                                                     m = length(maitra_mouse[[gsub("_.*", "", m)]][maitra_mouse[[gsub("_.*", "", m)]][,4] == ct, ]$mouse_gene),
                                                     n = 32285 - length(maitra_mouse[[gsub("_.*", "", m)]][maitra_mouse[[gsub("_.*", "", m)]][,4] == ct, ]$mouse_gene),
                                                     k = nrow(hicat[[gsub(".*_", "", m)]][[ct2]][[1]]),
                                                     lower.tail = FALSE)
      
    }
  }
}

intersect_maitra_gwas_hicat_ding <- list(female_female = as.data.frame(female_matrix), 
                                         male_male = as.data.frame(male_matrix),
                                         female_male = as.data.frame(fm_matrix), 
                                         male_female = as.data.frame(mf_matrix))
for (m in names(intersect_maitra_hicat)){
    for (ct in rownames(intersect_maitra_hicat[[m]])) {
      for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
        intersect_maitra_gwas_hicat_ding[[m]][ct, ct2] <- paste(intersect(
          strsplit(intersect_maitra_hicat[[m]][ct, ct2], split = ",")[[1]], gwas[["Ding2015"]]$mouse_gene), collapse = ",")
        
      }
    }
}


intersect_maitra_gwas_hicat_shorter <- list(female_female = as.data.frame(female_matrix), 
                                         male_male = as.data.frame(male_matrix),
                                         female_male = as.data.frame(fm_matrix), 
                                         male_female = as.data.frame(mf_matrix))
for (m in names(intersect_maitra_hicat)){
  for (ct in rownames(intersect_maitra_hicat[[m]])) {
    for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
      intersect_maitra_gwas_hicat_shorter[[m]][ct, ct2] <- paste(intersect(
        strsplit(intersect_maitra_hicat[[m]][ct, ct2], split = ",")[[1]], gwas[["Shorter2025_mixed"]]$mouse_gene), collapse = ",")
      
    }
  }
}

intersect_maitra_gwas_hicat_shorter_eo <- list(female_female = as.data.frame(female_matrix), 
                                            male_male = as.data.frame(male_matrix),
                                            female_male = as.data.frame(fm_matrix), 
                                            male_female = as.data.frame(mf_matrix))
for (m in names(intersect_maitra_hicat)){
  for (ct in rownames(intersect_maitra_hicat[[m]])) {
    for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
      intersect_maitra_gwas_hicat_shorter_eo[[m]][ct, ct2] <- paste(intersect(
        strsplit(intersect_maitra_hicat[[m]][ct, ct2], split = ",")[[1]], gwas[["Shorter2025_early"]]$mouse_gene), collapse = ",")
      
    }
  }
}


intersect_maitra_gwas_hicat_shorter_lo <- list(female_female = as.data.frame(female_matrix), 
                                            male_male = as.data.frame(male_matrix),
                                            female_male = as.data.frame(fm_matrix), 
                                            male_female = as.data.frame(mf_matrix))
for (m in names(intersect_maitra_hicat)){
  for (ct in rownames(intersect_maitra_hicat[[m]])) {
    for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
      intersect_maitra_gwas_hicat_shorter_lo[[m]][ct, ct2] <- paste(intersect(
        strsplit(intersect_maitra_hicat[[m]][ct, ct2], split = ",")[[1]], gwas[["Shorter2025_late"]]$mouse_gene), collapse = ",")
      
    }
  }
}

intersect_hicat_shorter_eo <- list(female = as.data.frame(female_matrix[1:2, ]),
                                   male = as.data.frame(male_matrix[1:2, ]))
for (s in c("male", "female")){
  rownames(intersect_hicat_shorter_eo[[s]]) <- c("up", "down")
  for (ct2 in colnames(intersect_maitra_hicat[[paste(s, s, sep = "_")]])) {
    if (length(hicat[[s]][[ct2]]) == 0){
      next
    }
    signed_genes <- hicat[[s]][[ct2]][[1]]
    signed_genes$mouse_gene <- rownames(hicat[[s]][[ct2]][[1]])
    signed_genes <- merge(signed_genes,  gwas[["Shorter2025_early"]], by = "mouse_gene", all = FALSE)
    ups <- signed_genes$mouse_gene[sign(signed_genes$ZSTAT) == sign(signed_genes$avg_log2FC)]
    downs <- signed_genes$mouse_gene[sign(signed_genes$ZSTAT) != sign(signed_genes$avg_log2FC)]
      intersect_hicat_shorter_eo[[s]]["up", ct2] <- paste(ups, collapse = ",")
      intersect_hicat_shorter_eo[[s]]["down", ct2] <- paste(downs, collapse = ",")
      
  }
}



intersect_hicat_shorter_lo <- list(female = as.data.frame(female_matrix[1:2, ]),
                                   male = as.data.frame(male_matrix[1:2, ]))
for (s in c("male", "female")){
  rownames(intersect_hicat_shorter_lo[[s]]) <- c("up", "down")
  for (ct2 in colnames(intersect_maitra_hicat[[paste(s, s, sep = "_")]])) {
    if (length(hicat[[s]][[ct2]]) == 0){
      next
    }
    signed_genes <- hicat[[s]][[ct2]][[1]]
    signed_genes$mouse_gene <- rownames(hicat[[s]][[ct2]][[1]])
    signed_genes <- merge(signed_genes,  gwas[["Shorter2025_late"]], by = "mouse_gene", all = FALSE)
    ups <- signed_genes$mouse_gene[sign(signed_genes$ZSTAT) == sign(signed_genes$avg_log2FC)]
    downs <- signed_genes$mouse_gene[sign(signed_genes$ZSTAT) != sign(signed_genes$avg_log2FC)]
    intersect_hicat_shorter_lo[[s]]["up", ct2] <- paste(ups, collapse = ",")
    intersect_hicat_shorter_lo[[s]]["down", ct2] <- paste(downs, collapse = ",")
    
  }
}


intersect_maitra_gwas_hicat_ding_count <- list(female_female = as.data.frame(female_matrix), 
                                         male_male = as.data.frame(male_matrix),
                                         female_male = as.data.frame(fm_matrix), 
                                         male_female = as.data.frame(mf_matrix))
for (m in names(intersect_maitra_hicat)){
  for (ct in rownames(intersect_maitra_hicat[[m]])) {
    for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
      intersect_maitra_gwas_hicat_ding_count[[m]][ct, ct2] <- length(strsplit(intersect_maitra_gwas_hicat_ding[[m]][ct, ct2], ",")[[1]])
      
    }
  }
}


intersect_maitra_gwas_hicat_shorter_count <- list(female_female = as.data.frame(female_matrix), 
                                            male_male = as.data.frame(male_matrix),
                                            female_male = as.data.frame(fm_matrix), 
                                            male_female = as.data.frame(mf_matrix))
for (m in names(intersect_maitra_hicat)){
  for (ct in rownames(intersect_maitra_hicat[[m]])) {
    for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
      intersect_maitra_gwas_hicat_shorter_count[[m]][ct, ct2] <- length(strsplit(intersect_maitra_gwas_hicat_shorter[[m]][ct, ct2], ",")[[1]])
      
    }
  }
}

library(SuperExactTest)

sig_maitra_gwas_hicat_ding <- list(female_female = as.data.frame(female_matrix), 
                                         male_male = as.data.frame(male_matrix),
                                         female_male = as.data.frame(fm_matrix), 
                                         male_female = as.data.frame(mf_matrix))
for (m in names(intersect_maitra_hicat)){
  for (ct in rownames(intersect_maitra_hicat[[m]])) {
    for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
      if (length(hicat[[gsub(".*_", "", m)]][[ct2]]) == 0){
        sig_maitra_gwas_hicat_ding[[m]][ct, ct2] <- 1.0
        next
      }
      sets <- list(maitra = maitra_mouse[[gsub("_.*", "", m)]][maitra_mouse[[gsub("_.*", "", m)]]$cluster_id == ct, "mouse_gene"], 
                   gwas = ding_mouse$mouse_gene, 
                   hicat = rownames(hicat[[gsub(".*_", "", m)]][[ct2]][[1]]))
      sig_maitra_gwas_hicat_ding[[m]][ct, ct2] <- supertest(sets, n = 32285)[["P.value"]][["111"]]
        
        
    }
  }
}


sig_maitra_gwas_hicat_shorter <- list(female_female = as.data.frame(female_matrix), 
                                            male_male = as.data.frame(male_matrix),
                                            female_male = as.data.frame(fm_matrix), 
                                            male_female = as.data.frame(mf_matrix))
for (m in names(intersect_maitra_hicat)){
  for (ct in rownames(intersect_maitra_hicat[[m]])) {
    for (ct2 in colnames(intersect_maitra_hicat[[m]])) {
      if (length(hicat[[gsub(".*_", "", m)]][[ct2]]) == 0){
        sig_maitra_gwas_hicat_shorter[[m]][ct, ct2] <- 1.0
        next
      }
      sets <- list(maitra = maitra_mouse[[gsub("_.*", "", m)]][maitra_mouse[[gsub("_.*", "", m)]]$cluster_id == ct, "mouse_gene"], 
                   gwas = nordic_mouse$mouse_gene, 
                   hicat = rownames(hicat[[gsub(".*_", "", m)]][[ct2]][[1]]))
      sig_maitra_gwas_hicat_shorter[[m]][ct, ct2] <- supertest(sets, n = 32285)[["P.value"]][["111"]]
      
      
    }
  }
}


sel_cells = c("Ast-TE NN", "OPC NN", "Oligo NN",
              "COP NN_1", "Microglia NN", "Peri NN")
library(tibble)

cell_type_triple_overlap_heatmap <- function(odf = intersect_maitra_gwas_hicat_ding_count, 
                                             pdf = sig_maitra_gwas_hicat_ding, 
                                             stat = "female_female", 
                                             study = "Ding2015",
                                             comp = "PFC_female_aged vs PFC_female_adult"){
  df <- odf[[stat]]
  
  row_hclust <- hclust(dist(df, method = "euclidean"), method = "ward.D")
  
  # Get the order of the clustered items
  row_order <- row_hclust$order
  
  col_hclust <- hclust(dist(t(df), method = "euclidean"), method = "ward.D")
  
  # Get the order of the clustered items
  col_order <- col_hclust$order
  df_long <- df %>% rownames_to_column(var = "maitra_cell_type") %>%
    pivot_longer(
      cols = colnames(df), # Select columns to reshape (can also use column numbers or names, e.g., c(3:5))
      names_to = "hicat_cell_type",   # Name of the new column to store old column names
      values_to = "overlap_size"        # Name of the new column to store the values
    )
  df_longp <- pdf[[stat]] %>% rownames_to_column(var = "maitra_cell_type") %>%
    pivot_longer(
      cols = colnames(df), # Select columns to reshape (can also use column numbers or names, e.g., c(3:5))
      names_to = "hicat_cell_type",   # Name of the new column to store old column names
      values_to = "pval"        # Name of the new column to store the values
    )
  df_long <- merge(df_long, df_longp, by = c("maitra_cell_type", "hicat_cell_type"), all = FALSE)
  df_long$hicat_cell_type <- factor(df_long$hicat_cell_type, levels = colnames(df)[col_order])
  df_long$maitra_cell_type <- factor(df_long$maitra_cell_type, levels = rownames(df)[row_order])
  df_long$pval <- p.adjust(df_long$pval, method = "fdr")
  df_long$pstar <- ifelse(df_long[, "pval"] < 0.0001, "***",
                     ifelse(df_long[, "pval"] < 0.001, "**",
                            ifelse(df_long[, "pval"] < 0.05, "*", "")))
  leg_name <- "overlap size"
  
  p <- (ggplot(df_long, aes(hicat_cell_type, maitra_cell_type)) +
            geom_tile(aes(fill = overlap_size)) +
            scale_fill_gradient(low = "white", high = "red", name = leg_name) +
            ggtitle(paste("Maitra2023", study, "Hicat DE overlap", stat, sep = ", "),
                    subtitle = comp) +
            geom_text(aes(label = pstar), color = "green", size=3, vjust = -0.5) +
            geom_text(aes(label = overlap_size), vjust = 0.7, size = 4, color = "black") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  panel.background = element_rect(fill = "grey", colour = NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
  
  ggsave(plot = p, filename = paste0("E:/BCB430/images/overlap_hm_", 
                                     stat, "_", study, "_", comp, "_3.png"), height = 8, width = 17, dpi = 600)
  return(p)
}

cell_type_triple_overlap_heatmap(odf = intersect_maitra_gwas_hicat_ding_count, 
                                             pdf = sig_maitra_gwas_hicat_ding, 
                                             stat = "female_female", 
                                             study = "Ding2015",
                                             comp = "PFC_female_aged vs PFC_female_adult")
cell_type_triple_overlap_heatmap(odf = intersect_maitra_gwas_hicat_ding_count, 
                                 pdf = sig_maitra_gwas_hicat_ding, 
                                 stat = "male_female", 
                                 study = "Ding2015",
                                 comp = "PFC_female_aged vs PFC_female_adult")
cell_type_triple_overlap_heatmap(odf = intersect_maitra_gwas_hicat_shorter_count, 
                                 pdf = sig_maitra_gwas_hicat_shorter, 
                                 stat = "female_female", 
                                 study = "Shorter2025",
                                 comp = "PFC_female_aged vs PFC_female_adult")
cell_type_triple_overlap_heatmap(odf = intersect_maitra_gwas_hicat_shorter_count, 
                                 pdf = sig_maitra_gwas_hicat_shorter, 
                                 stat = "male_female", 
                                 study = "Shorter2025",
                                 comp = "PFC_female_aged vs PFC_female_adult")

cell_type_triple_overlap_heatmap(odf = intersect_maitra_gwas_hicat_ding_count, 
                                 pdf = sig_maitra_gwas_hicat_ding, 
                                 stat = "male_male", 
                                 study = "Ding2015",
                                 comp = "PFC_male_aged vs PFC_male_adult")
cell_type_triple_overlap_heatmap(odf = intersect_maitra_gwas_hicat_ding_count, 
                                 pdf = sig_maitra_gwas_hicat_ding, 
                                 stat = "female_male", 
                                 study = "Ding2015",
                                 comp = "PFC_male_aged vs PFC_male_adult")
cell_type_triple_overlap_heatmap(odf = intersect_maitra_gwas_hicat_shorter_count, 
                                 pdf = sig_maitra_gwas_hicat_shorter, 
                                 stat = "male_male", 
                                 study = "Shorter2025",
                                 comp = "PFC_male_aged vs PFC_male_adult")
cell_type_triple_overlap_heatmap(odf = intersect_maitra_gwas_hicat_shorter_count, 
                                 pdf = sig_maitra_gwas_hicat_shorter, 
                                 stat = "female_male", 
                                 study = "Shorter2025",
                                 comp = "PFC_male_aged vs PFC_male_adult")


cell_type_maitra_hicat_overlap_heatmap <- function(odf = intersect_maitra_hicat_count, 
                                             pdf = sig_maitra_hicat_count, 
                                             stat = "female_female", 
                                             comp = "PFC_female_aged vs PFC_female_adult"){
  df <- odf[[stat]]
  
  row_hclust <- hclust(dist(df, method = "euclidean"), method = "ward.D")
  
  # Get the order of the clustered items
  row_order <- row_hclust$order
  
  col_hclust <- hclust(dist(t(df), method = "euclidean"), method = "ward.D")
  
  # Get the order of the clustered items
  col_order <- col_hclust$order
  df_long <- df %>% rownames_to_column(var = "maitra_cell_type") %>%
    pivot_longer(
      cols = colnames(df), # Select columns to reshape (can also use column numbers or names, e.g., c(3:5))
      names_to = "hicat_cell_type",   # Name of the new column to store old column names
      values_to = "overlap_size"        # Name of the new column to store the values
    )
  df_longp <- pdf[[stat]] %>% rownames_to_column(var = "maitra_cell_type") %>%
    pivot_longer(
      cols = colnames(df), # Select columns to reshape (can also use column numbers or names, e.g., c(3:5))
      names_to = "hicat_cell_type",   # Name of the new column to store old column names
      values_to = "pval"        # Name of the new column to store the values
    )
  df_long <- merge(df_long, df_longp, by = c("maitra_cell_type", "hicat_cell_type"), all = FALSE)
  df_long$hicat_cell_type <- factor(df_long$hicat_cell_type, levels = colnames(df)[col_order])
  df_long$maitra_cell_type <- factor(df_long$maitra_cell_type, levels = rownames(df)[row_order])
  df_long$pval <- p.adjust(df_long$pval, method = "fdr")
  df_long$pstar <- ifelse(df_long[, "pval"] < 0.0001, "***",
                          ifelse(df_long[, "pval"] < 0.001, "**",
                                 ifelse(df_long[, "pval"] < 0.05, "*", "")))
  leg_name <- "overlap size"
  
  p <- (ggplot(df_long, aes(hicat_cell_type, maitra_cell_type)) +
          geom_tile(aes(fill = overlap_size)) +
          scale_fill_gradient(low = "white", high = "red", name = leg_name) +
          ggtitle(paste("Maitra2023", "Hicat DE overlap", stat, sep = ", "),
                  subtitle = comp) +
          geom_text(aes(label = pstar), color = "green", size=3, vjust = -0.5) +
          geom_text(aes(label = overlap_size), vjust = 0.7, size = 4, color = "black") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                panel.background = element_rect(fill = "grey", colour = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()))
  
  ggsave(plot = p, filename = paste0("E:/BCB430/images/overlap_hm_", 
                                     stat, "_", "_", comp, "_2.png"), height = 8, width = 17, dpi = 600)
  return(p)
}

cell_type_maitra_hicat_overlap_heatmap(odf = intersect_maitra_hicat_count, 
                                                   pdf = sig_maitra_hicat_count, 
                                                   stat = "female_female", 
                                                   comp = "PFC_female_aged vs PFC_female_adult")
cell_type_maitra_hicat_overlap_heatmap(odf = intersect_maitra_hicat_count, 
                                       pdf = sig_maitra_hicat_count, 
                                       stat = "male_female", 
                                       comp = "PFC_female_aged vs PFC_female_adult")

cell_type_maitra_hicat_overlap_heatmap(odf = intersect_maitra_hicat_count, 
                                       pdf = sig_maitra_hicat_count, 
                                       stat = "male_male", 
                                       comp = "PFC_male_aged vs PFC_male_adult")
cell_type_maitra_hicat_overlap_heatmap(odf = intersect_maitra_hicat_count, 
                                       pdf = sig_maitra_hicat_count, 
                                       stat = "female_male", 
                                       comp = "PFC_male_aged vs PFC_male_adult")


library(UpSetR)
library(grid)

listInput <- lapply(female_de, function(x){
  if (length(x) == 0){
    return(NULL)
  }
  rownames(x[[1]])
})

upset(fromList(listInput), nsets = length(listInput), order.by="freq")
grid.text("PFC_female_aged vs PFC_female_adult \ncell type DE gene overlaps", x = 0.65, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))

listInput <- lapply(male_de, function(x){
  if (length(x) == 0){
    return(NULL)
  }
  rownames(x[[1]])
})

upset(fromList(listInput), nsets = length(listInput), order.by="freq")
grid.text("PFC_male_aged vs PFC_male_adult \ncell type DE gene overlaps", x = 0.65, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))



mat <- matrix(nrow = length(hicat$male), ncol = length(hicat$female))
colnames(mat) <- names(hicat$female)
rownames(mat) <- names(hicat$male)
intersect_male_female <- as.data.frame(mat)

for (ct in rownames(mat)){
  for (ct2 in colnames(mat)){
    if (length(hicat$male[[ct]]) == 0 | length(hicat$female[[ct2]]) == 0){
      intersect_male_female[ct, ct2] <- NA
      next
    }
    intersect_male_female[ct, ct2] <- length(intersect(rownames(hicat$male[[ct]][[1]]),
                                                rownames(hicat$female[[ct2]][[1]])))
  }
}


sig_male_female <- as.data.frame(mat)

for (ct in rownames(mat)){
  for (ct2 in colnames(mat)){
    if (length(hicat$male[[ct]]) == 0 | length(hicat$female[[ct2]]) == 0){
      sig_male_female[ct, ct2] <- 1
      next
    }
    sig_male_female[ct, ct2] <- phyper(q = intersect_male_female[ct, ct2] - 1,
                                       m = nrow(hicat$male[[ct]][[1]]),
                                       n = 32285 - nrow(hicat$male[[ct]][[1]]),
                                       k = nrow(hicat$female[[ct2]][[1]]),
                                       lower.tail = FALSE)
  }
}

intersect_male_female <- intersect_male_female[-which(is.na(intersect_male_female[, 1])), -which(is.na(intersect_male_female[1, ]))]
row_hclust <- hclust(dist(intersect_male_female, method = "euclidean"), method = "ward.D")

# Get the order of the clustered items
row_order <- row_hclust$order

col_hclust <- hclust(dist(t(intersect_male_female), method = "euclidean"), method = "ward.D")

# Get the order of the clustered items
col_order <- col_hclust$order
df_long <- intersect_male_female %>% rownames_to_column(var = "male_cell_type") %>%
  pivot_longer(
    cols = colnames(intersect_male_female), # Select columns to reshape (can also use column numbers or names, e.g., c(3:5))
    names_to = "female_cell_type",   # Name of the new column to store old column names
    values_to = "overlap_size"        # Name of the new column to store the values
  )
df_longp <- sig_male_female %>% rownames_to_column(var = "male_cell_type") %>%
  pivot_longer(
    cols = colnames(sig_male_female), # Select columns to reshape (can also use column numbers or names, e.g., c(3:5))
    names_to = "female_cell_type",   # Name of the new column to store old column names
    values_to = "pval"        # Name of the new column to store the values
  )
df_long <- merge(df_long, df_longp, by = c("male_cell_type", "female_cell_type"), all = FALSE)
df_long$female_cell_type <- factor(df_long$female_cell_type, levels = colnames(intersect_male_female)[col_order])
df_long$male_cell_type <- factor(df_long$male_cell_type, levels = rownames(intersect_male_female)[row_order])
df_long$pval <- p.adjust(df_long$pval, method = "fdr")
df_long$pstar <- ifelse(df_long[, "pval"] < 0.0001, "***",
                        ifelse(df_long[, "pval"] < 0.001, "**",
                               ifelse(df_long[, "pval"] < 0.05, "*", "")))
leg_name <- "overlap size"

p <- (ggplot(df_long, aes(female_cell_type, male_cell_type)) +
        geom_tile(aes(fill = overlap_size)) +
        scale_fill_gradient(low = "white", high = "red", name = leg_name) +
        ggtitle("Male, Female DE gene overlap",
                subtitle = "PFC_male_aged vs PFC_male_adult, \nPFC_female_aged vs PFC_female_adult") +
        geom_text(aes(label = pstar), color = "green", size=3, vjust = -0.5) +
        geom_text(aes(label = overlap_size), vjust = 0.7, size = 4, color = "black") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              panel.background = element_rect(fill = "grey", colour = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))
p
ggsave(plot = p, filename = paste0("E:/BCB430/images/overlap_hm_male_female_de.png"), height = 8, width = 17, dpi = 600)
# could also just compare each inter-cluster DE list to each cell type specific MDD signature
# (compare to female, male, and non-sex specific results)
# (do both one sided and two sided to determine directionality + significance of similarity)

# do rrho test with GWAS meta analysis from Ding et al.

