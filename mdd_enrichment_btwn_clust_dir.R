# intersect DE genes with MDD gene lists, evaluate with hypergeometric tests, RRHO tests,
# and do GO:BP enrichment on the genes in this overlap, but this time, 
# separate analyses by increase and decrease depressive burden

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

de_dir <- commandArgs(trailingOnly = TRUE)[1]
updown <- commandArgs(trailingOnly = TRUE)[2]
matches <- gregexpr("h5ad", de_dir, fixed = TRUE)
counts <- sum(lengths(matches))
btwn <- counts == 2
btwn_clust_de <- readRDS(file = paste0(de_dir, "/DE_results.rds"))
# read in ding meta analysis genes, nordic genes

ding_age <- readxl::read_excel("C:/Users/megan_ding/Desktop/BCB430/ding2015genes.xlsx", 
                                       sheet = "Sheet1", range = cell_rows(c(2, NA)))
ding_age <- ding_age[-(1:4), -9]
colnames(ding_age)[5:14] <- c("roP_OC_p",	"roP_OC_q",	"REM_P",	"REM_q", 
                              "mixed_p", "mixed_q",	"female_p", "female_q",
                              "male_p", "male_q")
ding_age[, 5:14] <- as.data.frame(apply(ding_age[, 5:14], MARGIN = 2, 
                                        FUN = function(x){as.numeric(x)}))

nordic <- readxl::read_excel("C:/Users/megan_ding/Desktop/BCB430/nordic_eo_lo_mdd_genes.xlsx", 
                                sheet = "S2 Genome wide sig genes", range = cell_rows(c(3, NA)))
nordic <- as.data.frame(apply(nordic, MARGIN = 2, FUN = function(x){sub(",", ".", x)}))
nordic$ZSTAT <- as.numeric(nordic$ZSTAT)
nordic$P <- as.numeric(nordic$P)
nordic$p.fdr <- as.numeric(nordic$p.fdr)

nordic_lo <- nordic[grep("Late-onset", nordic$Trait), ]
nordic_eo <- nordic[grep("Early-onset", nordic$Trait), ]

nordic_snps <- readxl::read_excel("C:/Users/megan_ding/Desktop/BCB430/nordic_eo_lo_mdd_genes.xlsx", 
                                  sheet = "S1 All significant SNPs", range = cell_rows(c(3, NA)))

directory_path <- "C:/Users/megan_ding/Desktop/BCB430/eQTLs"

# Get all file names with their full paths
all_files <- list.files(
  path = directory_path, 
  full.names = TRUE, 
  recursive = FALSE # Set to TRUE to include subdirectories
)

eqtl_list <- lapply(all_files, read.delim)
names(eqtl_list) <- basename(all_files)

# TODO: Merge LD and sig eQTL results
# for (i in 1:6){
#   new_df <- data.frame(eqtl_list[[i + 6]]$Gene.name,
#                        NA, NA,
#                        eqtl_list[[i + 6]]$sc.eQTL.SNP
#                        )
# }
# sapply(eqtl_list, function(x){sum(x$variantId %in% nordic_snps$SNP)})

# simple hypergeometric test
# RRHO by p-value
# per cell type and Ding all, Ding female, Ding male and nordic all, nordic eo, nordic lo
mouse_human <- read.delim(file = "C:/Users/megan_ding/Desktop/BCB430/human_mouse_homologs.txt")
mouse_human_m <- mouse_human[mouse_human$Common.Organism.Name == "mouse, laboratory", 
                             c("DB.Class.Key", "Symbol")]
mouse_human_h <- mouse_human[mouse_human$Common.Organism.Name == "human", 
                             c("DB.Class.Key", "Symbol")]
mouse2human <- merge(mouse_human_m, mouse_human_h, by = "DB.Class.Key")
colnames(mouse2human) <- c("key", "mouse_gene", "human_gene")

all_genes <- (read.delim(file = "C:/Users/megan_ding/Desktop/BCB430/all_genes.txt",
                        header = FALSE))[, 1]
all_genes_human <- mouse2human[which(mouse2human$mouse_gene %in% all_genes), ]

gene_enrichment <- function(gene_set, organism = org.Mm.eg.db,
                            background = NULL,
                            p_adj_method = "fdr",
                            min_gene_set_size = 5,
                            max_gene_set_size = 800,
                            p_cutoff = 0.05, ont_type = "BP",
                            gene_nom = "SYMBOL"){
  
  if(is.null(background)){
      enrich_results <-
        clusterProfiler::enrichGO(
          gene = gene_set,
          OrgDb = organism,
          keyType = gene_nom,
          ont = ont_type,
          pAdjustMethod = p_adj_method,
          qvalueCutoff = p_cutoff,
          readable = F)
    } else {
      enrich_results <-
        clusterProfiler::enrichGO(
          gene = gene_set,
          universe = as.character(sapply(as.vector(background),
                                         as.character)),
          OrgDb = organism,
          keyType = gene_nom,
          ont = ont_type,
          pAdjustMethod = p_adj_method,
          qvalueCutoff = p_cutoff,
          readable = F)
    }
  return(enrich_results)
}

ding_mouse <- mouse2human[which(mouse2human$human_gene %in% ding_age$SYMBOL), ]
ding_mouse <- merge(ding_mouse, ding_age, by.x = "human_gene", by.y = "SYMBOL", all.y = FALSE)
nordic_mouse <- mouse2human[which(mouse2human$human_gene %in% nordic$GENE), ]
nordic_eo_mouse <- mouse2human[which(mouse2human$human_gene %in% nordic_eo$GENE), ]
nordic_lo_mouse <- mouse2human[which(mouse2human$human_gene %in% nordic_lo$GENE), ]
nordic_mouse <- merge(nordic_mouse, nordic, by.x = "human_gene", by.y = "GENE", all.y = FALSE)
nordic_eo_mouse <- merge(nordic_eo_mouse, nordic, by.x = "human_gene", by.y = "GENE", all.y = FALSE)
nordic_lo_mouse <- merge(nordic_lo_mouse, nordic, by.x = "human_gene", by.y = "GENE", all.y = FALSE)


mdd_enrich <- list()
mdd_overlaps <- data.frame(cell_type = NA, contrast = NA, between = NA,
                           ding_overlap_size = NA,
                           nordic_mixed_overlap_size = NA,
                           nordic_early_overlap_size = NA,
                           nordic_late_overlap_size = NA)
if (updown == "up"){
  thresh <- 1
} else {
  thresh <- -1
}
for(ct in names(btwn_clust_de)){
  for (feat in names(btwn_clust_de[[ct]])){
    des <- btwn_clust_de[[ct]][[feat]]
    des_sig <- des[des$p_val_adj < 0.05, ] 
    des_sig$gene <- rownames(des_sig)
    ding_overlap <- merge(des_sig, ding_mouse, by.x = "gene", by.y = "mouse_gene", all = FALSE)
    degs_human <- mouse2human[which(mouse2human$mouse_gene %in% rownames(des_sig)), ]
    ding_overlap <- ding_overlap[sign(ding_overlap$avg_log2FC) * sign(ding_overlap$`Average Effect size`) == thresh, ]
    
    nordic_overlap_mixed <- merge(des_sig, nordic_mouse, by.x = "gene", by.y = "mouse_gene", all = FALSE)
    nordic_overlap_mixed <- nordic_overlap_mixed[sign(nordic_overlap_mixed$avg_log2FC) * sign(nordic_overlap_mixed$ZSTAT) == thresh, ]
    
    nordic_overlap_early <- merge(des_sig, nordic_eo_mouse, by.x = "gene", by.y = "mouse_gene", all = FALSE)
    nordic_overlap_early <- nordic_overlap_early[sign(nordic_overlap_early$avg_log2FC) * sign(nordic_overlap_early$ZSTAT) == thresh, ]
    
    nordic_overlap_late <- merge(des_sig, nordic_lo_mouse, by.x = "gene", by.y = "mouse_gene", all = FALSE)
    nordic_overlap_late <- nordic_overlap_late[sign(nordic_overlap_late$avg_log2FC) * sign(nordic_overlap_late$ZSTAT) == thresh, ]
    
    gene_enrichment_ding <- gene_enrichment(gene_set = ding_overlap$gene, 
                                            background = unique(all_genes_human$mouse_gene))
    gene_enrichment_nordic_mixed <- gene_enrichment(gene_set = nordic_overlap_mixed, 
                                            background = unique(all_genes_human$mouse_gene))
    gene_enrichment_nordic_early <- gene_enrichment(gene_set = nordic_overlap_early, 
                                            background = unique(all_genes_human$mouse_gene))
    gene_enrichment_nordic_late <- gene_enrichment(gene_set = nordic_overlap_late, 
                                            background = unique(all_genes_human$mouse_gene))
    gene_enrichs <- list(Ding2015 = gene_enrichment_ding, Shorter2025_mixed = gene_enrichment_nordic_mixed,
                         Shorter2025_early = gene_enrichment_nordic_early, 
                         Shorter2025_late = gene_enrichment_nordic_late)
  
    mdd_enrich[[feat]][[ct]] <- gene_enrichs
    
    new_df <- data.frame(cell_type = ct, contrast = feat, between = btwn,
                         ding_overlap_size = ifelse(is.null(nrow(ding_overlap)), 0, nrow(ding_overlap)),
                         nordic_mixed_overlap_size = ifelse(is.null(nrow(nordic_overlap_mixed)), 0, nrow(nordic_overlap_mixed)),
                         nordic_early_overlap_size = ifelse(is.null(nrow(nordic_overlap_early)), 0, nrow(nordic_overlap_early)),
                         nordic_late_overlap_size = ifelse(is.null(nrow(nordic_overlap_late)), 0, nrow(nordic_overlap_late)))
    mdd_overlaps <- rbind(mdd_overlaps, new_df)
  }
}
mdd_overlaps <- mdd_overlaps[-1, ]

saveRDS(mdd_enrich, file = paste0(de_dir, "/", "mdd_enrich_", updown, ".rds"))
saveRDS(mdd_overlaps, file = paste0(de_dir, "/", "mdd_overlaps_", updown, ".rds"))

goenrich_heatmap <- function(enrich_res, top_n = 5, x_label = "Cell Type",
                             plt_title = "GO Enrichment Heatmap", subtitle = "",
                             xcol = rep("grey40", length(enrich_res))){
  # Select top N GO terms per module
  filt <- sapply(enrich_res,
                 function(x) !is.null(x) && nrow(x) > 0)
  filtered_results <- enrich_res[filt]
  xcol <- xcol[filt]
  
  plot_data <- lapply(names(filtered_results), function(mod) {
    df <- filtered_results[[mod]]@result
    df$Module <- mod
    dplyr::arrange(df, p.adjust) %>%
      dplyr::slice_head(n = top_n)
  }) %>% dplyr::bind_rows()
  
  if (nrow(plot_data) == 0){
    return(ggplot(data.frame()))
  }
  # Truncate GO term names
  plot_data$Term <- stringr::str_trunc(plot_data$Description, 40)
  plot_data <-  plot_data[!duplicated(plot_data[c("Module","Term")]),]
  
  # Create module x GO term matrix for clustering modules
  mat_df <- plot_data %>%
    dplyr::mutate(logp = -log10(p.adjust)) %>%
    dplyr::select(Module, Term, logp) %>%
    tidyr::pivot_wider(names_from = Term, values_from = logp, values_fill = 0)
  
  mat <- as.matrix(mat_df[, -1])
  rownames(mat) <- mat_df$Module
  
  # Cluster modules
  row_order <- tryCatch(hclust(dist(mat))$order, error=function(e){seq(along = mat)})
  ordered_modules <- rownames(mat)[row_order]
  xcol <-  xcol[row_order]
  plot_data$Module <- factor(plot_data$Module, levels = ordered_modules)
  # Plot heatmap with GO terms on y-axis
  ggplot2::ggplot(plot_data, ggplot2::aes(x = Module, y = Term)) +
    ggplot2::geom_tile(ggplot2::aes(fill = -log10(p.adjust)), color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Count), size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "red", name = "-log10 adj p") +
    ggplot2::labs(x = x_label, y = "GO Term", title = plt_title, subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = xcol),
                   axis.text.y = ggplot2::element_text(size = 10),
                   panel.grid = ggplot2::element_blank())
}

if (btwn){
  titlecmp <- sub("\\.h5ad_", " vs. ", basename(de_dir))
  titlecmp <- gsub("\\.[^.]+$", "", titlecmp)
  
} else{
  titlecmp <- gsub("\\.[^.]+$", "", basename(de_dir))
}
if (updown == "up"){
  risk <- "increase burden"
} else {
  risk <- "decrease burden"
}
for(feat in names(mdd_enrich)){
  for(sc in names(mdd_enrich[[1]][[1]])){
    enrich_lis <- lapply(mdd_enrich[[feat]], function(x){
      x[[sc]]
    })
    if(all(sapply(enrich_lis, is.null))){
      next
    }
    p <- goenrich_heatmap(enrich_lis, 
                          plt_title = "GO:BP Enrichment Heatmap", 
                          subtitle = paste(titlecmp, feat, sc, risk, sep = ", "))
    
    print(p)
    ggsave(paste0(de_dir,"/enrich_heatmap_", feat, "_", sc, "_", updown, ".png"), p, width = 10, height = 13)
    print(sc)
  }
}
