# intersect DE genes with MDD gene lists, evaluate with hypergeometric tests, RRHO tests,
# and do GO:BP enrichment on the genes in this overlap

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
                           ding_phyper = NA, ding_overlap_size = NA,
                           ding_mixed_rho = NA,
                           ding_mixed_pval = NA,
                           ding_female_rho = NA,
                           ding_female_pval = NA,
                           ding_male_rho = NA,
                           ding_male_pval = NA,
                           nordic_mixed_phyper = NA, 
                           nordic_mixed_overlap_size = NA,
                           nordic_mixed_rho = NA,
                           nordic_mixed_pval = NA,
                           nordic_early_phyper = NA, 
                           nordic_early_overlap_size = NA,
                           nordic_early_rho = NA,
                           nordic_early_pval = NA,
                           nordic_late_phyper = NA, 
                           nordic_late_overlap_size = NA,
                           nordic_late_rho = NA,
                           nordic_late_pval = NA)
for(ct in names(btwn_clust_de)){
  for (feat in names(btwn_clust_de[[ct]])){
    des <- btwn_clust_de[[ct]][[feat]]
    des_sig <- des[des$p_val_adj < 0.05, ] 
    ding_overlap <- rownames(des_sig)[rownames(des_sig) %in% unique(ding_mouse$mouse_gene)]
    degs_human <- mouse2human[which(mouse2human$mouse_gene %in% rownames(des_sig)), ]
    # number genes in expression matrix 32285
    ding_hyper <- phyper(q = sum(!duplicated(ding_overlap)) - 1,
                          m = sum(!duplicated(ding_mouse$mouse_gene)),
                          n = sum(!duplicated(all_genes_human$mouse_gene)) - 
                                    sum(!duplicated(ding_mouse$mouse_gene)),
                          k = sum(!duplicated(degs_human$mouse_gene)),
                          lower.tail = FALSE)
    ding_overlap_size <- sum(!duplicated(ding_overlap))
    
    DE_df <- data.frame(gene = rownames(des), pval = -des$p_val, logFC = des$avg_log2FC)
    DE_df$rank <- -log10(-DE_df$pval) * sign(DE_df$logFC)
    ding_df <- ding_mouse[, c("mouse_gene", "mixed_p", "female_p", "male_p", "Average Effect size")]
    ding_df$rank <- -log10(ding_df$mixed_p) * sign(ding_df$`Average Effect size`)
    ding_df <- ding_df[!duplicated(ding_df$mouse_gene), ]
    ding_df$female_p <- -ding_df$female_p
    ding_df$male_p <- -ding_df$male_p
    DE_df <- DE_df[DE_df$gene %in% ding_df$mouse_gene, ]
    ding_df <- ding_df[ding_df$mouse_gene %in% DE_df$gene, ]
    ss <- ceiling(nrow(DE_df)/70)
    ding_rrho <- tryCatch(RRHO(DE_df[, c("gene", "rank")], ding_df[, c("mouse_gene", "rank")], labels = c(paste0(ct, " ", feat), "Ding 2015 mixed sex"),
                      alternative = "two.sided", plots = TRUE, outputdir = de_dir,
                      BY = FALSE, log10.ind = TRUE, stepsize = ss), 
                      error=function(e){0})
    ding_rrho <-tryCatch(RRHO(DE_df[, c("gene", "pval")], ding_df[, c("mouse_gene", "female_p")], labels = c(paste0(ct, " ", feat), "Ding 2015 female"),
                      alternative = "two.sided", plots = TRUE, outputdir = de_dir,
                      BY = FALSE, log10.ind = TRUE, stepsize = ss), 
                      error=function(e){0})
    ding_rrho <- tryCatch(RRHO(DE_df[, c("gene", "pval")], ding_df[, c("mouse_gene", "male_p")], labels = c(paste0(ct, " ", feat), "Ding 2015 male"),
                      alternative = "two.sided", plots = TRUE, outputdir = de_dir,
                      BY = FALSE, log10.ind = TRUE, stepsize = ss), 
                      error=function(e){0})
    cor_test_ding_mixed <- tryCatch(
      cor.test(DE_df$rank, ding_df$rank, method = "spearman"), 
      error=function(e){a <- list()
      a[["estimate"]][["rho"]] <- NA
      a[["p.value"]] <- NA
      return(a)})

    cor_test_ding_female <- tryCatch(cor.test(DE_df$pval, ding_df$female_p, method = "spearman"), 
                                     error=function(e){a <- list()
                                     a[["estimate"]][["rho"]] <- NA
                                     a[["p.value"]] <- NA
                                     return(a)})
    
    cor_test_ding_male <- tryCatch(cor.test(DE_df$pval, ding_df$male_p, method = "spearman"), 
                                   error=function(e){a <- list()
                                   a[["estimate"]][["rho"]] <- NA
                                   a[["p.value"]] <- NA
                                   return(a)})
    
    
    nordic_overlap_mixed <- rownames(des_sig)[rownames(des_sig) %in% unique(nordic_mouse$mouse_gene)]
    # number genes in expression matrix 32285
    nordic_hyper_mixed <- phyper(q = sum(!duplicated(nordic_overlap_mixed)) - 1,
                         m = sum(!duplicated(nordic_mouse$mouse_gene)),
                         n = sum(!duplicated(all_genes_human$mouse_gene)) - 
                           sum(!duplicated(nordic_mouse$mouse_gene)),
                         k = sum(!duplicated(degs_human$mouse_gene)),
                         lower.tail = FALSE)
    nordic_overlap_size_mixed <- sum(!duplicated(nordic_overlap_mixed))
    
    nordic_overlap_early <- rownames(des_sig)[rownames(des_sig) %in% unique(nordic_eo_mouse$mouse_gene)]
    # number genes in expression matrix 32285
    nordic_hyper_early <- phyper(q = sum(!duplicated(nordic_overlap_early)) - 1,
                           m = sum(!duplicated(nordic_eo_mouse$mouse_gene)),
                           n = sum(!duplicated(all_genes_human$mouse_gene)) - 
                             sum(!duplicated(nordic_eo_mouse$mouse_gene)),
                           k = sum(!duplicated(degs_human$mouse_gene)),
                           lower.tail = FALSE)
    nordic_overlap_size_early <- sum(!duplicated(nordic_overlap_early))
    
    nordic_overlap_late <- rownames(des_sig)[rownames(des_sig) %in% unique(nordic_lo_mouse$mouse_gene)]
    # number genes in expression matrix 32285
    nordic_hyper_late <- phyper(q = sum(!duplicated(nordic_overlap_late)) - 1,
                           m = sum(!duplicated(nordic_lo_mouse$mouse_gene)),
                           n = sum(!duplicated(all_genes_human$mouse_gene)) - 
                             sum(!duplicated(nordic_lo_mouse$mouse_gene)),
                           k = sum(!duplicated(degs_human$mouse_gene)),
                           lower.tail = FALSE)
    nordic_overlap_size_late <- sum(!duplicated(nordic_overlap_late))
    
    DE_df <- data.frame(gene = rownames(des), pval = des$p_val, logFC = des$avg_log2FC)
    DE_df$rank <- -log10(DE_df$pval) * sign(DE_df$logFC)
    nordic_df <- nordic_mouse[, c("mouse_gene", "ZSTAT", "p.fdr")]
    nordic_df$rank <- -log10(nordic_df$p.fdr) * sign(nordic_df$ZSTAT)
    nordic_df <- nordic_df[!duplicated(nordic_df$mouse_gene), ]
    
    DE_df <- DE_df[DE_df$gene %in% nordic_df$mouse_gene, ]
    nordic_df <- nordic_df[nordic_df$mouse_gene %in% DE_df$gene, ]
    ss <- ceiling(nrow(DE_df)/70)
    
    nordic_rrho <- tryCatch(RRHO(DE_df[, c("gene", "rank")], nordic_df[, c("mouse_gene", "rank")], labels = c(paste0(ct, " ", feat), "Shorter 2025 mixed onset"),
                      alternative = "two.sided", plots = TRUE, outputdir = de_dir,
                      BY = FALSE, log10.ind = TRUE, stepsize = ss), 
                      error=function(e){0})
    
    
    cor_test_nordic_mixed <- tryCatch(cor.test(DE_df$rank, nordic_df$rank, method = "spearman"), 
                                      error=function(e){a <- list()
                                      a[["estimate"]][["rho"]] <- NA
                                      a[["p.value"]] <- NA
                                      return(a)})
    
    
    DE_df <- data.frame(gene = rownames(des), pval = des$p_val, logFC = des$avg_log2FC)
    DE_df$rank <- -log10(DE_df$pval) * sign(DE_df$logFC)
    nordic_df <- nordic_eo_mouse[, c("mouse_gene", "ZSTAT", "p.fdr")]
    nordic_df$rank <- -log10(nordic_df$p.fdr) * sign(nordic_df$ZSTAT)
    nordic_df <- nordic_df[!duplicated(nordic_df$mouse_gene), ]
    
    DE_df <- DE_df[DE_df$gene %in% nordic_df$mouse_gene, ]
    nordic_df <- nordic_df[nordic_df$mouse_gene %in% DE_df$gene, ]
    nordic_rrho <- tryCatch(RRHO(DE_df[, c("gene", "rank")], nordic_df[, c("mouse_gene", "rank")], labels = c(paste0(ct, " ", feat), "Shorter 2025 early onset"),
                        alternative = "two.sided", plots = TRUE, outputdir = de_dir,
                        BY = FALSE, log10.ind = TRUE, stepsize = ss), 
                        error=function(e){0})
    
    cor_test_nordic_early <- tryCatch(cor.test(DE_df$rank, nordic_df$rank, method = "spearman"), 
                                      error=function(e){a <- list()
                                      a[["estimate"]][["rho"]] <- NA
                                      a[["p.value"]] <- NA
                                      return(a)})
    
    DE_df <- data.frame(gene = rownames(des), pval = des$p_val, logFC = des$avg_log2FC)
    DE_df$rank <- -log10(DE_df$pval) * sign(DE_df$logFC)
    nordic_df <- nordic_lo_mouse[, c("mouse_gene", "ZSTAT", "p.fdr")]
    nordic_df$rank <- -log10(nordic_df$p.fdr) * sign(nordic_df$ZSTAT)
    nordic_df <- nordic_df[!duplicated(nordic_df$mouse_gene), ]
    
    DE_df <- DE_df[DE_df$gene %in% nordic_df$mouse_gene, ]
    nordic_df <- nordic_df[nordic_df$mouse_gene %in% DE_df$gene, ]
    nordic_rrho <- tryCatch(RRHO(DE_df[, c("gene", "rank")], nordic_df[, c("mouse_gene", "rank")], labels = c(paste0(ct, " ", feat), "Shorter 2025 late onset"),
                        alternative = "two.sided", plots = TRUE, outputdir = de_dir,
                        BY = FALSE, log10.ind = TRUE, stepsize = ss), 
                        error=function(e){0})
    cor_test_nordic_late <- tryCatch(cor.test(DE_df$pval, nordic_df$rank, method = "spearman"), 
                                     error=function(e){a <- list()
                                      a[["estimate"]][["rho"]] <- NA
                                      a[["p.value"]] <- NA
                                      return(a)})
    
    
    
    
    new_df <- data.frame(cell_type = ct, contrast = feat, between = btwn, 
                         ding_phyper = ding_hyper, ding_overlap_size = ding_overlap_size,
                         ding_mixed_rho = cor_test_ding_mixed[["estimate"]][["rho"]],
                         ding_mixed_pval = cor_test_ding_mixed[["p.value"]],
                         ding_female_rho = cor_test_ding_female[["estimate"]][["rho"]],
                         ding_female_pval = cor_test_ding_female[["p.value"]],
                         ding_male_rho = cor_test_ding_male[["estimate"]][["rho"]],
                         ding_male_pval = cor_test_ding_male[["p.value"]],
                         nordic_mixed_phyper = nordic_hyper_mixed, 
                         nordic_mixed_overlap_size = nordic_overlap_size_mixed,
                         nordic_mixed_rho = cor_test_nordic_mixed[["estimate"]][["rho"]],
                         nordic_mixed_pval = cor_test_nordic_mixed[["p.value"]],
                         nordic_early_phyper = nordic_hyper_early, 
                         nordic_early_overlap_size = nordic_overlap_size_early,
                         nordic_early_rho = cor_test_nordic_early[["estimate"]][["rho"]],
                         nordic_early_pval = cor_test_nordic_early[["p.value"]],
                         nordic_late_phyper = nordic_hyper_late, 
                         nordic_late_overlap_size = nordic_overlap_size_late,
                         nordic_late_rho = cor_test_nordic_late[["estimate"]][["rho"]],
                         nordic_late_pval = cor_test_nordic_late[["p.value"]])
    
    
    gene_enrichment_ding <- gene_enrichment(gene_set = ding_overlap, 
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
    mdd_overlaps <- rbind(mdd_overlaps, new_df)
  }
}
mdd_overlaps <- mdd_overlaps[-1, ]

saveRDS(mdd_enrich, file = paste0(de_dir, "/", "mdd_enrich.rds"))
saveRDS(mdd_overlaps, file = paste0(de_dir, "/", "mdd_overlaps.rds"))

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
for(feat in names(mdd_enrich)){
  for(sc in names(mdd_enrich[[1]][[1]])){
    enrich_lis <- lapply(mdd_enrich[[feat]], function(x){
      x[[sc]]
    })
    p <- goenrich_heatmap(enrich_lis, 
                          plt_title = "GO:BP Enrichment Heatmap", 
                          subtitle = paste(titlecmp, feat, sc, sep = ", "))
    
    print(p)
    ggsave(paste0(de_dir,"/enrich_heatmap_", feat, "_", sc, ".png"), p, width = 10, height = 13)
    
  }
}
