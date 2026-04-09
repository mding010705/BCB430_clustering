# heatmap of top GO:BP term enrichments for our fixed cells of interest

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

de_dir <- commandArgs(trailingOnly = TRUE)[1]
updown <- commandArgs(trailingOnly = TRUE)[2]

# if DEA was performed between different subsets or within one
matches <- gregexpr("h5ad", de_dir, fixed = TRUE)
counts <- sum(lengths(matches))
btwn <- counts == 2

mdd_enrich <- readRDS(file = paste0(de_dir, "/", "mdd_enrich_", updown, ".rds"))


goenrich_heatmap <- function(enrich_res, top_n = 5, x_label = "Cell Type",
                             plt_title = "GO Enrichment Heatmap", subtitle = "",
                             sel_cells = c("Astro-TE NN", "OPC NN", "Oligo NN",
                                           "COP NN_1", "Microglia NN", "Peri NN"), 
                             xcol = rep("grey40", length(enrich_res))){
  # Select top N GO terms per module
  if (!is.null(sel_cells)){
    enrich_res <- enrich_res[sel_cells]
  }
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
    ggsave(paste0(de_dir,"/enrich_heatmap_", feat, "_", sc, "_", updown, "_sel_cells.png"), p, width = 9, height = 10,
           dpi = 900)
    print(sc)
  }
}
