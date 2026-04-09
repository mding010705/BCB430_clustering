# overlap heatmap between the MDD gene lists and my age DEs

library(ggplot2)
library(tidyr)


cell_type_contr_rho_heatmap <- function(stat = "ding_overlap_size",
                                        pval = "ding_phyper", variable = "sex",
                                        contr_dir = "C:/Users/megan_ding/Desktop/BCB430/cluster"){
  all_dir <- list.files(path = contr_dir, full.names = TRUE, recursive = FALSE)
  all_dir <-  all_dir[!grepl("nonPFC", all_dir)]
  all_dir_name <- basename(all_dir)
  all_dir_name <- sub("\\.h5ad_", " vs ", all_dir_name)
  all_dir_name <- sub("\\.h5ad", "", all_dir_name)
  names(all_dir_name) <- all_dir
  df <- NULL
  for(d in all_dir){
    if (is.null(df)){
      df <- readRDS(file = paste0(d, "/mdd_overlaps.rds"))
      if(nrow(df) > 0){
        df$type <- all_dir_name[d]
      }
      
    } else {
      ndf <- readRDS(file = paste0(d, "/mdd_overlaps.rds"))
      if(nrow(ndf) > 0){
        ndf$type <- all_dir_name[d]
        df <- rbind(df, ndf)
      }
      
    }
  }
  if (variable == "sex"){
    row_mask <- df$contrast == variable | df$type == "PFC_female vs PFC_male" | df$type == "isocortex_female vs isocortex_male"
  } else if (variable == "age_cat"){
    row_mask <- df$contrast == variable | df$type == "PFC_aged vs PFC_adult" | df$type == "isocortex_aged vs isocortex_adult"
  } else {
    row_mask <- df$contrast == variable 
  }
  df <- df[row_mask, c("type", "cell_type", stat, pval)]
  df[[pval]] <- p.adjust(df[[pval]], method = "fdr")
  df_wide <- pivot_wider(
    df[, -4],
    names_from = cell_type,
    values_from = stat
  )
  if (grepl("rho", stat)){
    df_wide[is.na(df_wide)] <- 0
  } else {
    df_wide[is.na(df_wide)] <- -100
  }
  
  row_hclust <- hclust(dist(df_wide[, -c(1)], method = "euclidean"), method = "ward.D")
  
  # Get the order of the clustered items
  row_order <- row_hclust$order
  
  col_hclust <- hclust(dist(t(df_wide[, -c(1)]), method = "euclidean"), method = "ward.D")
  
  # Get the order of the clustered items
  col_order <- col_hclust$order
  df$cell_type <- factor(df$cell_type, levels = colnames(df_wide)[-1][col_order])
  df$type <- factor(df$type, levels = df_wide$type[row_order])
  
  df$pstar <- ifelse(df[, pval] < 0.0001, "***",
                     ifelse(df[, pval] < 0.001, "**",
                            ifelse(df[, pval] < 0.05, "*", "")))
  colnames(df)[colnames(df) == stat] <- "stat"
  if (grepl("overlap", stat)){
    leg_name <- "overlap size"
  } else {
    leg_name <- "Spearman's rho"
  }
  df <- df[!is.na(df$cell_type), ]
  if(grepl("rho", stat)){
    p <- (ggplot(df, aes(cell_type, type)) +
             geom_tile(aes(fill = stat)) +
             scale_fill_gradient2(low = "blue", high = "red", name = leg_name) +
             ggtitle(paste(stat, variable, sep = ", ")) +
             geom_text(aes(label = pstar), color = "green", size=3) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.background = element_rect(fill = "grey", colour = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()))
  } else{
    p <- (ggplot(df, aes(cell_type, type)) +
             geom_tile(aes(fill = stat)) +
             scale_fill_gradient(low = "white", high = "red", name = leg_name) +
             ggtitle(paste(stat, variable, sep = ", ")) +
             geom_text(aes(label = pstar), color = "green", size=3, vjust = -1) +
             geom_text(aes(label = stat), vjust = 0.5, size = 4, color = "black") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.background = element_rect(fill = "grey", colour = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()))
  }
  ggsave(plot = p, filename = paste0("C:/Users/megan_ding/Desktop/BCB430/images/overlap_hm_", 
                                     stat, "_", variable, ".png"), height = 8, width = 17, dpi = 600)
  return(p)
}

cell_type_contr_rho_heatmap(stat = "ding_overlap_size",
                            pval = "ding_phyper", variable = "sex")
cell_type_contr_rho_heatmap(stat = "nordic_early_overlap_size",
                            pval = "nordic_early_phyper", variable = "sex")
cell_type_contr_rho_heatmap(stat = "nordic_late_overlap_size",
                            pval = "nordic_late_phyper", variable = "sex")
cell_type_contr_rho_heatmap(stat = "nordic_mixed_overlap_size",
                            pval = "nordic_mixed_phyper", variable = "sex")

cell_type_contr_rho_heatmap(stat = "ding_overlap_size",
                            pval = "ding_phyper", variable = "age_cat")
cell_type_contr_rho_heatmap(stat = "nordic_early_overlap_size",
                            pval = "nordic_early_phyper", variable = "age_cat")
cell_type_contr_rho_heatmap(stat = "nordic_late_overlap_size",
                            pval = "nordic_late_phyper", variable = "age_cat")
cell_type_contr_rho_heatmap(stat = "nordic_mixed_overlap_size",
                            pval = "nordic_mixed_phyper", variable = "age_cat")

cell_type_contr_rho_heatmap(stat = "ding_overlap_size",
                            pval = "ding_phyper", variable = "roi")
cell_type_contr_rho_heatmap(stat = "nordic_early_overlap_size",
                            pval = "nordic_early_phyper", variable = "roi")
cell_type_contr_rho_heatmap(stat = "nordic_late_overlap_size",
                            pval = "nordic_late_phyper", variable = "roi")
cell_type_contr_rho_heatmap(stat = "nordic_mixed_overlap_size",
                            pval = "nordic_mixed_phyper", variable = "roi")


cell_type_contr_rho_heatmap(stat = "ding_mixed_rho",
                            pval = "ding_mixed_pval", variable = "sex")
cell_type_contr_rho_heatmap(stat = "nordic_early_rho",
                            pval = "nordic_early_pval", variable = "sex")
cell_type_contr_rho_heatmap(stat = "nordic_late_rho",
                            pval = "nordic_late_pval", variable = "sex")
cell_type_contr_rho_heatmap(stat = "nordic_mixed_rho",
                            pval = "nordic_mixed_pval", variable = "sex")

cell_type_contr_rho_heatmap(stat = "ding_mixed_rho",
                            pval = "ding_mixed_pval", variable = "age_cat")
cell_type_contr_rho_heatmap(stat = "nordic_early_rho",
                            pval = "nordic_early_pval", variable = "age_cat")
cell_type_contr_rho_heatmap(stat = "nordic_late_rho",
                            pval = "nordic_late_pval", variable = "age_cat")
cell_type_contr_rho_heatmap(stat = "nordic_mixed_rho",
                            pval = "nordic_mixed_pval", variable = "age_cat")

cell_type_contr_rho_heatmap(stat = "ding_mixed_rho",
                            pval = "ding_mixed_pval", variable = "roi")
cell_type_contr_rho_heatmap(stat = "nordic_early_rho",
                            pval = "nordic_early_pval", variable = "roi")
cell_type_contr_rho_heatmap(stat = "nordic_late_rho",
                            pval = "nordic_late_pval", variable = "roi")
cell_type_contr_rho_heatmap(stat = "nordic_mixed_rho",
                            pval = "nordic_mixed_pval", variable = "roi")

