# between subset clusters DEA
# based on code from Jiaqi

# DE between clusters of same name
# iso vs. pfc
# iso male vs. iso female, pfc male vs. pfc female
# iso adult vs. iso aged, pfc adult vs. pfc aged
# iso male adult vs. iso female adult, pfc male adult vs. pfc female adult
# iso male aged vs. iso female aged, pfc male aged vs. pfc female aged
# iso male adult vs. iso female aged, iso female adult vs. iso male aged
# pfc male adult vs. pfc female aged, pfc female adult vs. pfc male aged

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

expr1 <- commandArgs(trailingOnly = TRUE)[1]
expr2 <- commandArgs(trailingOnly = TRUE)[2]

output_dir <- "C:/Users/megan_ding/Desktop/BCB430/cluster/"
expr1 <- paste0("C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/", expr1)
expr2 <- paste0("C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/", expr2)

# for each data strat

SeuratObj1 <- anndataR::read_h5ad(expr1, as = "Seurat")

# Ensure default assay is RNA
DefaultAssay(SeuratObj1) <- "RNA"

# Assign layer X to both layer data and layer counts, since layer x is already normalized
LayerData(SeuratObj1, assay = "RNA", layer = "data") <- as.matrix(LayerData(SeuratObj1, assay = "RNA", layer = "X"))
LayerData(SeuratObj1, assay = "RNA", layer = "counts") <- as.matrix(LayerData(SeuratObj1, assay = "RNA", layer = "X"))
SeuratObj1@meta.data$region <- sub("[^a-zA-Z].*", "", basename(expr1))
SeuratObj2 <- anndataR::read_h5ad(expr2, as = "Seurat")

# Ensure default assay is RNA
DefaultAssay(SeuratObj2) <- "RNA"

# Assign layer X to both layer data and layer counts, since layer x is already normalized
LayerData(SeuratObj2, assay = "RNA", layer = "data") <- as.matrix(LayerData(SeuratObj2, assay = "RNA", layer = "X"))
LayerData(SeuratObj2, assay = "RNA", layer = "counts") <- as.matrix(LayerData(SeuratObj2, assay = "RNA", layer = "X"))
SeuratObj2@meta.data$region <- sub("[^a-zA-Z].*", "", basename(expr2))


all_seurat <-  merge(SeuratObj1, SeuratObj2, 
                    project = "All", merge.data = TRUE)
all_seurat <- JoinLayers(all_seurat)

all_seurat@meta.data$reg_sex_age <- paste0(all_seurat@meta.data$region, "_", all_seurat@meta.data$sex, "_", all_seurat@meta.data$age_cat)
if(length(unique(all_seurat$region)) > 1){
  voi <- "region"
  id1 <- grep("PFC", all_seurat$region, value = T)[1]
} else if(length(unique(all_seurat$sex)) > 1 & length(unique(all_seurat$age_cat)) > 1){
  voi <- "reg_sex_age"
  id1 <- grep("F_", all_seurat$reg_sex_age, value = T)[1]
} else if (length(unique(all_seurat$sex)) > 1) {
  voi <- "sex"
  id1 <- grep("F_", all_seurat$sex, value = T)[1]
} else {
  voi <- "age_cat"
  id1 <- grep("aged", all_seurat$age_cat, value = T)[1]
}


merge.result1 <- readRDS(paste0("C:/Users/megan_ding/Desktop/BCB430/clust_res_hicat/labeled/", sub(".h5ad", "", basename(expr1)), ".rds"))
names(merge.result1$cl) <- paste0(names(merge.result1$cl), "_1")
merge.result2 <- readRDS(paste0("C:/Users/megan_ding/Desktop/BCB430/clust_res_hicat/labeled/", sub(".h5ad", "", basename(expr2)), ".rds"))
names(merge.result2$cl) <- paste0(names(merge.result2$cl), "_2")

hicat_clusters <- c(merge.result1$cl, merge.result2$cl)
all_seurat <- AddMetaData(object = all_seurat,
                                  metadata = hicat_clusters,
                                  col.name = "hicat_cluster")

Idents(all_seurat) <- "hicat_cluster"
all_seurat <- FindVariableFeatures(all_seurat) %>%
  ScaleData() %>%
  RunPCA() 

target_subclasses <- unique(all_seurat$hicat_cluster)
variables_to_test <- c(voi)

# 2. Initialize a list to store results
# Structure: de_results$Subclass_Name$Variable_Name
de_results <- list()
num_genes_expr <- c()

# 3. Main Loop
for (subclass in target_subclasses) {
  
  message(paste0("\nProcessing Subclass: ", subclass, " ------------------"))
  
  # A. Subset the object to just this subclass
  # We use subset() to isolate the cells
  sub_obj <- subset(all_seurat, subset = hicat_cluster == subclass)
  counts_matrix <- GetAssayData(object = sub_obj, assay = "RNA", layer = "counts")
  num_cells_in_cluster <- ncol(counts_matrix)
  min_cells_threshold <- floor(0.01 * num_cells_in_cluster)
  genes_expressed_above_threshold <- Matrix::rowSums(counts_matrix > 0) >= min_cells_threshold
  num_genes_expr[[subclass]] <- sum(genes_expressed_above_threshold)
  # Initialize list for this subclass
  de_results[[subclass]] <- list()
  
  # B. Loop through the variables (Sex, Age, ROI)
  for (var in variables_to_test) {
    
    message(paste("Testing variable:", var))
    
    # Check 1: Does the variable exist in metadata?
    if (!var %in% colnames(sub_obj@meta.data)) {
      message(paste("Skipping: Variable", var, "not found."))
      next
    }
    
    # Set the identity to the variable of interest (e.g., set Idents to "sex")
    Idents(sub_obj) <- var
    # Check 2: Are there at least 2 groups?
    # (e.g., If a subclass only exists in Males, we can't do Male vs Female)
    unique_groups <- unique(Idents(sub_obj))
    unique_groups <- unique_groups[!is.na(unique_groups)] # Remove NAs
    
    if (length(unique_groups) < 2) {
      message(paste("Skipping: Only 1 group found for", var, "(", unique_groups, ")"))
      next
    }
    
    # C. Run DE Analysis
    # We use FindAllMarkers to handle both binary (M vs F) and multi-class (ROI A vs B vs C) cases automatically.
    # It returns markers for each group compared to the others in this subclass.
    tryCatch({
      markers <- FindMarkers(
        object = sub_obj,
        ident.1 = id1,
        only.pos = FALSE, 
        verbose = FALSE,
        logfc.threshold = 0.1,
        min.pct = 0.01,
        min.diff.pct = -Inf
      )
      
      # D. Add your Custom Metrics (p_FC, etc.)
      if (nrow(markers) > 0) {
        # Calculate p_FC
        markers$p_FC <- (1 - markers$p_val) * markers$avg_log2FC
        
        # Calculate difference in proportion (pct.1 - pct.2)
        # Note: FindAllMarkers returns pct.1 and pct.2 automatically
        markers$adj_diff <- abs(markers$pct.1 - markers$pct.2) / pmax(markers$pct.1, markers$pct.2)
        
        # Calculate combined score
        markers$diff_pfc <- markers$adj_diff * markers$p_FC
        
        # Save to list
        de_results[[subclass]][[var]] <- markers
        message(paste("Found", nrow(markers), "markers."))
        
      } else {
        message("No significant markers found.")
      }
      
    }, error = function(e) {
      message(paste("Error running DE for", var, ":", e$message))
    })
  }
}

dir.create(paste0(output_dir, basename(expr1), "_", basename(expr2)), showWarnings = FALSE, recursive = TRUE)

saveRDS(de_results, file = paste0(output_dir, basename(expr1), "_", basename(expr2), "/DE_results.rds"))

summary_list <- list()

# Loop through the nested list
for (subclass in names(de_results)) {
  for (variable in names(de_results[[subclass]])) {
    
    # Get the marker dataframe
    markers <- de_results[[subclass]][[variable]]
    markers <- markers[markers$p_val_adj < 0.05,]
    # Count rows (if NULL, count is 0)
    num_genes_up <- if (is.null(markers)) 0 else sum(markers$avg_log2FC > 0)
    num_genes_down <- if (is.null(markers)) 0 else sum(markers$avg_log2FC < 0)
    
    # Store in the list
    summary_list[[length(summary_list) + 1]] <- data.frame(
      Subclass = subclass,
      Variable = variable,
      Count_up = num_genes_up,
      Count_down = num_genes_down,
      Count_up_adj = num_genes_up/num_genes_expr[[subclass]],
      Count_down_adj = num_genes_down/num_genes_expr[[subclass]]
    )
  }
}

# 2. Combine into a clean Data Frame
plot_data <- do.call(rbind, summary_list)

# Optional: Clean up variable names for the plot (Capitalize)
plot_data$Variable <- as.factor(plot_data$Variable)

# 3. Create the Plot
text <- paste0(basename(gsub("\\..*$", "", expr1)), " vs. ", basename(gsub("\\..*$", "", expr2)), ": # Sig. DE Genes per Hicat Cluster")


axis_margin <- 5.5

p_up <- ggplot(plot_data, aes(x = Subclass, y = Count_up, fill = Variable)) +
  geom_col(position = "dodge") + # 'dodge' puts bars side-by-side
  coord_flip() + # Flip coordinates to make subclass names readable
  labs(
    title = "",
    y = "Number of DE Upreg. Genes",
    fill = "Criteria"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 7), # Adjust text size if you have many subclasses
    legend.position = "none",
    plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
    axis.text.y.left = element_text(margin = margin(0, axis_margin, 0, axis_margin))
  ) +
  geom_text(aes(label = Count_up), position = position_dodge(width = 0.9), hjust = -0.2, size = 3) # Add numbers

p_down <- ggplot(plot_data, aes(x = Subclass, y = Count_down, fill = Variable)) +
  geom_col(position = "dodge") + # 'dodge' puts bars side-by-side
  coord_flip() + # Flip coordinates to make subclass names readable
  scale_y_reverse() +
  labs(
    title = text,
    y = "Number of DE Downreg. Genes",
    fill = "Criteria"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(), # Adjust text size if you have many subclasses
    legend.position = "left",
    plot.margin = margin(axis_margin, 0, axis_margin, axis_margin)
  ) +
  geom_text(aes(label = Count_down), position = position_dodge(width = 0.9), hjust = 1, size = 3) # Add numbers


p <- ggarrange(p_down, p_up,
               ncol = 2,nrow = 1,heights = c(0.5,5))
# 4. Display and Save
print(p)
ggsave(paste0(output_dir, basename(expr1), "_", basename(expr2),"/DE_Counts_Summary_Hicat.png"), p, width = 17, height = 8)

text <- paste0(basename(gsub("\\..*$", "", expr1)), " vs. ", basename(gsub("\\..*$", "", expr2)), ": Prop. Sig. DE Genes per Hicat Cluster")


axis_margin <- 5.5

p_up <- ggplot(plot_data, aes(x = Subclass, y = Count_up_adj, fill = Variable)) +
  geom_col(position = "dodge") + # 'dodge' puts bars side-by-side
  coord_flip() + # Flip coordinates to make subclass names readable
  labs(
    title = "",
    y = "Prop. DE Upreg. Genes",
    fill = "Criteria"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 7), # Adjust text size if you have many subclasses
    legend.position = "none",
    plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
    axis.text.y.left = element_text(margin = margin(0, axis_margin, 0, axis_margin))
  ) +
  geom_text(aes(label = signif(Count_up_adj, digits = 3)), position = position_dodge(width = 0.9), hjust = -0.2, size = 3) # Add numbers

p_down <- ggplot(plot_data, aes(x = Subclass, y = Count_down_adj, fill = Variable)) +
  geom_col(position = "dodge") + # 'dodge' puts bars side-by-side
  coord_flip() + # Flip coordinates to make subclass names readable
  scale_y_reverse() +
  labs(
    title = text,
    y = "Prop. DE Downreg. Genes",
    fill = "Criteria"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(), # Adjust text size if you have many subclasses
    legend.position = "left",
    plot.margin = margin(axis_margin, 0, axis_margin, axis_margin)
  ) +
  geom_text(aes(label = signif(Count_down_adj, digits = 3)), position = position_dodge(width = 0.9), hjust = 1, size = 3) # Add numbers


p <- ggarrange(p_down, p_up,
               ncol = 2,nrow = 1,heights = c(0.5,5))
# 4. Display and Save
print(p)
ggsave(paste0(output_dir, basename(expr1), "_", basename(expr2),"/DE_Prop_Summary_Hicat.png"), p, width = 17, height = 8)
gc()

