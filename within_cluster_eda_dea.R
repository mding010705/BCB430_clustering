# EDA and DEA within each subset and cell type
# thank you Jiaqi for the starter code

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

iso_expression_h5_dir <- commandArgs(trailingOnly = TRUE)[1]
output_dir <- "C:/Users/megan_ding/Desktop/BCB430/cluster/"
hitcat_clust_file <- paste0("C:/Users/megan_ding/Desktop/BCB430/clust_res_hicat/labeled/", sub(".h5ad", "", basename(iso_expression_h5_dir)), ".rds")

SeuratObj <- anndataR::read_h5ad(iso_expression_h5_dir, as = "Seurat")

# Ensure default assay is RNA
DefaultAssay(SeuratObj) <- "RNA"

# Assign layer X to both layer data and layer counts, since layer x is already normalized
LayerData(SeuratObj, assay = "RNA", layer = "data") <- LayerData(SeuratObj, assay = "RNA", layer = "X")
LayerData(SeuratObj, assay = "RNA", layer = "counts") <- LayerData(SeuratObj, assay = "RNA", layer = "X")

`iso_seurat` <- SeuratObj
# expression matrix
`iso_expression_matrix` <- LayerData(SeuratObj, assay = "RNA", layer = "data")

############
# EDA
# 1. Define the variables you want to plot
vars_to_plot <- c("subclass_label", "roi", "sex", "age_cat")

# 2. Initialize a list to store the plots
pie_plots <- list()

# 3. Loop through each variable to create a pie chart
for (var in vars_to_plot) {
  
  # A. Summarize the data (Count and Percent)
  # We use .data[[var]] to dynamically select the column inside the loop
  plot_data <- `iso_seurat`@meta.data %>%
    group_by(.data[[var]]) %>%
    summarise(count = n()) %>%
    mutate(
      prop = count / sum(count),
      percentage = round(prop * 100, 1),
      # Create a label that looks like: "150 (10.5%)"
      label_text = paste0(count, "\n(", percentage, "%)")
    ) %>%
    arrange(desc(.data[[var]])) # Fix order for plotting
  
  # B. Create the Pie Chart
  # Pie charts in ggplot are just Bar Charts with polar coordinates
  p <- ggplot(plot_data, aes(x = "", y = prop, fill = .data[[var]])) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    theme_void() + # Removes axes and grey background
    labs(
      title = paste("Distribution by", var),
      fill = var
    ) +
    # Add text labels in the middle of each slice
    geom_text(aes(label = label_text),
              position = position_stack(vjust = 0.5),
              size = 3.5, fontface = "bold") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
  
  # Save to list
  pie_plots[[var]] <- p
}

# 4. Combine the 4 plots into one image
# (Requires 'patchwork' library, or use cowplot::plot_grid / gridExtra::grid.arrange)
final_plot <- wrap_plots(pie_plots, ncol = 2)

ggsave(paste0(output_dir, basename(iso_expression_h5_dir),"/Cell_Demographics_PieCharts.png"),
       final_plot, width = 12, height = 10, bg = "white")
# Display in RStudio
print(final_plot)

merge.result <- readRDS(hicat_clust_file)

hicat_clusters <- merge.result$cl
`iso_seurat` <- AddMetaData(object = `iso_seurat`,
                                  metadata = hicat_clusters,
                                  col.name = "hicat_cluster")
# Set the identity to be hicat cluster result
Idents(`iso_seurat`) <- "hicat_cluster"

# Calculate Umap
`iso_seurat` <- FindVariableFeatures(`iso_seurat`) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30, reduction.name = "umap_dim30")

# Visualize Umap
# 1. Define the list of variables to plot
group_vars <- c("hicat_cluster", "bigcat_cluster", "class_label", "subclass_label",
                "supertype_label", "cluster_label", "roi", "sex", "age_cat")

# 2. Loop through each variable
for (var in group_vars) {
  if (basename(iso_expression_h5_dir) == "isocortex.h5ad" && var == "bigcat_cluster"){
    next
  }
  # Print progress to console
  message(paste("Plotting and saving:", var))
  
  # Generate the plot
  p <- DimPlot(`iso_seurat`, group.by = var, label = TRUE, repel = TRUE, reduction = "umap_dim30") +
    ggtitle(paste(var, "on Seurat UMAP")) + NoLegend()
  
  # 3. Save the plot
  plot(p)
  filename <- paste0(output_dir, basename(iso_expression_h5_dir),"/UMAP_dim30_", var, ".png")
  
  ggsave(filename = filename, plot = p, width = 10, height = 7.5, dpi = 320)
}

`iso_marker` <- FindAllMarkers(object = `iso_seurat`)
`iso_marker` <- split(`iso_marker`, `iso_marker`$cluster)

variables_to_test <- c("sex", "age_cat", "roi")

# 2. Initialize a list to store results
# New Structure: de_results$Variable_Name (No longer nested by subclass)
global_de_results <- list()

# 3. Main Loop (Iterate through variables only)
for (var in variables_to_test) {
  
  message(paste0("\nProcessing Variable Globally: ", var, " ------------------"))
  
  # Use the full object directly
  # (No subsetting happens here)
  object_to_test <- `iso_seurat`
  
  # Check 1: Does the variable exist in metadata?
  if (!var %in% colnames(object_to_test@meta.data)) {
    message(paste("Skipping: Variable", var, "not found."))
    next
  }
  
  # Set the identity to the variable of interest
  Idents(object_to_test) <- var
  
  # Check 2: Are there at least 2 groups?
  unique_groups <- unique(Idents(object_to_test))
  unique_groups <- unique_groups[!is.na(unique_groups)] # Remove NAs
  
  if (length(unique_groups) < 2) {
    message(paste("Skipping: Only 1 group found for", var, "(", unique_groups, ")"))
    next
  }
  
  # C. Run DE Analysis
  tryCatch({
    markers <- FindAllMarkers(
      object = object_to_test,
      only.pos = FALSE, 
      verbose = FALSE,
      logfc.threshold = 0.1,
      min.pct = 0.01,
      min.diff.pct = -Inf  # Keep Seurat defaults or your specific params
    )
    
    # D. Add Custom Metrics
    if (nrow(markers) > 0) {
      # Calculate p_FC
      markers$p_FC <- (1 - markers$p_val) * markers$avg_log2FC
      
      # Calculate difference in proportion
      markers$adj_diff <- abs(markers$pct.1 - markers$pct.2) / pmax(markers$pct.1, markers$pct.2)
      
      # Calculate combined score
      markers$diff_pfc <- markers$adj_diff * markers$p_FC
      
      # Save to list (Key is just the variable name now)
      global_de_results[[var]] <- markers
      message(paste("Found", nrow(markers), "markers."))
      
    } else {
      message("No significant markers found.")
    }
    
  }, error = function(e) {
    message(paste("Error running DE for", var, ":", e$message))
  })
}

# Optional: Save the global results

# plot number of DE
# 1. Initialize empty list for summary data
plot_data_list <- list()

# 2. Loop through the results to extract counts
# 'de_results' is the list you created in the previous step
for (var_name in names(global_de_results)) {
  
  # Get the results dataframe for this variable
  df <- global_de_results[[var_name]]
  df <- df[df$p_val_adj < 0.05, ]
  if (!is.null(df) && nrow(df) > 0) {
    # Count the number of markers per group (stored in 'cluster' column)
    summary_df <- df %>%
      group_by(cluster) %>%
      summarise(Count = n()) %>%
      mutate(Variable = var_name) # Add label for the variable type
    
    plot_data_list[[var_name]] <- summary_df
  }
}
saveRDS(global_de_results, file = paste0(output_dir, basename(iso_expression_h5_dir), "/global_DE_results.rds"))

# 3. Combine into one master dataframe
plot_data <- do.call(rbind, plot_data_list)

# 4. Create the Plot
# We use 'facet_wrap' to create separate panels for Sex, Age, and ROI
p <- ggplot(plot_data, aes(x = cluster, y = Count, fill = cluster)) +
  geom_col() +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) + # Add numbers on top
  facet_wrap(~Variable, scales = "free_x", strip.position = "top") + # Separate panels
  theme_classic() +
  labs(
    title = "Number of Significant DE Genes by Criteria",
    subtitle = "Non-neuronal DE result",
    x = "Group",
    y = "Number of DE Genes",
    fill = "Group"
  ) +
  theme(
    legend.position = "none", # Hide legend (colors are self-explanatory)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.background = element_rect(fill = "lightgrey"), # Style the panel headers
    strip.text = element_text(face = "bold", size = 12)
  )

# 5. Save and Print
print(p)
ggsave(paste0(output_dir, basename(iso_expression_h5_dir),"/Global_DE_Counts.png"), p, width = 10, height = 6)

#################
# Subclass-Specific DE Analysis
# Goal: Find DE genes for Sex, Age, and ROI within each Subclass

# 1. Define the groups we want to loop through
target_subclasses <- unique(`iso_seurat`$hicat_cluster)
variables_to_test <- c("sex", "age_cat", "roi")

# 2. Initialize a list to store results
# Structure: de_results$Subclass_Name$Variable_Name
de_results <- list()
num_genes_expr <- c()
# 3. Main Loop
for (subclass in target_subclasses) {
  
  message(paste0("\nProcessing Subclass: ", subclass, " ------------------"))
  
  # A. Subset the object to just this subclass
  # We use subset() to isolate the cells
  sub_obj <- subset(`iso_seurat`, subset = hicat_cluster == subclass)
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
      markers <- FindAllMarkers(
        object = sub_obj,
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

saveRDS(de_results, file = paste0(output_dir, basename(iso_expression_h5_dir),"/DE_results.rds"))

#################
# 1. Extract the number of DE genes from the 'de_results' list
# Initialize an empty list to store the counts

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
plot_data$Variable <- factor(plot_data$Variable,
                             levels = c("sex", "age_cat", "roi"),
                             labels = c("Sex", "Age", "ROI"))

# 3. Create the Plot
text <- paste0(basename(gsub("\\..*$", "", iso_expression_h5_dir)), ": # Sig. DE Genes per Hicat Cluster")


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
ggsave(paste0(output_dir, basename(iso_expression_h5_dir),"/DE_Counts_Summary_Hicat.png"), p, width = 17, height = 8)

text <- paste0(basename(gsub("\\..*$", "", iso_expression_h5_dir)), ": Prop. Sig. DE Genes per Hicat Cluster")


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
ggsave(paste0(output_dir, basename(iso_expression_h5_dir),"/DE_Prop_Summary_Hicat.png"), p, width = 17, height = 8)