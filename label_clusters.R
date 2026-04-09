# heuristically transfer cell type labels from metadata to subset clusters

library(dplyr)


stratum <- commandArgs(trailingOnly = TRUE)[1]

cls <- readRDS(stratum)
cl <- cls$cl
metadata_iso <- read.csv("C:/Users/megan_ding/Desktop/BCB430/metadata_isocortex.csv")

Mode <- function(x) {
  # Remove NAs if present
  x <- na.omit(x)
  # Calculate frequency table
  tab <- table(x)
  # Find the name of the value with the maximum frequency
  mode_val <- names(tab)[which.max(tab)]
  
  # Convert back to numeric if the original input was numeric
  if (is.numeric(x)) {
    return(as.numeric(mode_val))
  } else {
    return(mode_val)
  }
}

tab <- sort(table(cl), decreasing = T)
for (c in names(tab)){
  cells <- names(cl)[cl == c]
  ct <- metadata_iso[metadata_iso$X %in% cells, c("subclass_label", "supertype_label", 
                                                  "cluster_label")]
  newc <- NULL
  consensus <- ct %>%
      summarise(across(everything(), Mode))
    if(is.null(newc)){
      newc <- paste0(consensus$cluster_label, "_new_", c)
    }
    if (!("subclass_label" %in% colnames(consensus))){
      consensus$subclass_label <- NA
    }
    if (!("supertype_label" %in% colnames(consensus))){
      consensus$supertype_label <- NA
    }
    if (!("cluster_label" %in% colnames(consensus))){
      cls$cl[cls$cl == c] <- newc
    }
    
    for (l in consensus){
      if (is.na(l)){
        next
      }
      if (!(l %in% cls$cl)){
        cls$cl[cls$cl == c] <- l
      }
    }
  if (any(cls$cl == c)){
    cls$cl[cls$cl == c] <- newc
  }
  
}

saveRDS(cls, file = paste0("C:/Users/megan_ding/Desktop/BCB430/clust_res_hicat/labeled/", 
                           basename(gsub("\\..*$", "", stratum)), ".rds"))
