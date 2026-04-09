if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("MatrixExtra", quietly = TRUE)) install.packages("MatrixExtra")
if (!requireNamespace("bluster", quietly = TRUE)) BiocManager::install("bluster")
if (!requireNamespace("ggbeeswarm", quietly = TRUE)) BiocManager::install("ggbeeswarm")

library(tidyverse)
library(anndata)
library(MatrixExtra)
library(bluster)
library(irlba)

findElbowPoint <- function(variance) {
  if (is.unsorted(-variance)) {
    stop("'variance' should be sorted in decreasing order")
  }
  
  # Finding distance from each point on the curve to the diagonal.
  dy <- -diff(range(variance))
  dx <- length(variance) - 1
  l2 <- sqrt(dx^2 + dy^2)
  dx <- dx/l2
  dy <- dy/l2
  
  dy0 <- variance - variance[1]
  dx0 <- seq_along(variance) - 1
  
  parallel.l2 <- sqrt((dx0 * dx)^2 + (dy0 * dy)^2)
  normal.x <- dx0 - dx * parallel.l2
  normal.y <- dy0 - dy * parallel.l2
  normal.l2 <- sqrt(normal.x^2 + normal.y^2)
  
  #Picking the maximum normal that lies below the line.
  #If the entire curve is above the line, we just pick the last point.
  below.line <- normal.x < 0 & normal.y < 0
  if (!any(below.line)) {
    length(variance)
  } else {
    which(below.line)[which.max(normal.l2[below.line])]
  }
}
library(ggplot2)
library(ggbeeswarm)

subset_data_dir <-  "C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/"
clust_dir <- "C:/Users/megan_ding/Desktop/BCB430/clust_res_hicat/labeled/"
plot_out_dir <- "C:/Users/megan_ding/Desktop/BCB430/cluster/"

silhouette_beeswarm <- function(file_prefix = "isocortex"){
  big_data <- read_h5ad(paste0(subset_data_dir, file_prefix, ".h5ad"))
  ex <- as.matrix(t(big_data$X))
  max.dim <- 50
  pca_dat <- prcomp_irlba(ex, n = min(max.dim, dim(ex) - 1));
  
  dim.elbow = findElbowPoint(pca_dat$sdev^2)
  dim.elbow = min(dim.elbow, ncol(pca_dat$rotation))
  
  sil.approx <- approxSilhouette(pca_dat$rotation[, 1:dim.elbow], clusters = clust[[file_prefix]]$cl)
  sil.data <- as.data.frame(sil.approx)
  sil.data$closest <- factor(ifelse(sil.data$width > 0, clust[[file_prefix]]$cl, sil.data$other))
  sil.data$cluster <- clust[[file_prefix]]$cl
  
  
  means <- aggregate(width ~ cluster, sil.data, mean)
  
  p <- ggplot(sil.data, aes(x=as.factor(cluster), y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley") +
    stat_summary(fun=mean, colour="darkred", geom="point", 
                 shape=18, size=3, show.legend=FALSE) + 
    stat_summary(fun=mean, colour="black", geom="text", show.legend=FALSE,
                 vjust=-0.7, aes( label=round(..y.., digits=3))) +
    ggtitle(paste(file_prefix, "silhouette scores"))
  ggsave(paste0(plot_out_dir, file_prefix, 
                ".h5ad/silhouette_score_beeswarm_hicat.png"),
         p, width = 15, height = 8, bg = "white")
  
  # table(Cluster=clust[[file_prefix]]$cl, sil.data$closest)
  return(sil.data)
}



files <- list.files(path = clust_dir,
                    full.names = TRUE)


clust <- list()
for (f in files){
  clust[[basename(gsub(pattern = "\\..*", replacement = "", f))]] <- readRDS(f)
}

sil_scores <- list()
for (n in names(clust)){
  sil_scores[[n]] <- silhouette_beeswarm(n)
}
saveRDS(sil_scores, file = "C:/Users/megan_ding/Desktop/BCB430/sil_scores_hicat.rds")
