# cluster subset of isocortex data using hicat

library(reticulate)
library(anndata)
library(Matrix)
library(matrixStats)
library(sparseMatrixStats)

setAs("dgRMatrix", to = "dgCMatrix", function(from){
  as(as(from, "CsparseMatrix"), "dgCMatrix")
})

read_h5ad_bigdat <- function(h5ad.file = "C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC_female.h5ad") {
  sce <- anndata::read_h5ad(h5ad.file)
  
  list(
    expr = as(t(sce$X), "dgCMatrix"),
    row_id = sce$var_names,
    col_id = rownames(sce$X)
  )
}

intermediate_step_dir <- "C:/Users/megan_ding/Desktop/BCB430/clust_inter/"
output_dir <- "C:/Users/megan_ding/Desktop/BCB430/clust_res_hicat/"

file_names <- list.files(path = "./scrattch.hicat/R", pattern = "\\.[Rr]$", full.names = TRUE)

# Use lapply to apply the source function to each file in the list
lapply(file_names, source)

args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
print(paste0(timestamp(), fn))


dat <- read_h5ad_bigdat(fn)

de.param = de_param(q1.th=0.4, 
                    q.diff.th = 0.7, 
                    de.score.th=300, 
                    min.cells=50)
# issues with heatmap3 script from the hicat package, so we define our own
source("./new_heatmap3.R")

set.seed(430)
png(filename = paste0(output_dir, 
                      basename(fn), "_cluster_marker_heatmap.png"), 
                      width = 2400, height = 2400, units = "px")

res <- iter_clust(norm.dat = dat$expr, 
                  select.cells = dat$col_id, 
                  de.param =de.param, 
                  prefix=paste0(intermediate_step_dir,  basename(fn)), 
                  max.cl.size=3000, 
                  split.size = 50, 
                  verbose=1, 
                  sampleSize=50000)
dev.off()
saveRDS(res, paste0(output_dir,"/unmerged/", 
                    basename(fn), ".rds"))

# Sample cells for merging
sampled.cells = sample_cells(res$cl, min(50000, length(dat$col_id)))


merge.res = merge_cl(norm.dat=dat$expr[, sampled.cells], 
                        cl=res$cl, 
                        rd.dat.t=dat$expr[res$markers, sampled.cells], 
                        de.param=de.param, 
                        verbose=TRUE)

saveRDS(merge.res, paste0(output_dir, "/merged/", 
                    basename(fn), ".rds"))

print(paste0(timestamp(), fn))
