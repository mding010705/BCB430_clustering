#!/usr/bin/env bash

# run cluster cell type labelling

filepath='C:/Users/megan_ding/Desktop/BCB430'
N=27
(
for i in ${filepath}/clust_res_hicat/merged/*.rds; do
  ((j=j%N)); ((j++==0)) && wait
  echo $i
  Rscript ${filepath}/label_clusters.R $i &
done
)