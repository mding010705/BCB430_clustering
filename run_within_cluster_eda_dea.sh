#!/usr/bin/env bash

# run EDA and DEA within cluster script

filepath='C:/Users/megan_ding/Desktop/BCB430'
N=27
(
for i in ${filepath}/split_isocortex_data/*.h5ad; do
  ((j=j%N)); ((j++==0)) && wait
  echo $i
  mkdir -p ${filepath}/cluster/$(basename "$i")
  Rscript ${filepath}/within_cluster_eda_dea.R $i &
done
)