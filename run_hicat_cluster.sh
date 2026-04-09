#!/usr/bin/env bash

# run hicat clustering on all isocortex subsets

filepath='C:/Users/megan_ding/Desktop/BCB430'
N=27
(
for i in ${filepath}/split_isocortex_data/*.h5ad; do
  ((j=j%N)); ((j++==0)) && wait
  echo $i
  Rscript ${filepath}/hicat_cluster_job.R $i &
done
)