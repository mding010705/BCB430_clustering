#!/usr/bin/env bash

# run directionally unstratified MDD gene list overlap enrichment

filepath='C:/Users/megan_ding/Desktop/BCB430'
N=45
(
for i in ${filepath}/cluster/*.h5ad; do
  ((j=j%N)); ((j++==0)) && wait
  echo $i
  Rscript ${filepath}/mdd_enrichment_btwn_clust.R $i &
done
)