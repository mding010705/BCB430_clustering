#!/usr/bin/env bash

# run directionally stratified MDD gene list overlap enrichment for selected cell types

filepath='C:/Users/megan_ding/Desktop/BCB430'
N=90
(
for i in ${filepath}/cluster/*.h5ad; do
  for k in {1..2}; do
    ((j=j%N)); ((j++==0)) && wait
    echo $i
    updown="up"
    if (( k % 2 == 1 )); then
        updown="down"
    fi
    Rscript ${filepath}/enrichment_heatmap_sel_cells.R $i $updown &
  done
done
)