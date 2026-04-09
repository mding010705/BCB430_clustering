#!/usr/bin/env bash

# run DEA script between subset clusters

filepath='C:/Users/megan_ding/Desktop/BCB430'
N=17
(
# iso vs. pfc
# iso male vs. iso female, pfc male vs. pfc female
# iso adult vs. iso aged, pfc adult vs. pfc aged
# iso male adult vs. iso female adult, pfc male adult vs. pfc female adult
# iso male aged vs. iso female aged, pfc male aged vs. pfc female aged
# iso male adult vs. iso female aged, iso female adult vs. iso male aged
# pfc male adult vs. pfc female aged, pfc female adult vs. pfc male aged
lis1=( "PFC.h5ad" "isocortex_female.h5ad" "PFC_female.h5ad" "isocortex_aged.h5ad" \
"PFC_aged.h5ad" "isocortex_male_aged.h5ad" "PFC_male_aged.h5ad" \
"isocortex_female_aged.h5ad" "PFC_female_aged.h5ad" \
"isocortex_female_aged.h5ad" "PFC_female_aged.h5ad" \
"isocortex_female_adult.h5ad" "PFC_female_adult.h5ad" \
"isocortex_female_adult.h5ad" "PFC_female_adult.h5ad" \
"isocortex_female_aged.h5ad" "PFC_female_aged.h5ad" )
lis2=( "isocortex.h5ad" "isocortex_male.h5ad" "PFC_male.h5ad" "isocortex_adult.h5ad" \
"PFC_adult.h5ad" "isocortex_male_adult.h5ad" "PFC_male_adult.h5ad" \
"isocortex_female_adult.h5ad" "PFC_female_adult.h5ad" \
"isocortex_male_aged.h5ad" "PFC_male_aged.h5ad" \
"isocortex_male_adult.h5ad" "PFC_male_adult.h5ad" \
"isocortex_male_aged.h5ad" "PFC_male_aged.h5ad" \
"isocortex_male_adult.h5ad" "PFC_male_adult.h5ad" )

for index in "${!lis1[@]}"; do
  ((j=j%N)); ((j++==0)) && wait
  echo $index
  Rscript ${filepath}/between_cluster_dea.R ${lis1[$index]} ${lis2[$index]}
done
)