# subset the total isocortex into non-neuronal, sex, age, and region specific files
'''
 it would be a better idea to save the cell ids for each of these subsets instead
 of their entire expression matrices, but I'm lazy and storage is not an issue for
 me right now :)
'''
import numpy as np
import pandas as pd
import anndata as ad
import h5py as h5
from scipy import sparse

f = ad.read_h5ad('C:/Users/megan_ding/Desktop/BCB430/Isocortex.h5ad', 'r')
iso_data = f[(np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
pfc_data = f[(np.isin(f.obs.roi, ["PL-ILA-ORB", "AI-CLA", "ACA"])) & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
female_data = f[(f.obs.sex == "F") & (np.isin(f.obs.roi, ["PL-ILA-ORB", "AI-CLA", "ACA"])) & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
male_data = f[(f.obs.sex == "M") & (np.isin(f.obs.roi, ["PL-ILA-ORB", "AI-CLA", "ACA"])) & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
adult_data = f[(f.obs.age_cat == "adult") & (np.isin(f.obs.roi, ["PL-ILA-ORB", "AI-CLA", "ACA"])) & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
aged_data = f[(f.obs.age_cat == "aged") & (np.isin(f.obs.roi, ["PL-ILA-ORB", "AI-CLA", "ACA"])) & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
female_adult_data = f[(f.obs.sex == "F") & (f.obs.age_cat == "adult") & (np.isin(f.obs.roi, ["PL-ILA-ORB", "AI-CLA", "ACA"])) & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
female_aged_data = f[(f.obs.sex == "F") & (f.obs.age_cat == "aged") & (np.isin(f.obs.roi, ["PL-ILA-ORB", "AI-CLA", "ACA"])) & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
male_adult_data = f[(f.obs.sex == "M") & (f.obs.age_cat == "adult") & (np.isin(f.obs.roi, ["PL-ILA-ORB", "AI-CLA", "ACA"])) & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
male_aged_data = f[(f.obs.sex == "M") & (f.obs.age_cat == "aged") & (np.isin(f.obs.roi, ["PL-ILA-ORB", "AI-CLA", "ACA"])) & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
iso_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/isocortex.h5ad', compression="gzip")
pfc_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC.h5ad', compression="gzip")
pfc_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC.h5ad', compression="gzip")
female_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC_female.h5ad', compression="gzip")
male_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC_male.h5ad', compression="gzip")
adult_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC_adult.h5ad', compression="gzip")
aged_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC_aged.h5ad', compression="gzip")
female_adult_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC_female_adult.h5ad', compression="gzip")
female_aged_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC_female_aged.h5ad', compression="gzip")
male_adult_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC_male_adult.h5ad', compression="gzip")
male_aged_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/PFC_male_aged.h5ad', compression="gzip")

female_data = f[(f.obs.sex == "F") & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
male_data = f[(f.obs.sex == "M") & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
adult_data = f[(f.obs.age_cat == "adult") & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
aged_data = f[(f.obs.age_cat == "aged") & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
female_adult_data = f[(f.obs.sex == "F") & (f.obs.age_cat == "adult") & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
female_aged_data = f[(f.obs.sex == "F") & (f.obs.age_cat == "aged") & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
male_adult_data = f[(f.obs.sex == "M") & (f.obs.age_cat == "adult") & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]
male_aged_data = f[(f.obs.sex == "M") & (f.obs.age_cat == "aged") & (np.isin(f.obs.class_label, ["Astro-Epen", "Vascular", "Immune", "OPC-Oligo"]))]

female_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/isocortex_female.h5ad', compression="gzip")
male_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/isocortex_male.h5ad', compression="gzip")
adult_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/isocortex_adult.h5ad', compression="gzip")
aged_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/isocortex_aged.h5ad', compression="gzip")
female_adult_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/isocortex_female_adult.h5ad', compression="gzip")
female_aged_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/isocortex_female_aged.h5ad', compression="gzip")
male_adult_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/isocortex_male_adult.h5ad', compression="gzip")
male_aged_data.write('C:/Users/megan_ding/Desktop/BCB430/split_isocortex_data/isocortex_male_aged.h5ad', compression="gzip")
