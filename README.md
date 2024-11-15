# Code to reproduce differential expression analysis in *Molecular profiling of frontal and occipital subcortical white matter hyperintensities in Alzheimer's disease*

This repository includes the code needed to reproduce our differential expression analysis and other robustness analyses for our preprint, [*Molecular profiling of frontal and occipital subcortical white matter hyperintensities in Alzheimer's disease*](https://www.biorxiv.org/content/10.1101/2024.06.13.598845).
All scripts are written in R and designed to be run from the command line with the argument `--study_path [\path\to\salmon\processed\data]`, with generated plots saved within a `plots/` folder in the designated `study_path`. 
Users should quantify transcripts from the corresponding RNA-seq data using the [Salmon software](https://combine-lab.github.io/salmon/).

## Running DESeq2 

The first step to reproduce our differential expression analysis is to run the script `All_DESeq_Analysis.R`, which will run the `DESeq` function with default parameters for each of our contrasts of interest: 

1. AD vs. CTRL: The model is defined as `~ Group`, where Group is a metadata column corresponding to 'AD' or 'CTRL' 
2. AWS vs. PWS: The model is defined as `~ White_Matter_Region`, where White_Matter_Region is a metadata column corresponding to 'AWS' or 'PWS' 
3. High vs. Low WMH: The model is defined as `~ WMH_Level`, where WMH_Level is a metadata column corresponding to 'High' (>=20000mm3) or 'Low' (<20000mm3) 
2. Bulk tissue vs. blood vessels: The model is defined as `~ Tissue`, where Tissue is a metadata column corresponding to 'Bulk' or 'Blood_Vessel' 

## Generating visuals

The core visual results from our differential expression analysis can be reproduced as follows: 

* AD vs. CTRL: `AD_vs_Control_DEG_Visualization.R` and `PCA_for_WMH_RNAseq.R`
* AWS vs. PWS: `AWS_vs_PWS_DEG_Visualization.R` 
* High vs. Low WMH: `Compare_WMH_volumes.R` 
* Bulk tissue vs. blood vessels: `Bulk_vs_BV_CTRL_DEG_Visualization_and_EWCE.R` 

Note that the last one also includes [expression-weighted cell type enrichment (EWCE) analysis](https://github.com/NathanSkene/EWCE). 


## Robustness analyses 

To contextualize our findings with those described in prior vascular tissue sequencing analysis in AD, we compared our top AD-associated differentially expressed genes (DEGs) with those identified by Bryant et al. (2023) and Mitroi et al. (2022).
The DEG data from these two studies is included in this repository with the file names `Bryant_2023_Extended_Data_Table_7_2.csv`, `Mitroi_2022_VaD_vs_NC.csv`, and `Mitroi_2022_VaDadj_vs_NC.csv`. 
We visualized overlap among these gene sets in `Compare_our_DEGs_with_literature.R`.