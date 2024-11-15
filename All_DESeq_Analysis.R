# All DESeq2 analysis for WMH AD transcriptome project
# Author: Annie G. Bryant
# Updated: 15 November 2024

library(tidyverse)
library(DESeq2)
library(cowplot)
theme_set(cowplot::theme_cowplot())
library(argparse)

# Set up command line arguments
parser <- ArgumentParser()
parser$add_argument("--study_path", help="Path to study directory")

# Parse the arguments
args <- parser$parse_args()
study_path <- args$study_path

# Output directory for plots
plot_dir <- paste0(study_path, "plots/")

# Load Salmon-processed data
load(paste0(study_path, "human_WMH_RNAseq_salmon.Rdata"))

# Read in metadata
metadata <- read.csv(paste0(study_path, "human_sample_info.csv")) %>%
  arrange(Sample_Name) %>%
  mutate(Group = factor(Group, levels=c("CTRL", "AD")))
rownames(metadata) <- colnames(human_WMH_RNAseq_salmon$counts)

# Helper function that runs DESeq2 for a given tissue type (bulk or blood vessel)
# And brain region (AWS or PWS)
run_deseq2_for_region_tissue <- function(region, tissue, dds_data, metadata, salmon_RNAseq_data) {
  tissue_full <- ifelse(tissue=="bulk", "Bulk", "Blood_Vessel")
  
  # Extract metadata
  region_tissue_meta <- metadata %>% filter(Tissue==tissue_full, 
                                            White_Matter_Region==region)
  
  # Get the gene names
  region_tissue_rows <- colnames(salmon_RNAseq_data$counts) %in% rownames(region_tissue_meta)
  
  # Only use non-missing data
  dds_use <- dds_data[!rownames(dds_data) == "", region_tissue_rows]

  # Fit DESeq with all default parameters, the design formula is already set in the dds_by_group object
  dds <- DESeq(dds_use)
  
  return(dds)
}

################### AD vs. CTRL ###################

# Apply the helper function to each region/tissue combination
if (!file.exists(paste0(study_path, "AD_vs_CTRL_dds_res_BH.Rdata"))) {
  # Set up gene expression data for differential expression analysis by group (AD vs. CTRL)
  dds_by_group <- DESeqDataSetFromTximport(human_WMH_RNAseq_salmon, metadata, ~ Group)

  # Iterate over regions and tissue types
  for (region in c("AWS", "PWS")) {
    for (tissue in c("bulk", "bv")) {
      # Run DESeq2 for the corresponding region/tissue combination
      dds <- run_deseq2_for_region_tissue(region=region, tissue=tissue, dds_data=dds_by_group, metadata=metadata, salmon_RNAseq_data=human_WMH_RNAseq_salmon)
      assign(paste0(region, "_", tissue, "_dds"), dds)
      # Extract the results of the differential expression analysis with default (Benjamini-Hochberg) correction
      dds_res <- results(dds, independentFiltering=FALSE, pAdjustMethod="BH")
      assign(paste0(region, "_", tissue, "_dds_res"), dds_res)
    }
  }

  # Save the DESeq2 datasets and results
  save(AWS_bulk_dds_res, AWS_bulk_dds, 
       AWS_bv_dds_res, AWS_bv_dds,
       PWS_bulk_dds_res, PWS_bulk_dds, 
       PWS_bv_dds_res, PWS_bv_dds,
       file=paste0(study_path, "AD_vs_CTRL_dds_res.Rdata"))
}

################### Frontal vs. Occipital White Matter ###################

# Helper function that runs DESeq2 for a given tissue type (bulk or blood vessel)
# And diagnostic group (AD or CTRL)
run_deseq2_for_dx_tissue <- function(dx, tissue, dds_data, metadata, salmon_RNAseq_data) {
  tissue_full <- ifelse(tissue=="bulk", "Bulk", "Blood_Vessel")
  
  # Extract metadata
  dx_tissue_meta <- metadata %>% filter(Tissue==tissue_full, 
                                            Group==dx)
  
  # Get the gene names
  dx_tissue_rows <- colnames(salmon_RNAseq_data$counts) %in% rownames(dx_tissue_meta)
  
  # Only use non-missing data
  dds_use <- dds_data[!rownames(dds_data) == "", dx_tissue_rows]

  # Fit DESeq with all default parameters, the design formula is already set in the dds_data object
  # In this case, it would correspond to ~ Region
  dds <- DESeq(dds_use)
  
  return(dds)
}

# Apply the helper function to each dx/tissue combination
if (!file.exists(paste0(study_path, "AWS_vs_PWS_dds_res_BH.Rdata"))) {
  # Set up gene expression data for differential expression analysis by group (AD vs. CTRL)
  dds_by_group <- DESeqDataSetFromTximport(human_WMH_RNAseq_salmon, metadata, ~ White_Matter_Region)
  # Iterate over regions and tissue types
  for (dx in c("AD", "CTRL")) {
    for (tissue in c("bulk", "bv")) {
      # Run DESeq2 for the corresponding region/tissue combination
      dds <- run_deseq2_for_dx_tissue(dx=dx, tissue=tissue, dds_data=dds_by_group, metadata=metadata, salmon_RNAseq_data=human_WMH_RNAseq_salmon)
      assign(paste0(dx, "_", tissue, "_dds"), dds)
      # Extract the results of the differential expression analysis with BH correction
      dds_res <- results(dds, independentFiltering=FALSE, pAdjustMethod="BH")
      assign(paste0(dx, "_", tissue, "_dds_res"), dds_res)
    }
  }
  # Save the DESeq2 datasets and results
  save(AD_bulk_dds_res, AD_bulk_dds, 
       AD_bv_dds_res, AD_bv_dds,
       CTRL_bulk_dds_res, CTRL_bulk_dds, 
       CTRL_bv_dds_res, CTRL_bv_dds,
       file=paste0(study_path, "AWS_vs_PWS_dds_res.Rdata"))
}

################### High vs. Low WMH Differential Expression ###################

# Set up gene expression data for differential expression analysis by WMH Level (High vs. Low)
dds_by_group <- DESeqDataSetFromTximport(human_WMH_RNAseq_salmon_volumes, metadata_WMH, ~ WMH_Level)

run_deseq2_for_region_tissue_WMH <- function(region, tissue, dds_data, metadata, salmon_RNAseq_data) {
  tissue_full <- ifelse(tissue=="bulk", "Bulk", "Blood_Vessel")
  
  # Extract metadata, only for samples with WMH volumes
  region_tissue_meta <- metadata %>% filter(Tissue==tissue_full, 
                                            White_Matter_Region==region)
  
  # Get the gene names
  region_tissue_WMH_rows <- colnames(salmon_RNAseq_data$counts) %in% rownames(region_tissue_meta)
  
  # Only use non-missing data
  dds_use <- dds_data[!rownames(dds_data) == "", region_tissue_WMH_rows]
  
  # Fit DESeq with all default parameters, the design formula is already set in the dds_by_group object
  dds <- DESeq(dds_use)
  
  return(dds)
}

# Apply the helper function to each region/tissue combination
if (!file.exists(paste0(study_path, "High_vs_Low_WMH_dds_res.Rdata"))) {
  # Filter metadata to only include samples with WMH data
  metadata_WMH <- metadata %>%   
    filter(!is.na(WMH_mm3)) %>%
    mutate(WMH_Level = case_when(WMH_mm3 > 20000 ~ "High",
                                 T ~ "Low")) %>% 
    mutate(WMH_Level = factor(WMH_Level, levels = c("Low", "High")))
  
  # Filter salmon object
  samples_without_WMH <- metadata %>% pull(WMH_mm3) %>% is.na()
  human_WMH_RNAseq_salmon_volumes <- human_WMH_RNAseq_salmon
  human_WMH_RNAseq_salmon_volumes$abundance <- human_WMH_RNAseq_salmon_volumes$abundance[,!samples_without_WMH]
  human_WMH_RNAseq_salmon_volumes$counts <- human_WMH_RNAseq_salmon_volumes$counts[,!samples_without_WMH]
  human_WMH_RNAseq_salmon_volumes$length <- human_WMH_RNAseq_salmon_volumes$length[,!samples_without_WMH]
  
  # Set up gene expression data for differential expression analysis by WMH Level (High vs. Low)
  dds_by_group <- DESeqDataSetFromTximport(human_WMH_RNAseq_salmon_volumes, metadata_WMH, ~ WMH_Level)
  
  # Iterate over regions and tissue types
  for (region in c("AWS", "PWS")) {
    for (tissue in c("bulk", "bv")) {
      # Run DESeq2 for the corresponding region/tissue combination
      dds <- run_deseq2_for_region_tissue(region=region, tissue=tissue, dds_data=dds_by_group, metadata=metadata_WMH, salmon_RNAseq_data=human_WMH_RNAseq_salmon_volumes)
      assign(paste0(region, "_", tissue, "_dds_high_low"), dds)
      # Extract the results of the differential expression analysis with default (Benjamini-Hochberg) correction
      dds_res <- results(dds, independentFiltering=FALSE, pAdjustMethod="BH")
      assign(paste0(region, "_", tissue, "_dds_res_high_low"), dds_res)
    }
  }
  
  # Save the DESeq2 datasets and results
  save(AWS_bulk_dds_res_high_low, AWS_bulk_dds_high_low, AWS_bv_dds_res_high_low, AWS_bv_dds_high_low,
       PWS_bulk_dds_res_high_low, PWS_bulk_dds_high_low, PWS_bv_dds_res_high_low, PWS_bv_dds_high_low,
       file=paste0(study_path, "High_vs_Low_WMH_dds_res.Rdata"))
}

################### Bulk tissue vs. Blood Vessel Analysis ###################

# Helper function that runs DESeq2 for a given brain region (AWS or PWS)
# And diagnostic group (AD or CTRL)
run_deseq2_for_dx_region <- function(dx, region, dds_data, metadata, salmon_RNAseq_data) {
  # Extract metadata
  dx_region_meta <- metadata %>% filter(White_Matter_Region==region, Group==dx)
  
  # Get the gene names
  dx_region_rows <- colnames(salmon_RNAseq_data$counts) %in% rownames(dx_region_meta)
  
  # Only use non-missing data
  dds_use <- dds_data[!rownames(dds_data) == "", dx_region_rows]
  dds <- DESeq(dds_use)
  
  return(dds)
}

# Apply the helper function to each dx/region combination
if (!file.exists(paste0(study_path, "Bulk_vs_BV_dds_res.Rdata"))) {
  # Set up gene expression data for differential expression analysis by group (AD vs. CTRL)
  dds_by_group <- DESeqDataSetFromTximport(human_WMH_RNAseq_salmon, metadata, ~ Tissue)
  # Iterate over regions and tissue types
  for (dx in c("AD", "CTRL")) {
    for (region in c("AWS", "PWS")) {
      # Run DESeq2 for the corresponding dx/region combination
      dds <- run_deseq2_for_dx_region(dx=dx, region=region, dds_data=dds_by_group, metadata=metadata, salmon_RNAseq_data=human_WMH_RNAseq_salmon)
      assign(paste0(dx, "_", region, "_dds"), dds)
      # Extract the results of the differential expression analysis with BH correction
      dds_res <- results(dds, independentFiltering=FALSE, pAdjustMethod="BH")
      assign(paste0(dx, "_", region, "_dds_res"), dds_res)
    }
  }

  # Save the DESeq2 datasets and results
  save(AD_AWS_dds_res, AD_AWS_dds, 
       AD_PWS_dds_res, AD_PWS_dds,
       CTRL_AWS_dds_res, CTRL_AWS_dds, 
       CTRL_PWS_dds_res, CTRL_PWS_dds,
       file=paste0(study_path, "Bulk_vs_BV_dds_res.Rdata"))
} 