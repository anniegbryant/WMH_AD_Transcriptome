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

# Load the DESeq2 results
load(paste0(study, "AWS_vs_PWS_dds_res.Rdata"))

# Load the metadata
metadata <- read.csv(paste0(study, "human_sample_info.csv")) %>%
  arrange(Sample_Name) %>%
  mutate(Group = factor(Group, levels=c("AD", "CTRL")))

# Filter DESeq2 results to BH-adjusted p-value < 0.05 and |lgo2FC| >= 0.5
AD_bulk_dds_df <- as.data.frame(AD_bulk_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Group = "AD",
         Tissue = "Bulk")

CTRL_bulk_dds_df <- as.data.frame(CTRL_bulk_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Group = "CTRL",
         Tissue = "Bulk")

AD_bv_dds_df <- as.data.frame(AD_bv_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Group = "AD",
         Tissue = "BV")

CTRL_bv_dds_df <- as.data.frame(CTRL_bv_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Group = "CTRL",
         Tissue = "BV")

# Combine the results
AWS_vs_PWS_res <- do.call(plyr::rbind.fill, list(AD_bulk_dds_df,
                                                      CTRL_bulk_dds_df,
                                                      AD_bv_dds_df,
                                                      CTRL_bv_dds_df))

# Count the number of significant DEGs in each region and tissue
significant_DEGs <- AWS_vs_PWS_res %>% 
  filter(abs(log2FoldChange)>0.5, padj < 0.05)

significant_DEGs %>% 
  mutate(Direction = ifelse(log2FoldChange>0, "Up", "Down")) %>% 
  group_by(Region, Tissue, Direction) %>%
  summarise(count = n())

# Save the significant DEGs
for (region in c("AWS", "PWS")) {
  for (tissue in c("Bulk", "BV")) {
    AD_vs_CTRL_res %>%
      filter(Region==region, Tissue==tissue) %>%
      dplyr::select(Gene, log2FoldChange, pvalue, padj) %>%
      arrange(desc(log2FoldChange)) %>%
      write.csv(., paste0(study_path, "AD_vs_CTRL_", region, "_", tissue, "_DE_Genes.csv"),
                row.names=F)
  }
}