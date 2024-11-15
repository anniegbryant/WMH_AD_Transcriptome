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
load(paste0(study, "AD_vs_CTRL_dds_res.Rdata"))

# Load the metadata
metadata <- read.csv(paste0(study, "human_sample_info.csv")) %>%
  arrange(Sample_Name) %>%
  mutate(Group = factor(Group, levels=c("AD", "CTRL")))

# Filter DESeq2 results to BH-adjusted p-value < 0.05 and |lgo2FC| >= 0.5
AWS_bulk_dds_df <- as.data.frame(AWS_bulk_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Region = "AWS",
         Tissue = "Bulk")
PWS_bulk_dds_df <- as.data.frame(PWS_bulk_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Region = "PWS",
         Tissue = "Bulk")
AWS_bv_dds_df <- as.data.frame(AWS_bv_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Region = "AWS",
         Tissue = "BV")
PWS_bv_dds_df <- as.data.frame(PWS_bv_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Region = "PWS",
         Tissue = "BV")

# Combine the results
AD_vs_CTRL_res <- do.call(plyr::rbind.fill, list(AWS_bulk_dds_df,
                                                      PWS_bulk_dds_df,
                                                      AWS_bv_dds_df,
                                                      PWS_bv_dds_df))

# Count the number of significant DEGs in each region and tissue
significant_DEGs <- AD_vs_CTRL_res %>% 
  filter(abs(log2FoldChange)>0.5, padj < 0.05)

significant_DEGs %>% 
  mutate(Direction = ifelse(log2FoldChange>0, "Up", "Down")) %>% 
  group_by(Region, Tissue, Direction) %>%
  summarise(count = n())

# Plot Euler diagrams:

# Upregulated genes
AD_vs_CTRL_venn_pos <- AD_vs_CTRL_res %>%
  filter(log2FoldChange > 0) %>%
  mutate(RegionTissue = paste(Region, Tissue, sep="_")) %>%
  distinct(RegionTissue, Gene)

AD_vs_CTRL_venn_pos <- split(AD_vs_CTRL_venn_pos, AD_vs_CTRL_venn_pos$RegionTissue)
AD_vs_CTRL_venn_pos <- lapply(AD_vs_CTRL_venn_pos, function(x) x %>% pull(Gene))
AD_vs_CTRL_venn_pos <- AD_vs_CTRL_venn_pos[c("AWS_Bulk", "AWS_BV", "PWS_BV", "PWS_Bulk")]

venn_pos <- Venn(AD_vs_CTRL_venn_pos)
data_pos <- process_data(venn_pos, shape_id = "401f")

ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = count, group = id), 
          data = venn_regionedge(data_pos)) +
  # 2. set edge layer
  geom_path(aes(X, Y, group = id), 
          data = venn_setedge(data_pos), 
          show.legend = FALSE) +
  # 3. set label layer
  geom_text(aes(X, Y, label = name), 
               data = venn_setlabel(data_pos)) +
  # 4. region label layer
  geom_label(aes(X, Y, label = count), 
                data = venn_regionlabel(data_pos)) +
  coord_equal() +
  theme_void() +
    scale_fill_distiller(palette = "Reds", direction = 1) +
  labs(fill="# Genes") +
  theme(legend.position="none")

ggsave(paste0(plot_dir, "AD_vs_CTRL_DEGs_Positive.png"),
       width=4.5, height=4.5, units="in", dpi=300)

# Downregulated genes
AD_vs_CTRL_venn_neg <- AD_vs_CTRL_res %>%
  filter(log2FoldChange < 0) %>%
  mutate(RegionTissue = paste(Region, Tissue, sep="_")) %>%
  distinct(RegionTissue, Gene)

AD_vs_CTRL_venn_neg <- split(AD_vs_CTRL_venn_neg, AD_vs_CTRL_venn_neg$RegionTissue)
AD_vs_CTRL_venn_neg <- lapply(AD_vs_CTRL_venn_neg, function(x) x %>% pull(Gene))
AD_vs_CTRL_venn_neg <- AD_vs_CTRL_venn_neg[c("AWS_Bulk", "AWS_BV", "PWS_BV", "PWS_Bulk")]

venn_neg <- Venn(AD_vs_CTRL_venn_neg)
data_neg <- process_data(venn_neg, shape_id = "401f")

ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = count, group = id), 
          data = venn_regionedge(data_neg)) +
  # 2. set edge layer
  geom_path(aes(X, Y, group = id), 
          data = venn_setedge(data_neg), 
          show.legend = FALSE) +
  # 3. set label layer
  geom_text(aes(X, Y, label = name), 
               data = venn_setlabel(data_neg)) +
  # 4. region label layer
  geom_label(aes(X, Y, label = count), 
                data = venn_regionlabel(data_neg)) +
  coord_equal() +
  theme_void() +
    scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(fill="# Genes") +
  theme(legend.position="none")

ggsave(paste0(plot_dir, "AD_vs_CTRL_DEGs_Negative.png"),
       width=4.5, height=4.5, units="in", dpi=300)


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