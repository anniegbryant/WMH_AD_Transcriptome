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

# Load AD vs. Control DESeq2 results
load(paste0(study_path, "AD_vs_CTRL_dds_res.Rdata"))

metadata <- read.csv(paste0(study_path, "human_sample_info.csv")) %>%
  arrange(Sample_Name) %>%
  mutate(Group = factor(Group, levels=c("AD", "CTRL")))

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

AD_vs_CTRL_BV_res <- do.call(plyr::rbind.fill, list(AWS_bv_dds_df,
                                                    PWS_bv_dds_df))


# How many unique AD-associated DEGs are there in the blood vessels?
AD_vs_CTRL_BV_res %>% distinct(Gene) %>% nrow()

# Load in extended data table 7-2 from Bryant et al. (2023) for comparison:
bryant_2023_endothelial_DEGs <- read.csv("Bryant_2023_Extended_Data_Table_7_2.csv", check.names = FALSE) %>% janitor::clean_names()

# Venn Diagram for upregulated genes
# Upregulated genes
AWS_BV_upregulated_genes <- AD_vs_CTRL_BV_res %>% filter(log2FoldChange>0, Region=="AWS") %>% distinct(Gene) %>% pull(Gene)
PWS_BV_upregulated_genes <- AD_vs_CTRL_BV_res %>% filter(log2FoldChange>0, Region=="PWS") %>% distinct(Gene) %>% pull(Gene)
Bryant_2023_upregulated_genes <- bryant_2023_endothelial_DEGs %>% filter(average_log2_fc>0) %>% distinct(gene) %>% pull(gene)
AD_vs_CTRL_venn_pos <- list( "Bryant_2023" = Bryant_2023_upregulated_genes, "AWS_BV" = AWS_BV_upregulated_genes, "PWS_BV" = PWS_BV_upregulated_genes)

venn_pos <- Venn(AD_vs_CTRL_venn_pos)
data_pos <- process_data(venn_pos, shape_id = "301")

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
ggsave(paste0(plot_dir, "BV_vs_Bryant2023_upregulated.svg"),
       width=4.5, height=4.5, units="in", dpi=300)

# Venn Diagram for downregulated genes
AWS_BV_downregulated_genes <- AD_vs_CTRL_BV_res %>% filter(log2FoldChange<0, Region=="AWS") %>% distinct(Gene) %>% pull(Gene)
PWS_BV_downregulated_genes <- AD_vs_CTRL_BV_res %>% filter(log2FoldChange<0, Region=="PWS") %>% distinct(Gene) %>% pull(Gene)
Bryant_2023_downregulated_genes <- bryant_2023_endothelial_DEGs %>% filter(average_log2_fc<0) %>% distinct(gene) %>% pull(gene)
AD_vs_CTRL_venn_neg <- list( "Bryant_2023" = Bryant_2023_downregulated_genes, 
                             "AWS_BV" = AWS_BV_downregulated_genes, 
                             "PWS_BV" = PWS_BV_downregulated_genes)

venn_neg <- Venn(AD_vs_CTRL_venn_neg)
data_neg <- process_data(venn_neg, shape_id = "301")

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
ggsave(paste0(plot_dir, "BV_vs_Bryant2023_downregulated.svg"),
       width=4.5, height=4.5, units="in", dpi=300)

# What are the 10 upregulated genes that overlap between the current study and with Bryant 2023?
shared_upregulation <- Reduce(intersect, list(AWS_BV_upregulated_genes,
                                              PWS_BV_upregulated_genes,
                                              Bryant_2023_upregulated_genes))
shared_upregulation

# What is the one downregulated gene that overlaps between the current study and with Bryant 2023?
shared_downregulation <- Reduce(intersect, list(AWS_BV_downregulated_genes,
                                              PWS_BV_downregulated_genes,
                                              Bryant_2023_downregulated_genes))
shared_downregulation

# Read in upregulated endothelial data from Mitroi et al. (2022)
mitroi_vad_genes <- read.csv("Mitroi_2022_VaD_vs_NC.csv", check.names = FALSE) %>% janitor::clean_names()
mitroi_vadadj_genes <-  read.csv("Mitroi_2022_VaDadj_vs_NC.csv", check.names = FALSE) %>% janitor::clean_names()

# Recreate 5-circle Euler diagram that now includes Mitroi VaD and VaDadj genes
mitroi_vad_upregulated_genes <- mitroi_vad_genes %>% filter(avg_log2fc>0) %>% distinct(gene) %>% pull(gene)
mitroi_vadadj_upregulated_genes <- mitroi_vadadj_genes %>% filter(avg_log2fc>0) %>% distinct(gene) %>% pull(gene)
AD_vs_CTRL_venn_pos <- list("Mitroi_VaD" = mitroi_vad_upregulated_genes, "Mitroi_VaDadj" = mitroi_vadadj_upregulated_genes, "Bryant_2023" = Bryant_2023_upregulated_genes, "AWS_BV" = AWS_BV_upregulated_genes, "PWS_BV" = PWS_BV_upregulated_genes)

venn_pos <- Venn(AD_vs_CTRL_venn_pos)
data_pos <- process_data(venn_pos, shape_id = "502")

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
  geom_text(aes(X, Y, label = count), 
                data = venn_regionlabel(data_pos)) +
  coord_equal() +
  theme_void() +
    scale_fill_distiller(palette = "Reds", direction = 1) +
  labs(fill="# Genes") +
  theme(legend.position="none")
ggsave(paste0(plot_dir, "BV_vs_Bryant2023_vs_Mitroi2022.svg"),
       width=4.5, height=4.5, units="in", dpi=300)

# What are the three upregulated genes that overlap across all five gene sets?
shared_upregulation_all5 <- Reduce(intersect, AD_vs_CTRL_venn_pos)
shared_upregulation_all5