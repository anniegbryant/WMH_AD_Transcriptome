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
load(paste0(study, "Bulk_vs_BV_dds_res.Rdata.Rdata"))

# Load the metadata
metadata <- read.csv(paste0(study, "human_sample_info.csv")) %>%
  arrange(Sample_Name) %>%
  mutate(Group = factor(Group, levels=c("AD", "CTRL")))

# Filter DESeq2 results to BH-adjusted p-value < 0.05 and |lgo2FC| >= 0.5
AD_AWS_dds_df <- as.data.frame(AD_AWS_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Group = "AD",
         Region = "AWS")

CTRL_AWS_dds_df <- as.data.frame(CTRL_AWS_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Group = "CTRL",
         Region = "AWS")

AD_PWS_dds_df <- as.data.frame(AD_PWS_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Group = "AD",
         Region = "PWS")

CTRL_PWS_dds_df <- as.data.frame(CTRL_PWS_dds_res) %>%
  filter(padj < 0.05,
         abs(log2FoldChange) >= 0.5) %>%
  rownames_to_column(var="Gene") %>%
  mutate(Group = "CTRL",
         Region = "PWS")

# Combine the results
Bulk_vs_BV_res <- do.call(plyr::rbind.fill, list(AD_AWS_dds_df,
                                                 CTRL_AWS_dds_df,
                                                 AD_PWS_dds_df,
                                                 CTRL_PWS_dds_df))

# Count the number of significant DEGs in each region and diagnostic group
Bulk_vs_BV_res %>% 
  mutate(Direction = ifelse(log2FoldChange>0, "Up", "Down")) %>% 
  group_by(Group, Region, Direction) %>%
  summarise(count = n())

# Euler diagrams
#### Downregulated (i.e., BV > Bulk)
BV_vs_bulk_venn_pos <- Bulk_vs_BV_res %>%
  filter(log2FoldChange > 0) %>%
  mutate(RegionDx = paste(Region, Group, sep="_")) %>%
  distinct(RegionDx, Gene)

BV_vs_bulk_venn_pos <- split(BV_vs_bulk_venn_pos, BV_vs_bulk_venn_pos$RegionDx)
BV_vs_bulk_venn_pos <- lapply(BV_vs_bulk_venn_pos, function(x) x %>% pull(Gene))
BV_vs_bulk_venn_pos <- BV_vs_bulk_venn_pos[c("AWS_AD", "AWS_CTRL", "PWS_CTRL", "PWS_AD")]

venn_pos <- Venn(BV_vs_bulk_venn_pos)
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

ggsave(paste0(plot_dir, "BV_vs_Bulk_DEGs_Positive.svg"),
       width=4.5, height=4.5, units="in", dpi=300)

#### Downregulated (i.e., Bulk > BV)
BV_vs_bulk_venn_neg <- Bulk_vs_BV_res %>%
  filter(log2FoldChange < 0) %>%
  mutate(RegionDx = paste(Region, Group, sep="_")) %>%
  distinct(RegionDx, Gene)

BV_vs_bulk_venn_neg <- split(BV_vs_bulk_venn_neg, BV_vs_bulk_venn_neg$RegionDx)
BV_vs_bulk_venn_neg <- lapply(BV_vs_bulk_venn_neg, function(x) x %>% pull(Gene))
BV_vs_bulk_venn_neg <- BV_vs_bulk_venn_neg[c("AWS_AD", "AWS_CTRL", "PWS_CTRL", "PWS_AD")]

venn_neg <- Venn(BV_vs_bulk_venn_neg)
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

ggsave(paste0(plot_dir, "BV_vs_Bulk_DEGs_negative.svg"),
       width=4.5, height=4.5, units="in", dpi=300)

### Expression-Weighted Celltype Enrichment (EWCE) Analysis
# Prepare the CTD dataset from mouse cortex+hippocampus from the Karolinska Institute, provided with EWCE package
cortex_mrna <- ewceData::cortex_mrna()

# Fix bad gene symbols
cortex_mrna$exp = fix_bad_mgi_symbols(cortex_mrna$exp)

# Normalize for different number of reads found across each cell
cortex_mrna$exp_scT_normed <- EWCE::sct_normalize(cortex_mrna$exp) 

# Generate cell type data for just the cortex/hippocampus data  
exp_CortexOnly_DROPPED <- EWCE::drop_uninformative_genes(
  exp = cortex_mrna$exp, 
  input_species = "mouse",
  output_species = "human",
  level2annot = cortex_mrna$annot$level2class) 

# Generate CellTypeDataset:
annotLevels <- list(level1class=cortex_mrna$annot$level1class,
                    level2class=cortex_mrna$annot$level2class)

fNames_CortexOnly <- EWCE::generate_celltype_data(
  exp = exp_CortexOnly_DROPPED,
  annotLevels = annotLevels,
  groupName = "kiCortexOnly") 
ctd_CortexOnly <- EWCE::load_rdata(fNames_CortexOnly)

# Look at the different cell types present:
cortex_mRNA_annot <- cortex_mrna$annot
cortex_mRNA_annot %>%
  distinct(level1class, level2class)

# Perform enrichment analysis with the top 500 genes that are upregulated across blood vessels throughout the brain in AD and controls:
# Our 'hits' are the BV genes upregulated in both AWS + PWS in CTRL participants
top500_BV_genes <- Bulk_vs_BV_res %>% 
  group_by(Gene) %>% 
  filter(log2FoldChange > 0 & Group == "CTRL") %>%
  filter(n_distinct(Region) == 2 & all(c("AWS", "PWS") %in% Region)) %>% 
  summarise(mean_log2FC = mean(log2FoldChange, na.rm=T)) %>% 
  arrange(desc(mean_log2FC)) %>%
  slice_max(order_by=mean_log2FC, n=500) %>%
  pull(Gene)

# Repeat 10000 bootstrap permutations
reps <- 10000

# Bootstrap significance test, no control for transcript length and GC content 
# Level 1 annotation results
top500_BV_genes_results_annotlevel_1 <- EWCE::bootstrap_enrichment_test(sct_data = ctd_CortexOnly,
                                                        sctSpecies = "mouse",
                                                        sctSpecies_origin = "mouse",
                                                        genelistSpecies = "human",
                                                        output_species = "human",
                                                        hits = top500_BV_genes, 
                                                        mtc_method = "BH",
                                                        reps = reps,
                                                        no_cores = 2,
                                                        annotLevel = 1)

# Level 2 annotation results
top500_BV_genes_results_annotlevel_2 <- EWCE::bootstrap_enrichment_test(sct_data = ctd_CortexOnly,
                                                        sctSpecies = "mouse",
                                                        sctSpecies_origin = "mouse",
                                                        genelistSpecies = "human",
                                                        output_species = "human",
                                                        hits = top500_BV_genes, 
                                                        mtc_method = "BH",
                                                        reps = reps,
                                                        no_cores = 2,
                                                        annotLevel = 2)

# Let's visualize the SD from the mean in the significant level 2 annotation results
top500_BV_genes_results_annotlevel_2$results %>% 
  filter(q<0.05) %>%
  rename("CellType" = "level2class") %>%
  left_join(., cortex_mRNA_annot %>% distinct(level1class, level2class)) %>%
  mutate(level2class = fct_reorder(level2class, sd_from_mean, .fun="max", .desc=T),
         sd_label = round(sd_from_mean, 1)) %>%
  ggplot(data=., mapping=aes(x=level2class, y=sd_from_mean)) +
  xlab("Level 2 Cell Type") +
  labs(fill="") +
  scale_fill_manual(values=c("endothelial-mural" = "#EB806E",
                             "microglia" = "#4B9698")) +
  scale_color_manual(values=c("endothelial-mural" = "#EB806E",
                             "microglia" = "#4B9698")) +
  ylab("SD from Mean Across Cells") +
  geom_text(aes(label=sd_label, color=level1class), 
            nudge_y = 1.5, show.legend = FALSE) +
  geom_bar(stat="identity", aes(fill=level1class)) +
  scale_y_continuous(expand=c(0,0,0.05,0)) +
  theme(legend.position="bottom")
ggsave(paste0(plot_dir, "EWCE_enrichment.svg"), width=6, height=3.25, units="in", dpi=300)