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

# Load the WMH RNAseq results
load(paste0(study_path, "human_WMH_RNAseq_salmon.Rdata"))

# Read in metadata
metadata <- read.csv(paste0(study_path, "human_sample_info.csv")) %>%
  arrange(Sample_Name) %>%
  mutate(White_Matter_Region = factor(White_Matter_Region, levels=c("AWS", "PWS")))
rownames(metadata) <- colnames(human_WMH_RNAseq_salmon$counts)

# Extract the counts data, in transcripts per million (TPM)
RNAseq_counts <- t(human_WMH_RNAseq_salmon$counts)

# Separate out into bulk AWS and bulk PWS
bulk_AWS_RNAseq_ind <- metadata %>% mutate(index = row_number()) %>% filter(Tissue=="Bulk" & White_Matter_Region=="AWS") %>% pull(index)
bulk_PWS_RNAseq_ind <- metadata %>% mutate(index = row_number()) %>% filter(Tissue=="Bulk" & White_Matter_Region=="PWS") %>% pull(index)
bulk_AWS_RNAseq_counts <- RNAseq_counts[bulk_AWS_RNAseq_ind,(-1)]
bulk_PWS_RNAseq_counts <- RNAseq_counts[bulk_PWS_RNAseq_ind,(-1)]

# Perform principal components analysis (PCA) to identify the first 10 PCs
# Fit the PCAs
bulk_AWS_pca_res <- PCA(bulk_AWS_RNAseq_counts, ncp=10, graph = FALSE, scale.unit = TRUE)
bulk_PWS_pca_res <- PCA(bulk_PWS_RNAseq_counts, ncp=10, graph = FALSE, scale.unit = TRUE)

# Extract the PCs
bulk_AWS_pca_scores <- as.data.frame(bulk_AWS_pca_res$ind$coord) %>% cbind(., metadata[bulk_AWS_RNAseq_ind,])
bulk_PWS_pca_scores <- as.data.frame(bulk_PWS_pca_res$ind$coord) %>% cbind(., metadata[bulk_PWS_RNAseq_ind,])

# Scree plots 
# Extract eigenvalues and percentage of variance
bulk_aws_scree <- as.data.frame(bulk_AWS_pca_res$eig) %>% 
  rownames_to_column(var="Component_Name") %>%
  mutate(pc_number=row_number()) %>%
  filter(pc_number <= 10) %>% 
  janitor::clean_names() %>%
  distinct() %>%
  mutate(line_group = "all", pc_number = as.factor(pc_number)) %>%
  ggplot(data=., mapping=aes(x=pc_number)) +
  geom_line(aes(y=cumulative_percentage_of_variance, group=line_group), linetype=2) +
  geom_bar(aes(y=percentage_of_variance, fill=pc_number), stat="identity") +
  xlab("Principal Component") +
  ylab("% Variance Explained") +
  ggtitle("Bulk AWS Tissue") +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5, face="bold")) +
  scale_y_continuous(expand=c(0,0))
bulk_pws_scree <- as.data.frame(bulk_PWS_pca_res$eig) %>% 
  rownames_to_column(var="Component_Name") %>%
  mutate(pc_number=row_number()) %>%
  filter(pc_number <= 10) %>% 
  janitor::clean_names() %>%
  distinct() %>%
  mutate(line_group = "all", pc_number = as.factor(pc_number)) %>%
  ggplot(data=., mapping=aes(x=pc_number)) +
  geom_line(aes(y=cumulative_percentage_of_variance, group=line_group), linetype=2) +
  geom_bar(aes(y=percentage_of_variance, fill=pc_number), stat="identity") +
  xlab("Principal Component") +
  ylab("% Variance Explained") +
  ggtitle("Bulk PWS Tissue") +
  theme(legend.position="none", 
        plot.title=element_text(hjust=0.5, face="bold")) +
  scale_y_continuous(expand=c(0,0))

bulk_aws_scree / bulk_pws_scree
ggsave(paste0(plot_dir, "PCA_scree_plot.png"),
       width=6, height=6, units="in", dpi=300)

# Plot the first two PC scores per sample and see if AD and CTRL samples are separated in AWS and/or PWS:
# AWS
AWS_biplot <- bulk_AWS_pca_scores %>%
  ggplot(data=., mapping=aes(x=`Dim.1`, y=`Dim.2`, color=Group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  stat_chull(aes(color = Group, fill = Group), 
             alpha = 0.1, geom = "polygon") +
  ylab("PC2") +
  xlab("PC1") +
  ggtitle("AWS Bulk Tissue") +
  scale_color_manual(values=c("CTRL"="#67913c", "AD"="#a864d2")) +
  scale_fill_manual(values=c("CTRL"="#67913c", "AD"="#a864d2")) +
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold"),
        plot.title=element_text(hjust=0.5, face="bold"),
        strip.background = element_blank())

# PWS
PWS_biplot <- bulk_PWS_pca_scores %>%
  ggplot(data=., mapping=aes(x=`Dim.2`, y=`Dim.3`, color=Group)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  stat_chull(aes(color = Group, fill = Group), 
             alpha = 0.1, geom = "polygon") +
  ylab("PC3") +
  xlab("PC2") +
  ggtitle("PWS Bulk Tissue") +
  scale_color_manual(values=c("CTRL"="#67913c", "AD"="#a864d2")) +
  scale_fill_manual(values=c("CTRL"="#67913c", "AD"="#a864d2")) +
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold"),
        plot.title=element_text(hjust=0.5, face="bold"),
        strip.background = element_blank())

AWS_biplot + PWS_biplot + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(paste0(plot_dir, "PCA_biplot.png"),
       width=8, height=4, units="in", dpi=300)

# Compare these biplots with ones colored by other potential confounders:
# AWS
AWS_biplot_sex <- bulk_AWS_pca_scores %>%
  ggplot(data=., mapping=aes(x=`Dim.1`, y=`Dim.2`, color=Sex)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  stat_chull(aes(color = Sex, fill = Sex), 
             alpha = 0.1, geom = "polygon") +
  ylab("PC2") +
  xlab("PC1") +
  ggtitle("AWS Bulk Tissue") +
  scale_color_manual(values=c("M"="#1A59E3", "F"="#BD2752")) +
  scale_fill_manual(values=c("M"="#1A59E3", "F"="#BD2752")) +
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold"),
        plot.title=element_text(hjust=0.5, face="bold"),
        strip.background = element_blank())

# PWS
PWS_biplot_sex <- bulk_PWS_pca_scores %>%
  ggplot(data=., mapping=aes(x=`Dim.2`, y=`Dim.3`, color=Sex)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  stat_chull(aes(color = Sex, fill = Sex), 
             alpha = 0.1, geom = "polygon") +
  ylab("PC3") +
  xlab("PC2") +
  ggtitle("PWS Bulk Tissue") +
  scale_color_manual(values=c("M"="#1A59E3", "F"="#BD2752")) +
  scale_fill_manual(values=c("M"="#1A59E3", "F"="#BD2752")) +
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold"),
        plot.title=element_text(hjust=0.5, face="bold"),
        strip.background = element_blank())

AWS_biplot_sex + PWS_biplot_sex + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(paste0(plot_dir, "PCA_biplot_sex.png"),
       width=8, height=4, units="in", dpi=300)

# As a robustness analysis, we can also look at PC1 vs PC2 in PWS colored by sex
bulk_PWS_pca_scores %>%
  ggplot(data=., mapping=aes(x=`Dim.1`, y=`Dim.2`, color=Sex)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  stat_chull(aes(color = Sex, fill = Sex), 
             alpha = 0.1, geom = "polygon") +
  ylab("PC2") +
  xlab("PC1") +
  ggtitle("PWS Bulk Tissue") +
  scale_color_manual(values=c("M"="#1A59E3", "F"="#BD2752")) +
  scale_fill_manual(values=c("M"="#1A59E3", "F"="#BD2752")) +
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold"),
        plot.title=element_text(hjust=0.5, face="bold"),
        strip.background = element_blank())

# Create another set of biplots colored by age
# AWS
AWS_biplot_age <- bulk_AWS_pca_scores %>%
  rowwise() %>%
  mutate(Age_Label = ifelse(Age==">90", 90, as.numeric(Age))) %>%
  ggplot(data=., mapping=aes(x=`Dim.1`, y=`Dim.2`, color=Age_Label)) +
  scale_color_viridis_c() +
  labs(color="Age") +
  geom_point(size=3) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  ylab("PC2") +
  xlab("PC1") +
  ggtitle("AWS Bulk Tissue") +
  theme(legend.position = "bottom",
        legend.key.width= unit(1, 'cm'),
        strip.text = element_text(face="bold"),
        plot.title=element_text(hjust=0.5, face="bold"),
        strip.background = element_blank())

# PWS
PWS_biplot_age <- bulk_PWS_pca_scores %>%
  rowwise() %>%
  mutate(Age_Label = ifelse(Age==">90", 90, as.numeric(Age))) %>%
  ggplot(data=., mapping=aes(x=`Dim.2`, y=`Dim.3`, color=Age_Label)) +
  scale_color_viridis_c() +
  labs(color="Age") +
  geom_point(size=3) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  ylab("PC3") +
  xlab("PC2") +
  ggtitle("PWS Bulk Tissue") +
  theme(legend.position = "bottom",
        legend.key.width= unit(1, 'cm'),
        strip.text = element_text(face="bold"),
        plot.title=element_text(hjust=0.5, face="bold"),
        strip.background = element_blank())

AWS_biplot_age + PWS_biplot_age + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(paste0(plot_dir, "PCA_biplot_age.png"), width=8, height=4, units="in", dpi=300)