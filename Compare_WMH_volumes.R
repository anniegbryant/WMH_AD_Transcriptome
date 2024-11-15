library(tidyverse)
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

# Read in metadata
metadata <- read.csv(paste0(study_path, "human_sample_info.csv")) %>%
  arrange(Sample_Name) %>%
  mutate(White_Matter_Region = factor(White_Matter_Region, levels=c("AWS", "PWS")))
rownames(metadata) <- colnames(human_WMH_RNAseq_salmon$counts)

# Plot the WMH volumes as a bar chart
human_sample_info %>% 
    mutate(WMH_Group = ifelse(WMH_mm3 > 20000, "High", "Low")) %>% 
    distinct(ADRC, Group, WMH_mm3, WMH_Group) %>%
    filter(!is.na(WMH_mm3)) %>%
    mutate(ADRC = fct_reorder(ADRC, WMH_mm3)) %>%
    ggplot(data=., mapping=aes(x=ADRC, y=WMH_mm3, fill=WMH_Group)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("Low"="#efa3b0", "High"="#76170b")) +
    xlab("Donors") +
    ylab("WMH volume (mm3)") +
    labs(fill="WMH Group") +
    scale_y_continuous(expand=c(0,0)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) 

ggsave(paste0(plot_dir, "WMH_volumes.svg"), width=6, height=3, units="in", dpi=300)