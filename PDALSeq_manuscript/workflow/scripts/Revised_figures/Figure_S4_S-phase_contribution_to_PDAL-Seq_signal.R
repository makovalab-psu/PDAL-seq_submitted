library(tidyverse)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)

df = readr::read_csv("results/25_Sphase_contribution/Technical_replicate_correlation.csv")

Sample_levels = unique(df$Sample_1)
df$Window_size = as.character(df$Window_size)
Window_levels = c("10000", "1e+05", "1e+06")
Window_labels = c("10000", "100000", "1000000")
df$Window_size = factor(df$Window_size, levels = Window_levels, labels = Window_labels)   

df$Sample_1 = factor(df$Sample_1, levels = Sample_levels)
df$Sample_2 = factor(df$Sample_2, levels = Sample_levels)

P = ggplot(df, aes(x = Sample_1, y = Sample_2, fill = Spearman)) +
	facet_wrap(~Window_size) +	
	geom_tile() +
	theme_classic() +
	theme(axis.text.x = element_text(angle =90),
	      axis.text.y = element_blank(),
	      axis.title = element_blank()) +
	scale_fill_viridis(limits = c(0,1), option = "turbo") +
	coord_fixed()
	
ggsave("revised_plots/Figure_S4_S-phase_contribution_to_PDAL-Seq_signal.svg", P, width = 5, height = 10, scale = 3)
