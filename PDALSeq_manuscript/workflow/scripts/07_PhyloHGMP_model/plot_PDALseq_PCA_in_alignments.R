library(tidyverse)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)
#args = c("data/07_PhyloHGMP_model/alignment_clustering/hominini/PDALseq_PCA_in_alignment.csv", "data/07_PhyloHGMP_model/alignment_clustering/hominini/PDALseq_PCA_plot.svg")

df.PCA = read.csv(args[1])

df.PCA$Species = factor(df.PCA$Species, levels = unique(df.PCA$Species))

P.PCA = ggplot(df.PCA, aes(x = PC1, y = PC2, color = Species, shape = Technique)) +
	geom_point() +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(df.PCA$Species))+1, option = "turbo"))
	
ggsave(args[2], P.PCA, width = 5.5, height = 4, units = "in", scale = 1)

