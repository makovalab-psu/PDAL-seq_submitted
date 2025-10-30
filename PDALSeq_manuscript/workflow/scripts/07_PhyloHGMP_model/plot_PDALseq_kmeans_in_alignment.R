library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
#args = c("data/07_PhyloHGMP_model/alignment_clustering/hominidae/PDALseq_kmeans_scree_plot_in_alignment.csv", "data/07_PhyloHGMP_model/alignment_clustering/hominidae/PDALseq_kmeans_scree_plot_in_alignment.png", "hominidae")

df = read.csv(args[1])

P = ggplot(df, aes(x = n_clusters, y = wss)) +
	geom_point() +
	geom_line() +
	theme_classic() +
	xlab("number of clusters (k)") +
	ylab("Total within sum of squares")
	
ggsave(args[2], P, width = 5, height = 4, units = "in", scale = 1)

