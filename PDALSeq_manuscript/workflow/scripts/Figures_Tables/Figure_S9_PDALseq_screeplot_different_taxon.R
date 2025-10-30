library(tidyverse)
library(cowplot)

args = commandArgs(trailingOnly = TRUE)

taxon = c("hominoidae", "hominidae", "homininae", "hominini")

####Kmeans####

l.df = {}
for (i in 1:length(taxon)){
	file_path = paste("data/07_PhyloHGMP_model/alignment_clustering/", taxon[i], "/PDALseq_kmeans_scree_plot_in_alignment.csv", sep = "")
	df.i = read.csv(file_path)
	df.i$taxon = taxon[i]
	l.df[[i]] = df.i
}

df = bind_rows(l.df) 

df$taxon = factor(df$taxon, levels = taxon)

P.Kmeans = ggplot(df, aes(x = n_clusters, y = wss)) +
	geom_point() +
	geom_line() +
	theme_classic() +
	facet_wrap(~taxon, scales = "free") +
	xlab("cluster number (k)")+
	ylab("Within-cluster-Sum\nof Square Error (WSS)")

####HMM BIC####

l.df = {}
for (i in 1:length(taxon)){
	file_path = paste("data/07_PhyloHGMP_model/Multivariate_guassian_HMM/", taxon[i], "/BIC.csv", sep = "")
	df.i = read.csv(file_path)
	df.i$taxon = taxon[i]
	l.df[[i]] = df.i
}

df = bind_rows(l.df) 

df$taxon = factor(df$taxon, levels = taxon)

P.BIC = ggplot(df, aes(x =States, y = BIC)) +
	geom_point() +
	geom_line() +
	theme_classic() +
	facet_wrap(~taxon, scales = "free") +
	xlab("Hidden States (n)")+
	ylab("Bayesian Information Criterion (BIC)")

####Compile and save plots####

P.out = plot_grid(P.Kmeans, P.BIC, ncol = 1, labels = c("A", "B"))

ggsave(args[1], P.out, height = 7, width = 6, scale = 1, limitsize = FALSE)
