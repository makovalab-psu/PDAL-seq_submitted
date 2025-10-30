library(tidyverse)
library(viridis)

args = commandArgs(trailingOnly = TRUE)

taxon = c("hominoidae", "hominidae", "homininae", "hominini")

l.df = {}
for (i in 1:length(taxon)){
	file_path = paste("data/07_PhyloHGMP_model/alignment_clustering/", taxon[i], "/PDALseq_PCA_in_alignment.csv", sep = "")
	df.i = read.csv(file_path)
	df.i$taxon = taxon[i]
	l.df[[i]] = df.i
}

df = bind_rows(l.df) 

SPECIES = c("Hsapien", "Ptroglodytes", "Ppaniscus", "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")
SPECIES.label = c("H. sapiens", "P. troglodytes", "P. paniscus", "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")
df$Species = factor(df$Species, levels = SPECIES.label) 
df$taxon = factor(df$taxon, levels = taxon)

P.PCA = ggplot(df, aes(x = PC1, y = PC2, color = Species, shape = Technique)) +
	geom_point() +
	theme_classic() +
	facet_wrap(~taxon) +
	scale_color_manual(values = viridis(length(unique(df$Species))+1, option = "turbo"))


ggsave(args[1], P.PCA, height = 5, width = 6, scale = 1, limitsize = FALSE)
