library(tidyverse)
library(viridis)

df.unfiltered = read.delim("results/05_nonB_motif_cluster_analysis/Cluster_count_vs_cluster_size.txt", header = F, sep = "")
df.filtered = read.delim("results/05_nonB_motif_cluster_analysis/Filtered_cluster_count_vs_cluster_size.txt", header = F, sep = "")

colnames(df.unfiltered) = c("Count", "File")
colnames(df.filtered) = c("Count", "File")

df.unfiltered$Filtered = "Before filtering"
df.filtered$Filtered = "After filtering"

df = bind_rows(df.unfiltered, df.filtered) %>% filter(File != "total")

l.file = strsplit(df$File, split = "/")

Species = c()
Cluster.distance = c()

for (i in 1:nrow(df)){
	Species[i] = gsub(".bed", "", l.file[[i]][5])
	Cluster.distance[i] = l.file[[i]][4]
}

df$Species = Species
df$Cluster.distance = as.numeric(Cluster.distance)

df$Species = factor(df$Species, levels = c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")) 

P = ggplot(df, aes(x = Cluster.distance, y = Count, color = Species)) +
	geom_point() +
	geom_line() +
	theme_classic() +
	scale_x_continuous(trans = "log10", breaks = floor(10^c(1:8/2))) +
	scale_y_continuous(trans = "log10", breaks = floor(10^c(1:16/2))) +
	facet_wrap(~Filtered) +
	theme(legend.position = "top") +
	scale_color_manual(values = viridis(length(unique(df$Species)), option = "turbo")) +
	xlab("Cluster distance (bp) - log10 tranformed") +
	ylab("non-B motif clusters (count)\nlog10 tranformed")
	
ggsave("plots/Figure_S6_Number_of_non-B_DNA_cluster_motifs_versus_clustering_distance/Figure_S6.svg", P, width = 6, height = 4, units = "in", scale = 2)
