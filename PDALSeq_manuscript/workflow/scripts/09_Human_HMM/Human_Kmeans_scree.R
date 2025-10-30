library(tidyverse)
library(cowplot)
library(viridis)

df = read.csv("data/09_Human_HMM/02_Normalize_human_window_coverage/Normalized_human_window_coverage.csv")

####Correlation####

PDALseq = colnames(df)
Experiments = PDALseq
df.cor = data.frame(expand.grid(PDALseq, PDALseq))
colnames(df.cor) = c("Sample_1", "Sample_2")
pearson = c()
spearman = c()
for (i in 1:nrow(df.cor)){
	print(i)
	x = which(Experiments == df.cor$Sample_1[i])
	y = which(Experiments == df.cor$Sample_2[i])
	t.pearson = cor.test(df[[x]], df[[y]], method = "pearson")
	t.spearman = cor.test(df[[x]], df[[y]], method = "spearman")
	pearson[i] = t.pearson$estimate 
	spearman[i] = t.spearman$estimate
	print(paste(x, y, spearman[i]))
}

df.cor$pearson = pearson
df.cor$spearman = spearman

write.csv(df.cor, "data/09_Human_HMM/03_Human_Kmeans_scree/Correlation_data.csv", row.names = F)

df.cor$Sample_1 = factor(df.cor$Sample_1, levels = PDALseq)
df.cor$Sample_2 = factor(df.cor$Sample_2, levels = PDALseq)

P.pearson = ggplot(df.cor, aes(x = Sample_1, y = Sample_2, fill = pearson)) +
	geom_tile() +
	scale_fill_viridis(option = "turbo", limits = c(-1, 1)) +
	theme_classic() +
	coord_fixed() +
	theme(axis.text.y = element_blank(),
	      axis.text.x = element_text(angle = 90))
	      
P.spearman = ggplot(df.cor, aes(x = Sample_1, y = Sample_2, fill = spearman)) +
	geom_tile() +
	scale_fill_viridis(option = "turbo", limits = c(-1, 1)) +
	theme_classic() +
	coord_fixed() +
	theme(axis.text.y = element_blank(),
	      axis.text.x = element_text(angle = 90))
  
P.cor = plot_grid(P.pearson, P.spearman, nrow = 1)

ggsave("data/09_Human_HMM/03_Human_Kmeans_scree/Correlation_data.png", P.cor, width = 5, height = 4, scale = 2)

####Scree plot####

# Decide how many clusters to look at
n_clusters = seq(2, 20, by = 2)

# Initialize total within sum of squares error: wss
wss = rep(NA, length(n_clusters))

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:length(n_clusters)) {
  print(i)
  # Fit the model: km.out
  km.out = kmeans(df, centers = n_clusters[i], nstart = 20)
  # Save the within cluster sum of squares
  wss[i] = km.out$tot.withinss
  print(wss[i])
}

# Produce a scree plot
df.scree = data.frame(n_clusters, wss)

#Write the data
write.csv(df.scree, "data/09_Human_HMM/03_Human_Kmeans_scree/Kmeans_scree_data.csv", row.names = F)

#Make the plot

P.SCREE = ggplot(df.scree, aes(x = n_clusters, y = wss)) +
	geom_point() +
	geom_line() +
	theme_classic() +
	xlab("number of clusters (k)") +
	ylab("Total within sum of squares")
	
ggsave("data/09_Human_HMM/03_Human_Kmeans_scree/Kmeans_scree_plot.svg", P.SCREE, width = 5, height = 4, units = "in", scale = 1)

