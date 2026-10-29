library(tidyverse)
library(cowplot)
library(viridis)

df = read.csv("data/26_Ape_PCA_and_Human_HMM_with_K562/05_normalized_coverage/Normalized_human_window_coverage_with_K562.csv")

PDALseq = colnames(df)

df.cor = data.frame(expand.grid(PDALseq, PDALseq))
colnames(df.cor) = c("Sample_1", "Sample_2")
pearson  = c()
spearman = c()

for (i in 1:nrow(df.cor)){
	print(i)
	x = which(PDALseq == df.cor$Sample_1[i])
	y = which(PDALseq == df.cor$Sample_2[i])
	t.pearson  = cor.test(df[[x]], df[[y]], method = "pearson")
	t.spearman = cor.test(df[[x]], df[[y]], method = "spearman")
	pearson[i]  = t.pearson$estimate
	spearman[i] = t.spearman$estimate
}

df.cor$pearson  = pearson
df.cor$spearman = spearman

dir.create("results/26_Ape_PCA_and_Human_HMM_with_K562", recursive = TRUE, showWarnings = FALSE)
write.csv(df.cor, "results/26_Ape_PCA_and_Human_HMM_with_K562/Correlation_data.csv", row.names = FALSE)

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
ggsave("results/26_Ape_PCA_and_Human_HMM_with_K562/Correlation_data.png",
       P.cor, width = 5, height = 4, scale = 2)

n_clusters = seq(2, 20, by = 2)
wss = rep(NA, length(n_clusters))
set.seed(123)

for (i in seq_along(n_clusters)){
	print(i)
	km.out = kmeans(df, centers = n_clusters[i], nstart = 20)
	wss[i] = km.out$tot.withinss
	print(wss[i])
}

df.scree = data.frame(n_clusters, wss)
write.csv(df.scree, "results/26_Ape_PCA_and_Human_HMM_with_K562/Kmeans_scree_data.csv", row.names = FALSE)

P.SCREE = ggplot(df.scree, aes(x = n_clusters, y = wss)) +
	geom_point() +
	geom_line() +
	theme_classic() +
	xlab("number of clusters (k)") +
	ylab("Total within sum of squares")

ggsave("results/26_Ape_PCA_and_Human_HMM_with_K562/Kmeans_scree_plot.svg",
       P.SCREE, width = 5, height = 4, units = "in")
