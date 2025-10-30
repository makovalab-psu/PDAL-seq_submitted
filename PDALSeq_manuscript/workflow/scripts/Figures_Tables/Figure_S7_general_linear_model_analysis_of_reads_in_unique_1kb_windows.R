library(dplyr)
library(viridis)
library(ggbeeswarm)
library(ggpubr)

args = commandArgs(trailingOnly = TRUE)

df.glm = read.csv("results/16_lm_analysis_of_nonbmotifs/glm_unique_regions_1000_nucleotides.csv")  
df.lm = df.lm %>% filter(Fit == "Deduplicated parameters")

df.lm$Species = factor(df.lm$Species, levels = c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")) 
df.lm$X = factor(df.lm$X, levels = c("(Intercept)", "length.N", "GC.content", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ"))

P.lm = ggplot(data = df.lm, mapping = aes(x = X, y = Estimate)) +
	facet_wrap(~Technique, ncol = 1) +
	geom_bar(stat = "summary", fun = mean, fill = "grey") +
	stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black") +
	geom_beeswarm(mapping = aes(color = Species)) +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(df.lm$Species)), option = "turbo")) +
	geom_hline(yintercept = 0) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) +
	      ylab("Reads per motif\nper million reads")

ggsave(args[1], P.lm, width = 5, height = 4, units = "in", scale = 2)

df.lm = read.csv("results/05_nonB_motif_cluster_analysis/log10_linear_models.csv") %>% filter(Fit == "Deduplicated parameters")

df.lm$Species = factor(df.lm$Species, levels = c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")) 
df.lm$X = factor(df.lm$X, levels = c("(Intercept)", "length.N", "GC.content", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ"))

P.log10_lm = ggplot(data = df.lm, mapping = aes(x = X, y = Estimate)) +
	facet_wrap(~Technique, ncol = 1) +
	geom_bar(stat = "summary", fun = mean, fill = "grey") +
	stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black") +
	geom_beeswarm(mapping = aes(color = Species)) +
	theme_classic() +
	scale_color_manual(values = viridis(length(unique(df.lm$Species)), option = "turbo")) +
	geom_hline(yintercept = 0) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)) +
	      ylab("Reads per motif\nper million reads")

ggsave("plots/Figure_S7_linear_model_analysis_of_nonB_motif_clusters/Figure_S7_log10.svg", P.log10_lm, width = 5, height = 4, units = "in", scale = 2)
