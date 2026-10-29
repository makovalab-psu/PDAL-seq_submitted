library(dplyr)
library(viridis)
library(ggbeeswarm)
library(ggpubr)

args = commandArgs(trailingOnly = TRUE)

df.glm = read.csv("results/21_add_K562_to_analysis/glm_unique_regions_1000_nucleotides.csv")  
df.lm = df.glm %>% filter(Fit %in% c("dedup_GC", "dedup_no_GC"))

df.lm$Species = factor(df.lm$Species, levels = c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus"), labels = c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")) 
df.lm$X = factor(df.lm$Coef, levels = c("(Intercept)", "length.N", "GC_content", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS",  "GQ"))

P.lm = ggplot(data = df.lm, mapping = aes(x = X, y = Estimate, fill = Fit)) +
	facet_wrap(~Technique, nrow = 1) +
	geom_bar(stat = "summary", fun = mean, position = "dodge") +
	stat_summary(fun.data = mean_sd, geom = "errorbar", color = "black", position = "dodge") +
	theme_classic() +
	scale_fill_manual(values = viridis(4)[c(1,3)]) +
	geom_hline(yintercept = 0) +
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1),
	      legend.position = "bottom") +
	      ylab("Reads per motif\nper million reads")

ggsave("revised_plots/Figure_S6_glm_analysis_of_nonB_motif_clusters.svg", P.lm, width = 6.5, height = 4, units = "in")
