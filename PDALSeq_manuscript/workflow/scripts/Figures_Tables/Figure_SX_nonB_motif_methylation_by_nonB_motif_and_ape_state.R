
library(dplyr)
library(viridis)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)


nonb_levels = c("all", "APR", "IR", "CRU", "MR", "TRI", "Z", "STR", "DR", "SLS", "GQ")
State.levels = c("4", "7", "3", "1", "6", "0", "5", "2")


df.5mC = read.csv("results/10_Methylation_analysis/Methylation_by_state_nonB_motif_and_species.csv")

df.5mC$State = factor(as.character(df.5mC$State), levels = State.levels)
df.5mC$Species = factor(df.5mC$Species, levels = rev(c("Hsapien","Ptroglodytes", "Ppaniscus",    "Ggorilla", "Ppygmaeus", "Pabelii", "Ssyndactylus")), labels = rev(c("H. sapiens","P. troglodytes", "P. paniscus",    "G. gorilla", "P. pygmaeus", "P. abelii", "S. syndactylus")))
df.5mC$nonb = factor(df.5mC$nonb, levels = nonb_levels)
Min_max = c(0, max(df.5mC$Percent))
l.df = df.5mC %>% group_by(State) %>% group_split()

l.P = {}
for (i in 1:length(l.df)){
	df.i = l.df[[i]]
	l.P[[i]] = ggplot(df.i, aes(x = Range, y = Species, fill = Percent)) +
		ggtitle(paste("Ape state =", df.i$State[1])) + 
		facet_wrap(~nonb, ncol = 1, strip.position = "left") +
		geom_tile() +
		scale_fill_gradient(low = "white", high = "black", limits = Min_max) +
		theme_classic() +
		scale_y_discrete(position = "right") +
		theme(axis.text.x = element_text(angle = 0),
	      		legend.position = "right") +
		xlab("5mC rate")
}

P = plot_grid(l.P[[1]],l.P[[2]],l.P[[3]],l.P[[4]],l.P[[5]],l.P[[6]],l.P[[7]],l.P[[8]], ncol = 1)

ggsave("plots/Figure_X_CpG_methylation_by_PDALseq_state/Figure_X.svg", P, width = 6.5, height = 50, units = "in", limitsize = F)
