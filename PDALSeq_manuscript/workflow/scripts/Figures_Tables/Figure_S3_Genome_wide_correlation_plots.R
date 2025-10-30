####Figure 1D Correlation####

library(tidyverse)
library(viridis)
library(cowplot)

v.margins = c(0.1,0.1,0.1,0.1)

nonb = c("all", "APR", "IR", "CRU", "MR", "TRI", "STR", "DR", "SLS",  "GQ", "Z")

l.df.all = {}

####99% Winsorization####

for (k in 1:length(nonb)){
	df.windows = read.csv("results/04_Genome_wide_comparison_to_nonB_motifs/nonB_coverage_versus_read.csv") %>% filter(nonB == nonb[k])
	df.windows$Experiment = paste(df.windows$Species, df.windows$Bio_sample, df.windows$Technique, sep = "_")
	
	l.df = df.windows %>% group_by(Experiment) %>% group_split()
	l.df2 ={}
	
	for (i in 1:length(l.df)){
    		l.df2[[i]] = data.frame(l.df[[i]])
    		Read.Q99 = quantile(l.df[[i]]$reads, 0.99)
    		Motif.Q99 = quantile(l.df[[i]]$nonb.count, 0.99)
    		l.df[[i]]$reads[which(l.df[[i]]$reads > Read.Q99)] = as.integer(Read.Q99)
    		l.df[[i]]$nonb.count[which(l.df[[i]]$nonb.count > Motif.Q99)] = as.integer(Motif.Q99)
    		l.df[[i]] = data.frame(l.df[[i]])
    		l.df2[[i]]$reads = l.df[[i]]$reads/Read.Q99
    		l.df2[[i]]$nonb.count = l.df[[i]]$nonb.count/Motif.Q99
    		l.df2[[i]]$interval = 1:nrow(l.df2[[i]]) 
	}
	
	l.df.all[[k]] = bind_rows(l.df2)
}

df.windows = bind_rows(l.df.all)

df.windows$nonB = factor(df.windows$nonB, levels = nonb)

Group = rep("none", nrow(df.windows))
Group[which(grepl("chrX", df.windows$chr))] = "chrX"
Group[which(grepl("chrY", df.windows$chr))] = "chrY"
df.windows$Group = Group

Experiments = unique(df.windows$Experiment)

for (x in Experiments){
	df = df.windows %>% filter(Experiment == x)
	P = ggplot(df, aes(x = nonb.count, y = reads, color = Group)) +
		geom_point(alpha = 0.1) +
		theme_classic() +
	facet_wrap(~nonB, nrow = 1) +
	scale_color_manual(values = c("#3B528BFF", "#AADC32FF", "black")) +
	coord_fixed() +
	xlab("non-B DNA motifs") +
	ylab("PDAL-seq reads") +
	ggtitle(x) +
	theme(axis.text.x = element_text(angle = 90))
	ggsave(paste("plots/Figure_S3_Genome_wide_correlation_plots/", x, ".png", sep = ""), P, width = 8, height = 1.1, units = "in", scale = 2.5) 	
}

