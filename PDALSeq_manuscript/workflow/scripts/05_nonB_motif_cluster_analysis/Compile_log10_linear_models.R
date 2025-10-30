library(tidyverse)

df.index = read.delim("PDALseq_samples.txt")

Experiment = unique(df.index$Experiment)

l.Experiment = strsplit(Experiment, split = "_")

l.df = {}

for (i in 1:length(Experiment)){
	df.i = read.csv(paste("results/05_nonB_motif_cluster_analysis/log10_linear_models/", Experiment[i], ".csv", sep = ""))
	df.i$Experiment = Experiment[i]
	df.i$Species = l.Experiment[[i]][1]
	df.i$Bio_sample = l.Experiment[[i]][2]
	df.i$Technique = l.Experiment[[i]][3]
	l.df[[i]] = df.i
}

df = bind_rows(l.df)

write.csv(df, "results/05_nonB_motif_cluster_analysis/log10_linear_models.csv", row.names = F)
