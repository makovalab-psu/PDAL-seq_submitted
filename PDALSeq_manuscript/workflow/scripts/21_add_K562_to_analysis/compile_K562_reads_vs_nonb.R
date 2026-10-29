library(tidyverse)

df.index = read.delim("PDALseq_samples.txt")

K562_experiments = c("Hsapien_K562_PDALseq")

Experiment = c(unique(df.index$Experiment), K562_experiments)

nonB = c("all", "APR", "CRU", "DR", "GQ", "IR", "MR", "SLS", "STR", "TRI", "Z")
windows = c("1000000", "100000", "10000")

df.corr = data.frame(expand_grid(Experiment, nonB, windows))

l.Experiment = strsplit(df.corr$Experiment, split = "_")

get_cov_path = function(experiment, window){
	if (experiment %in% K562_experiments){
		paste("data/21_add_K562_to_analysis/04_nonB/04_genome_window_coverage/", window, "_nucleotides/", experiment, ".bg", sep = "")
	} else {
		paste("data/04_Genome_wide_comparison_to_nonB_motifs/04_genome_window_coverage/", window, "_nucleotides/", experiment, ".bg", sep = "")
	}
}

l.df.compare = {}
Species = c()
Bio_sample = c()
Technique = c()
nonB = c()
Window = c()
Spearman = c()
Pearson = c()

for (i in 1:nrow(df.corr)){
	print(i)
	Species[i] = l.Experiment[[i]][1]
	Bio_sample[i] = l.Experiment[[i]][2]
	Technique[i] = l.Experiment[[i]][3]
	nonB[i] = df.corr$nonB[i]
	Window[i] = df.corr$windows[i]
	df.nonb = read.delim(paste("data/04_Genome_wide_comparison_to_nonB_motifs/03_nonb_density/", df.corr$windows[i], "_nucleotides/", Species[i], "/", df.corr$nonB[i], ".bg", sep = ""), header = F)
	colnames(df.nonb) = c("chr", "start", "end", "count", "bases", "length", "base.density")
	df.reads = read.delim(get_cov_path(df.corr$Experiment[i], df.corr$windows[i]), header = F)
	colnames(df.reads) = c("chr", "start", "end", "reads")

	df.reads$nonb.count = df.nonb$count
	df.reads$Species = Species[i]
	df.reads$Bio_sample = Bio_sample[i]
	df.reads$Technique = Technique[i]
	df.reads$nonB = nonB[i]
	df.reads$Window = Window[i]

	l.df.compare[[i]] = df.reads

	df.nonb = df.nonb %>% filter(!chr %in% c("chrX", "chrY"))
	df.reads = df.reads %>% filter(!chr %in% c("chrX", "chrY"))
	t.pearson = cor.test(as.numeric(df.nonb$count), as.numeric(df.reads$reads), method = "pearson")
	Pearson[i] = t.pearson$estimate
	t.spearman = cor.test(as.numeric(df.nonb$count), as.numeric(df.reads$reads), method = "spearman")
	Spearman[i] = t.spearman$estimate
}

df.corr$Species = Species
df.corr$Bio_sample = Bio_sample
df.corr$Technique = Technique
df.corr$nonB = nonB
df.corr$Window = Window
df.corr$Spearman = Spearman
df.corr$Pearson = Pearson

df.compare = bind_rows(l.df.compare)

dir.create("results/21_add_K562_to_analysis", recursive = TRUE, showWarnings = FALSE)
write.csv(df.corr, "results/21_add_K562_to_analysis/nonB_correlation_coefficient.csv", row.names = F)
write.csv(df.compare %>% filter(Window == "1000000"), "results/21_add_K562_to_analysis/nonB_coverage_versus_read.csv", row.names = F)
