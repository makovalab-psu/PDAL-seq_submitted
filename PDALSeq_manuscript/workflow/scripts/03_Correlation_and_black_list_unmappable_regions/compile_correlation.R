library(tidyverse)

Windows = c("1000000", "100000", "10000")

read.bg = function(x){
	df = read.delim(x, sep = "", header = F)
	colnames(df) = c("chr", "start", "end", "reads")
	return(df)
}

df.index = read.delim("PDALseq_samples.txt")

l.df.index = df.index %>% group_by(Species) %>% group_split()

l.df.species = {}

for (i in 1:length(l.df.index)){
	df.i = expand.grid(l.df.index[[i]]$Sequencing_sample, l.df.index[[i]]$Sequencing_sample)
	colnames(df.i) = c("Sample_1", "Sample_2")
	df.i$Species = l.df.index[[i]]$Species[1]
	l.df.species[[i]] = df.i
}

df.species = bind_rows(l.df.species)

l.df.window = {}

for (i in 1:length(Windows)){
	df.i = df.species
	df.i$Window = Windows[i]
	l.df.window[[i]] = df.i
}

df = bind_rows(l.df.window)

Pearson = rep(NA, nrow(df))
Spearman = rep(NA, nrow(df))

for (i in 1:nrow(df)){
        print(i)
	Path1 = paste("data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/", df$Window[i], "_nucleotides/", df$Sample_1[i], ".bed", sep = "")
	Path2 = paste("data/03_Correlation_and_black_list_unmappable_regions/03_genome_window_coverage/", df$Window[i], "_nucleotides/", df$Sample_2[i], ".bed", sep = "")
	df1 = read.bg(Path1)
	df2 = read.bg(Path2)
	t.pearson = cor.test(df1$reads, df2$reads, method = "pearson")
	Pearson[i] = t.pearson$estimate
	t.spearman = cor.test(df1$reads, df2$reads, method = "spearman")
	Spearman[i] = t.spearman$estimate
}


df$Pearson = Pearson
df$Spearman = Spearman

write.csv(df, "results/03_Correlation_and_black_list_unmappable_regions/Technical_replicate_correlation.csv", row.names = F)

