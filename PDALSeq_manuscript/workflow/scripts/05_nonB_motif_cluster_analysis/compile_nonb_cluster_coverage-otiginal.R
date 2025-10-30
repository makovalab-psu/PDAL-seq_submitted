library(tidyverse)

read.m = function(x = "data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/200/Hsapien.bed"){
	df = data.frame(readr::read_delim(x, col_names = c("chr", "start", "end", "cluster.ID", "length.N", "MR", "TRI", "STR", "Z", "APR", "DR", "SLS", "GQ", "IR", "CRU")))
	return(df)
}

read.gc = function(x = "data/05_nonB_motif_cluster_analysis/10_nonB_cluster_matrix_GC_content/200/Hsapien.txt"){
	df = data.frame(readr::read_delim(x, col_names = c("GC.count"), delim = "\t"))
	return(df)
}

read.r = function(x = "data/05_nonB_motif_cluster_analysis/09_nonb_cluster_coverage/200/Hsapien_GM24385_PDALseq.bg"){
	df = data.frame(readr::read_delim(x, col_names = c("chr", "start", "end", "read.count"), delim = "\t"))
	return(df)
}

df.index = read.delim("PDALseq_samples.txt")

Experiment = unique(df.index$Experiment)

l.Experiment = strsplit(Experiment, split = "_")

l.df = {}


for (i in 1:length(l.Experiment)){
	print(i)
	print(Experiment[i])
	Species = l.Experiment[[i]][1]
	print(Species)
	df.m = read.m(paste("data/05_nonB_motif_cluster_analysis/08_sex_chromosome_blacklist_removed_sorted_nonB_cluster_matrix/200/", Species, ".bed", sep = ""))
	df.gc = read.gc(paste("data/05_nonB_motif_cluster_analysis/10_nonB_cluster_matrix_GC_content/200/", Species, ".txt", sep = ""))
	df.r = read.r(paste("data/05_nonB_motif_cluster_analysis/09_nonb_cluster_coverage/200/", Experiment[i], ".bg", sep = ""))
	df.m$GC.count = df.gc$GC.count
	df.m$GC.content = df.m$GC.count/df.m$length.N
	df.m$read.count = df.r$read.count
	df.m$Species = l.Experiment[[i]][1]
	df.m$Bio_sample = l.Experiment[[i]][2]
	df.m$Technique = l.Experiment[[i]][3]
	df.m$Experiment = Experiment[i]
	l.df[[i]] = df.m
}

df = bind_rows(l.df)

write.csv(df, "results/05_nonB_motif_cluster_analysis/nonb_cluster_coverage.csv", row.names = F)
